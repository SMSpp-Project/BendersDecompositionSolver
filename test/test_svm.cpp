/*--------------------------------------------------------------------------*/
/*---------------------------- File test_svm.cpp ---------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * The two dual ways of splitting a SVM training problem along the samples,
 * one against the other on the very same instance.
 *
 * The training problem is the sum over the samples of a loss plus one
 * regularisation term, so dealing the samples out to P chunks splits it, and
 * it does so in two opposite ways [see SVMBlock::set_structure()]:
 *
 * - the *consensus* one, in which each chunk holds a whole SVM with its own
 *   copy of the model and its share of the regularisation term, the copies
 *   being tied by consensus Constraint: relaxing those is the Lagrangian, or
 *   equivalently Dantzig-Wolfe, decomposition, and it is what
 *   LagrangianDualSolver does;
 *
 * - the *Benders* one, in which the model and the regularisation term stay in
 *   the master and each chunk holds only the slacks of its samples and their
 *   loss: projecting the slacks out leaves the loss of the chunk as a value
 *   function of the model, and approximating it from below is what
 *   BendersDecompositionSolver does.
 *
 * Both are exact reformulations of the same problem, hence they must agree
 * with each other and with the ad hoc SMOSolver, which ignores the structure
 * altogether: that they do is the test, and how they get there is what the
 * comparison is about.
 *
 * Note that the Benders side is built here as an AbstractBlock rendition of
 * the structure the SVMBlock itself builds, rather than by attaching the
 * Solver to the SVMBlock: BendersDecompositionSolver needs to add the
 * epigraph Variable and the cuts to the master, which only an AbstractBlock
 * lets a Solver do. The rendition goes away as soon as the core lets a Solver
 * own the master and exclude the sub-Block it must not see.
 *
 * \author Antonio Frangioni \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \author Donato Meoli \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \copyright &copy; by Antonio Frangioni, Donato Meoli
 */
/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include <chrono>
#include <cmath>
#include <iostream>
#include <random>

#include "AbstractBlock.h"
#include "BendersDecompositionSolver.h"
#include "BlockSolverConfig.h"
#include "DQuadFunction.h"
#include "FRealObjective.h"
#include "FRowConstraint.h"
#include "LinearFunction.h"
#include "OneVarConstraint.h"
#include "SMOSolver.h"
#include "SVCBlock.h"

using namespace SMSpp_di_unipi_it;

using Index = Block::Index;
using Subset = Block::Subset;
using doubleVec = SVMBlock::doubleVec;

/*--------------------------------------------------------------------------*/
/*------------------------------ THE INSTANCE ------------------------------*/
/*--------------------------------------------------------------------------*/

// a linearly separable two-class data set, the same generator the SVMBlock
// tester uses

static void make_data( Index n , Index m , doubleVec & X , doubleVec & y ,
                       unsigned seed )
{
 std::mt19937 rng( seed );
 std::normal_distribution< double > gauss( 0 , 1 );

 X.resize( std::size_t( n ) * m );
 y.resize( n );

 for( Index i = 0 ; i < n ; ++i ) {
  const double lbl = ( i % 2 ) ? 1 : -1;
  y[ i ] = lbl;
  for( Index j = 0 ; j < m ; ++j )
   X[ std::size_t( i ) * m + j ] = gauss( rng ) + ( j ? 0 : 4 * lbl );
  }
 }

/*--------------------------------------------------------------------------*/
/*--------------------------- BLOCK BUILDERS -------------------------------*/
/*--------------------------------------------------------------------------*/

/// the AbstractBlock rendition of the Benders structure of \p svm
/** Builds the very structure SVMBlock::set_structure() builds with kBenders,
 * i.e., the model and the regularisation term in the master and one sub-Block
 * per chunk holding the slacks of its samples, their margin Constraint and
 * their loss, out of the data of \p svm and of the partition it deals out. */

static AbstractBlock * build_benders( const SVMBlock * svm , Index P )
{
 const Index m = svm->get_NFeatures();
 const Index N = svm->get_NDual();
 const double rw = svm->get_reg_weight() / 2;
 const double C = svm->get_C();

 auto & s = svm->get_dual_signs();
 auto & di = svm->get_dual_samples();
 auto & q = svm->get_dual_costs();

 // the chunk of each dual index: a sample gives all of its dual indices to
 // the chunk it belongs to
 Subset smap( svm->get_NSamples() , 0 );
 for( Index p = 0 ; p < P ; ++p )
  for( auto i : svm->get_chunk( p ) )
   smap[ i ] = p;

 // ----- the master: the model and the regularisation term ---------------- #

 auto root = new AbstractBlock();

 auto w = new std::vector< ColVariable >( m );
 root->add_static_variable( *w , "w" );

 auto b = new ColVariable();
 root->add_static_variable( *b , "b" );

 DQuadFunction::v_coeff_triple triples( m + 1 );
 for( Index j = 0 ; j < m ; ++j )
  triples[ j ] = std::make_tuple( &(*w)[ j ] , double( 0 ) , rw );
 triples[ m ] = std::make_tuple( b , double( 0 ) ,
                                 svm->get_reg_bias() ? rw : double( 0 ) );

 auto robj = new FRealObjective( root ,
                                 new DQuadFunction( std::move( triples ) ) );
 robj->set_sense( Objective::eMin );
 root->set_objective( robj );

 // ----- one sub-Block per chunk: the slacks and the loss ----------------- #

 for( Index p = 0 ; p < P ; ++p ) {
  Subset dk;
  for( Index k = 0 ; k < N ; ++k )
   if( smap[ di[ k ] ] == p )
    dk.push_back( k );

  auto sub = new AbstractBlock( root );

  auto xi = new std::vector< ColVariable >( dk.size() );
  for( auto & xk : *xi )
   xk.is_positive( true );
  sub->add_static_variable( *xi , "xi" );

  // s_k ( < w , x_{ i( k ) } > + b ) + xi_k >= r_k
  auto cons = new std::vector< FRowConstraint >( dk.size() );
  for( Index t = 0 ; t < dk.size() ; ++t ) {
   const Index k = dk[ t ];
   const double sk = s[ k ];
   const double * xk = svm->get_x( di[ k ] );

   LinearFunction::v_coeff_pair cp( m + 2 );
   for( Index j = 0 ; j < m ; ++j )
    cp[ j ] = std::make_pair( &(*w)[ j ] , sk * xk[ j ] );
   cp[ m ] = std::make_pair( b , sk );
   cp[ m + 1 ] = std::make_pair( &(*xi)[ t ] , double( 1 ) );

   (*cons)[ t ].set_lhs( - q[ k ] );
   (*cons)[ t ].set_rhs( Inf< RowConstraint::RHSValue >() );
   (*cons)[ t ].set_function( new LinearFunction( std::move( cp ) ) );
   }
  sub->add_static_constraint( *cons , "cons" );

  LinearFunction::v_coeff_pair lp( dk.size() );
  for( Index t = 0 ; t < dk.size() ; ++t )
   lp[ t ] = std::make_pair( &(*xi)[ t ] , C );

  auto sobj = new FRealObjective( sub ,
                                  new LinearFunction( std::move( lp ) ) );
  sobj->set_sense( Objective::eMin );
  sub->set_objective( sobj );

  root->add_nested_Block( sub );
  }

 return( root );
 }

/*--------------------------------------------------------------------------*/
/*------------------------------ SOLVING -----------------------------------*/
/*--------------------------------------------------------------------------*/

// configure block out of the BlockSolverConfig file, solve it and return the
// lower bound together with the time it took

static double solve_from_config( Block * block , const std::string & fn ,
                                 int & status , double & time )
{
 auto cfg = Configuration::deserialize( fn );
 auto bsc = dynamic_cast< BlockSolverConfig * >( cfg );
 if( ! bsc ) {
  std::cerr << "Error: " << fn << " is not a BlockSolverConfig" << std::endl;
  std::exit( 1 );
  }

 bsc->apply( block );
 auto solver = block->get_registered_solvers().front();

 const auto start = std::chrono::steady_clock::now();
 status = solver->compute( false );
 const double lb = solver->get_lb();
 time = std::chrono::duration< double >(
                       std::chrono::steady_clock::now() - start ).count();

 bsc->clear();
 delete bsc;

 return( lb );
 }

/*--------------------------------------------------------------------------*/
/*--------------------------------- MAIN -----------------------------------*/
/*--------------------------------------------------------------------------*/

int main( int argc , char ** argv )
{
 // link anchor: the Solver are used through configuration files only, hence
 // no symbol of their libraries would be referenced [see test.cpp]
 delete new BendersDecompositionSolver();
 delete new SMOSolver();

 const Index n = ( argc > 1 ) ? std::stoi( argv[ 1 ] ) : 200;
 const Index m = ( argc > 2 ) ? std::stoi( argv[ 2 ] ) : 5;
 const Index P = ( argc > 3 ) ? std::stoi( argv[ 3 ] ) : 4;

 doubleVec X , y;
 make_data( n , m , X , y , 1 );

 std::cout << n << " samples, " << m << " features, " << P << " chunks"
           << std::endl;

 // ----- the reference: the ad hoc Solver on the whole problem ------------ #

 SVCBlock svm;
 svm.set_kernel( SVMBlock::kLinear );
 svm.set_C( 1 );
 svm.load( n , m , X , y );

 double t_smo;
 int st_smo;
 const double smo = solve_from_config( & svm , "BSPar_svm_smo.txt" , st_smo ,
                                       t_smo );

 std::cout << "SMOSolver        = " << smo << "  ( " << t_smo << " s )"
           << std::endl;

 // ----- the other yardstick: LIBSVM, if SVMBlock was built with it ------- #

 double lsvm = smo , t_lsvm = 0;
 bool has_lsvm = false;

 if( auto probe = Solver::new_Solver( "LIBSVMSolver" ) ) {
  delete probe;
  has_lsvm = true;

  SVCBlock lsv;
  lsv.set_kernel( SVMBlock::kLinear );
  lsv.set_C( 1 );
  lsv.load( n , m , X , y );

  int st_lsvm;
  lsvm = solve_from_config( & lsv , "BSPar_svm_libsvm.txt" , st_lsvm ,
                            t_lsvm );

  std::cout << "LIBSVMSolver     = " << lsvm << "  ( " << t_lsvm << " s )"
            << std::endl;
  }

 // ----- the consensus structure under a Lagrangian Solver ---------------- #

 SVCBlock cns;
 cns.set_kernel( SVMBlock::kLinear );
 cns.set_C( 1 );
 cns.load( n , m , X , y );

 SimpleConfiguration< std::pair< int , int > > ccfg(
  std::make_pair( int( SVMBlock::kConsensus ) , int( P ) ) );
 cns.set_structure( & ccfg );
 cns.generate_abstract_variables();
 cns.generate_abstract_constraints();
 cns.generate_objective();

 double t_ld;
 int st_ld;
 const double ld = solve_from_config( & cns , "BSPar_svm_ld.txt" , st_ld ,
                                      t_ld );

 // ----- the Benders structure under BendersDecompositionSolver ----------- #

 // the partition is the one the SVMBlock deals out, so that the two
 // decompositions split the very same samples the very same way
 SVCBlock ben;
 ben.set_kernel( SVMBlock::kLinear );
 ben.set_C( 1 );
 ben.load( n , m , X , y );

 SimpleConfiguration< std::pair< int , int > > bcfg(
  std::make_pair( int( SVMBlock::kBenders ) , int( P ) ) );
 ben.set_structure( & bcfg );

 auto root = build_benders( & ben , P );

 double t_bd;
 int st_bd;
 const double bd = solve_from_config( root , "BSPar_svm_benders.txt" , st_bd ,
                                      t_bd );

 // ----- compare ---------------------------------------------------------- #

 auto rel = []( double a , double b ) {
  return( std::abs( a - b )
          / std::max( 1.0 , std::max( std::abs( a ) , std::abs( b ) ) ) );
  };

 const double tol = 1e-5;
 const double e_ld = rel( smo , ld );
 const double e_bd = rel( smo , bd );

 std::cout << "Lagrangian dual  = " << ld << "  ( " << t_ld << " s , err "
           << e_ld << " , status " << st_ld << " )" << std::endl;
 std::cout << "Benders          = " << bd << "  ( " << t_bd << " s , err "
           << e_bd << " , status " << st_bd << " )" << std::endl;

 const bool ok = ( e_ld <= tol ) && ( e_bd <= tol ) &&
                 ( ( ! has_lsvm ) || ( rel( smo , lsvm ) <= tol ) );
 std::cout << ( ok ? "-> OK ( the two decompositions agree )"
                   : "-> FAIL" ) << std::endl;

 delete root;

 return( ok ? 0 : 1 );
 }

/*--------------------------------------------------------------------------*/
/*---------------------------- End File test_svm.cpp -----------------------*/
/*--------------------------------------------------------------------------*/
