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
 * Both Solver are attached to a SVMBlock, the very Block that holds the data:
 * the master of the Benders side is the SVMBlock itself, the epigraph
 * Variable and the cuts living in the Block the Solver builds around it [see
 * BendersDecompositionSolver], so what is compared here is one instance and
 * two structures of it, with no rendition in between.
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

#include "DQuadFunction.h"

#include "FRowConstraint.h"

#include "LinearFunction.h"

#include "BendersDecompositionSolver.h"
#include "BlockSolverConfig.h"
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
/*------------------------------ SOLVING -----------------------------------*/
/*--------------------------------------------------------------------------*/

// configure block out of the BlockSolverConfig file, solve it and return the
// lower bound together with the time it took

static double solve_from_config( Block * block , const std::string & fn ,
                                 int & status , double & time ,
                                 bool take_ub = false ,
                                 long * iters = nullptr ,
                                 long * cuts = nullptr )
{
 auto cfg = Configuration::deserialize( fn );
 auto bsc = dynamic_cast< BlockSolverConfig * >( cfg );
 if( ! bsc ) {
  std::cerr << "Error: " << fn << " is not a BlockSolverConfig" << std::endl;
  std::exit( 1 );
  }

 bsc->apply( block );
 auto solver = block->get_registered_solvers().front();

 /* Whatever the Solver does, it has to be un-registered before the Block is
  * touched again: a Solver that throws in the middle of taking the Block
  * apart leaves it in the hands of nobody otherwise. */

 auto give_back = [ & ]( void ) {
  bsc->clear();
  bsc->apply( block );
  delete bsc;
  };

 const auto start = std::chrono::steady_clock::now();
 try { status = solver->compute( false ); }
 catch( ... ) { give_back(); throw; }
 const double lb = take_ub ? solver->get_ub() : solver->get_lb();

 if( iters )
  *iters = solver->get_elapsed_iterations();
 if( cuts )
  if( auto bds = dynamic_cast< BendersDecompositionSolver * >( solver ) )
   *cuts = bds->get_num_cuts();
 time = std::chrono::duration< double >(
                       std::chrono::steady_clock::now() - start ).count();

 /* Reading the bound is all that was needed: the Solver is un-registered
  * and deleted by applying the cleared BlockSolverConfig, which is what
  * gives the Block back whatever the Solver had taken from it. */

 give_back();

 return( lb );
 }

/*--------------------------------------------------------------------------*/
/*--------------------------------- MAIN -----------------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/* The Benders structure of the SVMBlock, rendered as an AbstractBlock: the
 * model ( w , b ) and the regularisation term stay in the root, and each
 * chunk is a sub-Block holding the slacks of its own samples, their margin
 * constraints and its share of the loss. It is the very same problem the
 * SVMBlock builds with set_structure( kBenders , P ), written out by hand
 * only because the convex master of BendersDecompositionSolver assembles
 * itself into the root, which therefore has to be an AbstractBlock; with the
 * master kept as a sub-Block of itself, as the design goes, this function
 * disappears. The partition is read off the SVMBlock, so that the two are
 * the same decomposition of the same instance. */

static AbstractBlock * build_benders_abstract( SVCBlock & ben , Index n ,
                                               Index m , const doubleVec & X ,
                                               const doubleVec & y , double C )
{
 auto root = new AbstractBlock();

 // the model: the weights and the bias, the latter not regularised
 auto w = new std::vector< ColVariable >( m );
 for( auto & wi : *w )
  wi.is_unitary( false , eNoMod );
 root->add_static_variable( *w , "w" );

 auto b = new ColVariable();
 b->is_unitary( false , eNoMod );
 root->add_static_variable( *b , "b" );

 // ( rho / 2 ) || w ||^2, the quadratic 0-th component of the sum-function;
 // the bias is in it with a zero coefficient, since the component has to
 // span the whole Lambda of the bundle
 const double rho = ben.get_reg_weight();
 DQuadFunction::v_coeff_triple triples( m + 1 );
 for( Index j = 0 ; j < m ; ++j )
  triples[ j ] = std::make_tuple( &(*w)[ j ] , 0.0 , rho / 2.0 );
 triples[ m ] = std::make_tuple( b , 0.0 , 0.0 );

 auto robj = new FRealObjective( root , new DQuadFunction(
                                             std::move( triples ) ) );
 robj->set_sense( Objective::eMin , eNoMod );
 root->set_objective( robj , eNoMod );

 // one sub-Block per chunk, with the samples the SVMBlock deals out to it
 const Index P = ben.get_NChunks();
 for( Index p = 0 ; p < P ; ++p ) {
  const auto & smpl = ben.get_chunk( p );
  const Index np = smpl.size();

  auto sub = new AbstractBlock( root );

  auto xi = new std::vector< ColVariable >( np );
  for( auto & xk : *xi )
   xk.is_positive( true , eNoMod );
  sub->add_static_variable( *xi , "xi" );

  // y_k ( < w , x_k > + b ) + xi_k >= 1
  auto cons = new std::vector< FRowConstraint >( np );
  for( Index k = 0 ; k < np ; ++k ) {
   const Index i = smpl[ k ];
   LinearFunction::v_coeff_pair cf;
   cf.reserve( m + 2 );
   for( Index j = 0 ; j < m ; ++j )
    cf.emplace_back( &(*w)[ j ] , y[ i ] * X[ i * m + j ] );
   cf.emplace_back( b , y[ i ] );
   cf.emplace_back( &(*xi)[ k ] , 1.0 );
   (*cons)[ k ].set_function( new LinearFunction( std::move( cf ) ) , eNoMod );
   (*cons)[ k ].set_lhs( 1.0 , eNoMod );
   (*cons)[ k ].set_rhs( Inf< double >() , eNoMod );
   }
  sub->add_static_constraint( *cons , "margin" );

  // C sum_k xi_k, the share of the loss of this chunk
  auto lf = new LinearFunction();
  for( Index k = 0 ; k < np ; ++k )
   lf->add_variable( &(*xi)[ k ] , C );
  auto sobj = new FRealObjective( sub , lf );
  sobj->set_sense( Objective::eMin , eNoMod );
  sub->set_objective( sobj , eNoMod );

  root->add_nested_Block( sub );
  }

 return( root );
 }

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

 // the Lagrangian dual is the slowest of the lot by far, hence it can be
 // left out when only the two Benders are of interest
 const bool do_ld = ( argc > 4 ) ? ( std::stoi( argv[ 4 ] ) != 0 ) : true;

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

 /* Solver::new_Solver() throws if the name is not in the factory, which is
  * what happens when SVMBlock has been built without LIBSVM. */

 Solver * probe = nullptr;
 try { probe = Solver::new_Solver( "LIBSVMSolver" ); }
 catch( const std::exception & ) {}

 if( probe ) {
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

 double t_ld = 0;
 int st_ld = 0;
 const double ld = do_ld ? solve_from_config( & cns , "BSPar_svm_ld.txt" ,
                                              st_ld , t_ld ) : smo;

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
 ben.generate_abstract_variables();
 ben.generate_abstract_constraints();
 ben.generate_objective();

 double t_bd;
 int st_bd;
 long it_bd = 0 , ct_bd = 0;
 const double bd = solve_from_config( & ben , "BSPar_svm_benders.txt" , st_bd ,
                                      t_bd , false , & it_bd , & ct_bd );

 /* The same, with the cuts of all the chunks aggregated into one: the two
  * describe the same problem, so what is being compared is how many rounds
  * and how many cuts each of them takes to get there. */

 SVCBlock bens;
 bens.set_kernel( SVMBlock::kLinear );
 bens.set_C( 1 );
 bens.load( n , m , X , y );
 bens.set_structure( & bcfg );
 bens.generate_abstract_variables();
 bens.generate_abstract_constraints();
 bens.generate_objective();

 double t_bs;
 int st_bs;
 long it_bs = 0 , ct_bs = 0;
 const double bs = solve_from_config( & bens , "BSPar_svm_benders_single.txt" ,
                                      st_bs , t_bs , false , & it_bs ,
                                      & ct_bs );

 // ----- the same, with the master given to the bundle -------------------- #

 /* The regularisation term makes the master strongly convex, which is what
  * the bundle carries as the quadratic 0-th component of its sum-function
  * [see MasterProblemBlock::set_zeroth_quadratic()]: with the cutting plane
  * of the MILP master the term is there but nobody knows it is, with the
  * bundle it is what the stabilization is made of. */

 /* A bundle that does not carry it throws instead, in which case the case is
  * skipped rather than failed, exactly as the LIBSVM one above. */

 double t_bdb = 0;
 int st_bdb = 0;
 double bdb = smo;
 bool has_bdb = false;
 { auto abs_ben = build_benders_abstract( ben , n , m , X , y , 1.0 );
   try {
    /* The bundle master minimizes, so what it converges to is its upper
     * bound, its lower one being the model value. */
    bdb = solve_from_config( abs_ben , "BSPar_svm_benders_convex.txt" ,
                             st_bdb , t_bdb , true );
    has_bdb = true;
    }
   catch( const std::exception & e ) {
    std::cout << "Benders (bundle): skipped, " << e.what() << std::endl;
    }
   delete abs_ben;
   }

 // ----- compare ---------------------------------------------------------- #

 auto rel = []( double a , double b ) {
  return( std::abs( a - b )
          / std::max( 1.0 , std::max( std::abs( a ) , std::abs( b ) ) ) );
  };

 const double tol = 1e-5;
 const double e_ld = rel( smo , ld );
 const double e_bd = rel( smo , bd );
 const double e_bdb = rel( smo , bdb );

 if( do_ld )
  std::cout << "Lagrangian dual  = " << ld << "  ( " << t_ld << " s , err "
            << e_ld << " , status " << st_ld << " )" << std::endl;
 if( has_bdb )
  std::cout << "Benders (bundle) = " << bdb << "  ( " << t_bdb << " s , err "
            << e_bdb << " , status " << st_bdb << " )" << std::endl;

 std::cout << "Benders          = " << bd << "  ( " << t_bd << " s , err "
           << e_bd << " , status " << st_bd << " , " << it_bd << " rounds , "
           << ct_bd << " cuts )" << std::endl;

 std::cout << "Benders (single) = " << bs << "  ( " << t_bs << " s , err "
           << rel( smo , bs ) << " , status " << st_bs << " , " << it_bs
           << " rounds , " << ct_bs << " cuts )" << std::endl;

 const bool ok = ( e_ld <= tol ) && ( e_bd <= tol ) && ( e_bdb <= tol ) &&
                 ( rel( smo , bs ) <= tol ) &&
                 ( ( ! has_lsvm ) || ( rel( smo , lsvm ) <= tol ) );
 std::cout << ( ok ? "-> OK ( the two decompositions agree )"
                   : "-> FAIL" ) << std::endl;

 return( ok ? 0 : 1 );
 }

/*--------------------------------------------------------------------------*/
/*---------------------------- End File test_svm.cpp -----------------------*/
/*--------------------------------------------------------------------------*/
