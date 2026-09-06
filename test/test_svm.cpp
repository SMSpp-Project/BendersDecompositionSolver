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

 /* Reading the bound is all that was needed: the Solver is un-registered
  * and deleted by applying the cleared BlockSolverConfig, which is what
  * gives the Block back whatever the Solver had taken from it. */

 bsc->clear();
 bsc->apply( block );
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
 ben.generate_abstract_variables();
 ben.generate_abstract_constraints();
 ben.generate_objective();

 double t_bd;
 int st_bd;
 const double bd = solve_from_config( & ben , "BSPar_svm_benders.txt" , st_bd ,
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

 return( ok ? 0 : 1 );
 }

/*--------------------------------------------------------------------------*/
/*---------------------------- End File test_svm.cpp -----------------------*/
/*--------------------------------------------------------------------------*/
