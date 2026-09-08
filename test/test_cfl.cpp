/*--------------------------------------------------------------------------*/
/*---------------------------- File test_cfl.cpp ---------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Extirpation comparison: a real Capacitated Facility Location instance is
 * solved with BendersDecompositionSolver + BundleSolver on a solver-agnostic
 * 2-level AbstractBlock built from the instance data, and the result is
 * checked against the LP optimum of the same model (a monolithic solve). The
 * objectives must coincide (same model), and the solve times are reported.
 *
 * The instance is loaded through CapacitatedFacilityLocationBlock only to read
 * its data; the Benders model is then built generically (the problem Block
 * never builds a BendersBFunction).
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

#include <boost/multi_array.hpp>

#include "AbstractBlock.h"
#include "BendersDecompositionSolver.h"
#include "BlockSolverConfig.h"
#include "CapacitatedFacilityLocationBlock.h"
#include "ColVariable.h"
#include "FRealObjective.h"
#include "FRowConstraint.h"
#include "LinearFunction.h"
#include "Objective.h"

using namespace SMSpp_di_unipi_it;

using array_type = boost::multi_array< ColVariable , 2 >;
using CFLB = CapacitatedFacilityLocationBlock;

/*--------------------------------------------------------------------------*/
/*--------------------------- BLOCK BUILDERS -------------------------------*/
/*--------------------------------------------------------------------------*/

// big-M cost of an unserved-demand slack, large enough never to be used at the
// optimum, making the transportation subproblem always feasible
/* It has to dominate the cost of serving one customer, whichever facility
 * serves it, and nothing more: the *sum* of all the costs would do as well in
 * theory, but it grows with the instance, and a big-M orders of magnitude
 * larger than the optimum makes every relative tolerance meaningless, the
 * solvers declaring convergence on a value that is all slack. */

static double big_M( CFLB * B )
{
 double m = 0;
 for( CFLB::Index i = 0 ; i < B->get_NFacilities() ; ++i ) {
  m = std::max( m , B->get_Fixed_Cost( i ) );
  for( CFLB::Index j = 0 ; j < B->get_NCustomers() ; ++j )
   m = std::max( m , B->get_Transportation_Cost( i , j ) );
  }

 return( 10 * ( m + 1 ) );
 }

/*--------------------------------------------------------------------------*/

// add to block the transportation part of the CFL model ( the flow Variable x
// in [0,1], a slack, the demand and capacity Constraint -- the latter coupling
// the design Variable y -- ) and its costs to obj

static void add_transport( AbstractBlock * block , CFLB * B ,
			   std::vector< ColVariable > * y , double M ,
			   LinearFunction * obj )
{
 const int m = B->get_NFacilities() , n = B->get_NCustomers();

 boost::array< array_type::index , 2 > shape = { m , n };
 auto x = new array_type( shape );
 for( auto p = x->data() ; p < x->data() + x->num_elements() ; ++p )
  p->is_positive( true );
 block->add_static_variable( * x , "x" );
 auto sl = new std::vector< ColVariable >( n );
 for( auto & v : * sl ) v.is_positive( true );
 block->add_static_variable( * sl , "s" );

 auto dem = new std::vector< FRowConstraint >( n );
 for( int j = 0 ; j < n ; ++j ) {
  auto f = new LinearFunction();
  for( int i = 0 ; i < m ; ++i )
   f->add_variable( & ( * x )[ i ][ j ] , 1 );
  f->add_variable( & ( * sl )[ j ] , 1 );
  ( * dem )[ j ].set_function( f );
  ( * dem )[ j ].set_both( 1 );
  }
 block->add_static_constraint( * dem , "demand" );

 auto cap = new std::vector< FRowConstraint >( m );
 for( int i = 0 ; i < m ; ++i ) {
  auto f = new LinearFunction();
  for( int j = 0 ; j < n ; ++j )
   f->add_variable( & ( * x )[ i ][ j ] , B->get_Demand( j ) );
  f->add_variable( & ( * y )[ i ] , - B->get_Capacity( i ) );
  ( * cap )[ i ].set_function( f );
  ( * cap )[ i ].set_rhs( 0 );
  }
 block->add_static_constraint( * cap , "capacity" );

 for( int i = 0 ; i < m ; ++i )
  for( int j = 0 ; j < n ; ++j )
   obj->add_variable( & ( * x )[ i ][ j ] , B->get_Transportation_Cost( i , j ) );
 for( int j = 0 ; j < n ; ++j )
  obj->add_variable( & ( * sl )[ j ] , M );
 }

/*--------------------------------------------------------------------------*/

static std::vector< ColVariable > * make_y( CFLB * B , bool set_start )
{
 auto y = new std::vector< ColVariable >( B->get_NFacilities() );
 for( auto & y_i : * y ) {
  y_i.is_unitary( true );
  y_i.is_positive( true );
  if( set_start )
   y_i.set_value( 0.5 );
  }
 return( y );
 }

/*--------------------------------------------------------------------------*/

// the monolithic LP model ( y and x in one Block ): the reference optimum

static AbstractBlock * build_monolithic( CFLB * B , double M )
{
 const int m = B->get_NFacilities();
 auto block = new AbstractBlock();
 auto y = make_y( B , false );
 block->add_static_variable( * y , "y" );

 auto f = new LinearFunction();
 for( int i = 0 ; i < m ; ++i )
  f->add_variable( & ( * y )[ i ] , B->get_Fixed_Cost( i ) );
 add_transport( block , B , y , M , f );

 auto obj = new FRealObjective( block , f );
 obj->set_sense( Objective::eMin );
 block->set_objective( obj );
 return( block );
 }

/*--------------------------------------------------------------------------*/

// the solver-agnostic 2-level model: master Block with y, one nested LP
// transportation sub-Block whose capacity Constraint couple y

static AbstractBlock * build_structured( CFLB * B , double M )
{
 const int m = B->get_NFacilities();
 auto root = new AbstractBlock();
 auto y = make_y( B , true );
 root->add_static_variable( * y , "y" );

 auto df = new LinearFunction();
 for( int i = 0 ; i < m ; ++i )
  df->add_variable( & ( * y )[ i ] , B->get_Fixed_Cost( i ) );
 auto robj = new FRealObjective( root , df );
 robj->set_sense( Objective::eMin );
 root->set_objective( robj );

 auto sub = new AbstractBlock( root );
 auto sf = new LinearFunction();
 add_transport( sub , B , y , M , sf );
 auto sobj = new FRealObjective( sub , sf );
 sobj->set_sense( Objective::eMin );
 sub->set_objective( sobj );
 root->add_nested_Block( sub );

 return( root );
 }

/*--------------------------------------------------------------------------*/

static double solve( AbstractBlock * block , const std::string & cfg ,
		     double & seconds , int * status = nullptr ,
		     double * ub = nullptr , long * iters = nullptr ,
		     long * cuts = nullptr , int rounds = 0 )
{
 auto c = Configuration::deserialize( cfg );
 auto bsc = dynamic_cast< BlockSolverConfig * >( c );
 if( ! bsc ) { std::cerr << cfg << " not a BlockSolverConfig\n"; std::exit( 1 ); }
 bsc->apply( block );
 auto solver = block->get_registered_solvers().front();

 /* Capping the rounds turns the run into a fixed budget of cuts, which is
  * how the strength of a family of cuts is measured: whoever has the better
  * bound after the same number of rounds has the stronger cuts. */

 if( rounds > 0 )
  solver->set_par( solver->int_par_str2idx( "int_BDSlv_MaxRounds" ) , rounds );
 auto t0 = std::chrono::steady_clock::now();
 const int st = solver->compute( false );
 auto t1 = std::chrono::steady_clock::now();
 seconds = std::chrono::duration< double >( t1 - t0 ).count();

 /* The bound the Solver reports is only the optimum if it says it converged:
  * with the two of them, and the status, a value that is off tells whether
  * the Solver stopped early or believes what it says. */

 const double lb = solver->get_lb();
 if( status ) *status = st;
 if( ub ) *ub = solver->get_ub();
 if( iters ) *iters = solver->get_elapsed_iterations();
 if( cuts )
  if( auto bds = dynamic_cast< BendersDecompositionSolver * >( solver ) )
   *cuts = bds->get_num_cuts();

 /* Reading the bound is all that was needed: the Solver is un-registered
  * and deleted by applying the cleared BlockSolverConfig, which is what
  * gives the Block back whatever the Solver had taken from it. */

 bsc->clear();
 bsc->apply( block );
 delete bsc;
 return( lb );
 }

/*--------------------------------------------------------------------------*/

// solve the instance with CapacitatedFacilityLocationBlock's own ad-hoc Benders
// decomposition ( BenForm: a hidden BendersBFunction over an inner MCFBlock,
// hand-driven master and cut separation )

static double solve_benform( CFLB * B , double & seconds )
{
 auto bc = Configuration::deserialize( "BenForm_BCfg.txt" );
 dynamic_cast< BlockConfig * >( bc )->apply( B );
 B->generate_abstract_variables();
 auto sc = Configuration::deserialize( "BenForm_BSCfg.txt" );
 auto bsc = dynamic_cast< BlockSolverConfig * >( sc );
 bsc->apply( B );
 auto solver = B->get_registered_solvers().front();
 auto t0 = std::chrono::steady_clock::now();
 solver->compute( false );
 auto t1 = std::chrono::steady_clock::now();
 seconds = std::chrono::duration< double >( t1 - t0 ).count();
 const double lb = solver->get_lb();
 bsc->clear();
 delete bc;
 delete sc;
 return( lb );
 }

/*--------------------------------------------------------------------------*/
/*--------------------------------- MAIN -----------------------------------*/
/*--------------------------------------------------------------------------*/

int main( int argc , char ** argv )
{
 // link anchor for the factory self-registration ( see test.cpp )
 delete new BendersDecompositionSolver();

 const std::string fn = ( argc > 1 ) ? argv[ 1 ]
   : "../../CapacitatedFacilityLocationBlock/data/txt/ORLib/cap102.txt";

 auto B = new CFLB();
 B->Block::load( fn , 'C' );
 const double M = big_M( B );
 std::cout << "instance " << fn << ": " << B->get_NFacilities()
           << " facilities, " << B->get_NCustomers() << " customers"
           << std::endl;

 // ----- reference: monolithic LP ----------------------------------------- #
 auto mono = build_monolithic( B , M );
 double t_ref;
 const double ref = solve( mono , "BSPar_sub.txt" , t_ref );

 /* The ad hoc Benders of CapacitatedFacilityLocationBlock is the reference
  * this Solver is measured against, but on a large instance it takes tens of
  * minutes: any second argument skips it, so that the Solver alone can be
  * looked at without paying for it. */
 const bool with_benform = ( argc <= 2 );

 auto Bben = new CFLB();
 double t_ben = 0 , ben = 0;
 if( with_benform ) {
  Bben->Block::load( fn , 'C' );
  ben = solve_benform( Bben , t_ben );
  }

 // ----- BDS + BundleSolver on the solver-agnostic 2-level model ----------- #
 auto root = build_structured( B , M );
 double t_bds , ub_bds;
 int st_bds;
 const double bds = solve( root , "BSPar_benders_convex.txt" , t_bds ,
                           & st_bds , & ub_bds );

 /* ----- the MILP master, with and without the Pareto-optimal cuts -------- #
  *
  * Facility location is where the Pareto-optimal cuts were invented: the
  * transportation subproblem is massively dual degenerate, so which of its
  * optimal duals the Solver returns, and hence how strong the cut is, is
  * arbitrary. Both runs have to end at the same value; the figures to compare
  * are the rounds, the cuts and the time. */

 /* How many rounds the two capped runs get: enough to tell the cuts apart,
  * few enough that a pure cutting plane does not have to converge. */

 const int budget = ( argc > 3 ) ? std::stoi( argv[ 3 ] ) : 30;
 auto root_m = build_structured( B , M );
 auto root_p = build_structured( B , M );
 double t_m , t_p;
 int st_m , st_p;
 long it_m = 0 , it_p = 0 , ct_m = 0 , ct_p = 0;
 const double v_m = solve( root_m , "BSPar_benders_milp.txt" , t_m , & st_m ,
                           nullptr , & it_m , & ct_m , budget );
 const double v_p = solve( root_p , "BSPar_benders_milp_pareto.txt" , t_p ,
                           & st_p , nullptr , & it_p , & ct_p , budget );

 // ----- compare ---------------------------------------------------------- #
 const double tol = 1e-5;
 auto rel = [ ref ]( double v ) {
  return( std::abs( ref - v )
	  / std::max( 1.0 , std::max( std::abs( ref ) , std::abs( v ) ) ) );
  };
 const double e_ben = with_benform ? rel( ben ) : 0 , e_bds = rel( bds );

 /* The two capped runs are lower bounds, not optima: what they have to
  * satisfy is that they are below the optimum, and the interest is in which
  * of the two is closer to it. */

 const bool ok = ( e_ben <= tol ) && ( e_bds <= tol ) &&
                 ( v_m <= ref * ( 1 + tol ) ) && ( v_p <= ref * ( 1 + tol ) );
 std::cout.precision( 10 );
 std::cout << "monolithic LP   = " << ref << "  ( " << t_ref << " s )\n";
 if( with_benform )
  std::cout << "CFLB BenForm    = " << ben << "  ( " << t_ben << " s )  err "
            << e_ben << "\n";
 std::cout << "BDS+Bundle      = " << bds << "  ( " << t_bds << " s )  err "
           << e_bds << " , status " << st_bds << " , ub " << ub_bds << "\n";
 std::cout << "BDS+MILP master = " << v_m << "  ( " << t_m << " s )  gap "
           << rel( v_m ) << " , " << it_m << " rounds , " << ct_m << " cuts\n";
 std::cout << "BDS+MILP+Pareto = " << v_p << "  ( " << t_p << " s )  gap "
           << rel( v_p ) << " , " << it_p << " rounds , " << ct_p << " cuts\n";
 std::cout << ( ok ? "-> OK ( same model )" : "-> FAIL" ) << std::endl;

 delete root_m;
 delete root_p;
 delete root;
 delete Bben;
 delete mono;
 delete B;
 return( ok ? 0 : 1 );
 }

/*--------------------------------------------------------------------------*/
/*------------------------- End File test_cfl.cpp --------------------------*/
/*--------------------------------------------------------------------------*/
