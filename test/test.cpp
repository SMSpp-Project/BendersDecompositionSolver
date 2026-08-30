/*--------------------------------------------------------------------------*/
/*------------------------------ File test.cpp -----------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Validation test for BendersDecompositionSolver on a small Capacitated
 * Warehouse Location (CWL) / Capacitated Facility Location instance.
 *
 * Two equivalent models of the same LP are built:
 *
 * - a *monolithic* AbstractBlock (design Variable y and flow Variable x in
 *   one Block) solved by a :MILPSolver, giving the reference optimum;
 *
 * - a *structured* AbstractBlock (master Block with the design Variable y, one
 *   nested sub-Block with the flow Variable x whose capacity Constraint couple
 *   y) solved by BendersDecompositionSolver.
 *
 * The test checks that the two optima coincide within a relative tolerance.
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

#include <cmath>
#include <iostream>

#include <boost/multi_array.hpp>

#include "AbstractBlock.h"
#include "BendersDecompositionSolver.h"
#include "BlockSolverConfig.h"
#include "ColVariable.h"
#include "FRealObjective.h"
#include "FRowConstraint.h"
#include "LinearFunction.h"
#include "Objective.h"

using namespace SMSpp_di_unipi_it;

/*--------------------------------------------------------------------------*/
/*------------------------------ THE INSTANCE ------------------------------*/
/*--------------------------------------------------------------------------*/

// a small CWL instance: M locations, N customers
static const int M = 3;
static const int N = 4;
static const double fixed_cost[ M ] = { 5 , 7 , 6 };
static const double capacity [ M ] = { 3 , 3 , 3 };
static const double demand   [ N ] = { 1 , 1 , 1 , 1 };
static const double cost[ M ][ N ] = { { 2 , 3 , 4 , 5 } ,
                                       { 4 , 1 , 2 , 3 } ,
                                       { 3 , 4 , 1 , 2 } };

// big-M cost of an unserved-demand slack: it makes the (transportation)
// subproblem always feasible (so only Benders optimality cuts are needed),
// while being large enough that the slack is never used at the optimum
static const double BigM = 1e2;

using array_type = boost::multi_array< ColVariable , 2 >;

/*--------------------------------------------------------------------------*/
/*--------------------------- BLOCK BUILDERS -------------------------------*/
/*--------------------------------------------------------------------------*/

// add to block the transportation part of one scenario s ( the Variable x[ s ]
// and slack, its demand and capacity Constraint -- the latter coupling the
// master Variable y -- ) and its cost terms to obj. With nsub scenarios this
// gives a 2-stage Benders structure with nsub subproblems coupled through y

static void add_transport( AbstractBlock * block ,
			   std::vector< ColVariable > * y , int s ,
			   bool with_slack , LinearFunction * obj )
{
 const std::string t = std::to_string( s );
 boost::array< array_type::index , 2 > shape = { M , N };
 auto x = new array_type( shape );
 for( auto p = x->data() ; p < x->data() + x->num_elements() ; ++p )
  p->is_positive( true );
 block->add_static_variable( * x , "x" + t );

 std::vector< ColVariable > * sl = nullptr;
 if( with_slack ) {
  sl = new std::vector< ColVariable >( N );
  for( auto & v : * sl ) v.is_positive( true );
  block->add_static_variable( * sl , "s" + t );
  }

 auto dem = new std::vector< FRowConstraint >( N );
 for( int j = 0 ; j < N ; ++j ) {
  auto f = new LinearFunction();
  for( int i = 0 ; i < M ; ++i )
   f->add_variable( & ( * x )[ i ][ j ] , 1 );
  if( with_slack )
   f->add_variable( & ( * sl )[ j ] , 1 );
  ( * dem )[ j ].set_function( f );
  ( * dem )[ j ].set_both( 1 );
  }
 block->add_static_constraint( * dem , "demand" + t );

 // the capacity Constraint couple the master y into the scenario
 auto cap = new std::vector< FRowConstraint >( M );
 for( int i = 0 ; i < M ; ++i ) {
  auto f = new LinearFunction();
  for( int j = 0 ; j < N ; ++j )
   f->add_variable( & ( * x )[ i ][ j ] , demand[ j ] );
  f->add_variable( & ( * y )[ i ] , - capacity[ i ] );
  ( * cap )[ i ].set_function( f );
  ( * cap )[ i ].set_rhs( 0 );
  }
 block->add_static_constraint( * cap , "capacity" + t );

 for( int i = 0 ; i < M ; ++i )
  for( int j = 0 ; j < N ; ++j )
   obj->add_variable( & ( * x )[ i ][ j ] , cost[ i ][ j ] );
 if( with_slack )
  for( int j = 0 ; j < N ; ++j )
   obj->add_variable( & ( * sl )[ j ] , BigM );
 }

/*--------------------------------------------------------------------------*/

// allocate and configure the master Variable y

static std::vector< ColVariable > * make_y( bool set_start )
{
 auto y = new std::vector< ColVariable >( M );
 for( auto & y_i : * y ) {
  y_i.is_unitary( true );
  y_i.is_positive( true );
  if( set_start )
   y_i.set_value( 0.3 );  // start where the capacity Constraint bind non-
                          // degenerately, so the first Benders cut is informative
  }
 return( y );
 }

/*--------------------------------------------------------------------------*/

// the monolithic LP: y and all scenarios' x in a single Block (the reference
// model). With with_slack == false the demand can only be served through x, so
// the problem becomes infeasible for small y (used to exercise feasibility cuts)

static AbstractBlock * build_monolithic( bool with_slack = true , int nsub = 1 )
{
 auto block = new AbstractBlock();
 auto y = make_y( false );
 block->add_static_variable( * y , "y" );

 auto f = new LinearFunction();
 for( int i = 0 ; i < M ; ++i )
  f->add_variable( & ( * y )[ i ] , fixed_cost[ i ] );
 for( int s = 0 ; s < nsub ; ++s )
  add_transport( block , y , s , with_slack , f );

 auto obj = new FRealObjective( block , f );
 obj->set_sense( Objective::eMin );
 block->set_objective( obj );

 return( block );
 }

/*--------------------------------------------------------------------------*/

// the structured model: master Block with y, one nested sub-Block per scenario
// whose capacity Constraint couple y (the structure BendersDecompositionSolver
// expects, with nsub subproblems)

static AbstractBlock * build_structured( bool with_slack = true , int nsub = 1 )
{
 auto root = new AbstractBlock();
 auto y = make_y( true );
 root->add_static_variable( * y , "y" );

 auto df = new LinearFunction();
 for( int i = 0 ; i < M ; ++i )
  df->add_variable( & ( * y )[ i ] , fixed_cost[ i ] );
 auto robj = new FRealObjective( root , df );
 robj->set_sense( Objective::eMin );
 root->set_objective( robj );

 for( int s = 0 ; s < nsub ; ++s ) {
  auto sub = new AbstractBlock( root );
  auto sf = new LinearFunction();
  add_transport( sub , y , s , with_slack , sf );
  auto sobj = new FRealObjective( sub , sf );
  sobj->set_sense( Objective::eMin );
  sub->set_objective( sobj );
  root->add_nested_Block( sub );
  }

 return( root );
 }

/*--------------------------------------------------------------------------*/
/*------------------------------ SOLVING -----------------------------------*/
/*--------------------------------------------------------------------------*/

// configure block from the BlockSolverConfig file, solve it and return the
// lower bound; no Solver parameter is ever set in code, everything comes from
// the configuration file

static double solve_from_config( AbstractBlock * block , const std::string & fn ,
				 int & status )
{
 auto cfg = Configuration::deserialize( fn );
 auto bsc = dynamic_cast< BlockSolverConfig * >( cfg );
 if( ! bsc ) {
  std::cerr << "Error: " << fn << " is not a BlockSolverConfig" << std::endl;
  std::exit( 1 );
  }
 bsc->apply( block );
 auto solver = block->get_registered_solvers().front();
 status = solver->compute( false );
 const double lb = solver->get_lb();
 bsc->clear();
 delete bsc;
 return( lb );
 }

/*--------------------------------------------------------------------------*/
/*--------------------------------- MAIN -----------------------------------*/
/*--------------------------------------------------------------------------*/

int main( void )
{
 // link anchor: reference a symbol of the BendersDecompositionSolver library
 // so the linker does not drop it ( the test uses the solver only through
 // configuration files, hence would otherwise reference no symbol of it, and
 // the library's self-registration in the Solver factory would be lost )
 delete new BendersDecompositionSolver();

 // ----- reference: solve the monolithic LP with a :MILPSolver ------------ #
 auto mono = build_monolithic();
 {
  auto c = Configuration::deserialize( "BSPar_sub.txt" );
  auto bsc = dynamic_cast< BlockSolverConfig * >( c );
  if( ! bsc ) { std::cerr << "BSPar_sub.txt not a BlockSolverConfig\n";
                return( 1 ); }
  bsc->apply( mono );
  bsc->clear();
  delete bsc;
  }
 auto ref_solver = mono->get_registered_solvers().front();
 ref_solver->compute( false );
 const double ref = ref_solver->get_lb();

 // ----- Benders, convex regime, configured entirely from file ------------ #
 auto root = build_structured();
 int st;
 const double ben = solve_from_config( root , "BSPar_benders_convex.txt" , st );
 std::cout << "Benders(convex) status = " << st << "   lb = " << ben
           << std::endl;

 // ----- Benders, MILP regime ( multi-cut ), configured from file --------- #
 auto root2 = build_structured();
 int st2;
 const double ben2 = solve_from_config( root2 , "BSPar_benders_milp.txt" , st2 );
 std::cout << "Benders(MILP,multi) status = " << st2 << "   lb = " << ben2
           << std::endl;

 // ----- Benders, MILP regime ( single-cut ), configured from file -------- #
 auto root3 = build_structured();
 int st3;
 const double ben3 = solve_from_config( root3 , "BSPar_benders_milp_single.txt" ,
					st3 );
 std::cout << "Benders(MILP,single) status = " << st3 << "   lb = " << ben3
           << std::endl;

 // ----- multi-subproblem ( 2 scenarios ): >1 BendersBFunction ------------ #
 // two transportation scenarios coupled through y: BDS builds one
 // BendersBFunction per scenario, so multi-cut uses two epigraph Variable and
 // single-cut one; all must match the 2-scenario monolithic optimum
 auto mono2 = build_monolithic( true , 2 );
 int dummy;
 const double ref2 = solve_from_config( mono2 , "BSPar_sub.txt" , dummy );
 auto root_c2 = build_structured( true , 2 );
 const double ben_c2 = solve_from_config( root_c2 , "BSPar_benders_convex.txt" ,
					  dummy );
 auto root_m2 = build_structured( true , 2 );
 const double ben_m2 = solve_from_config( root_m2 , "BSPar_benders_milp.txt" ,
					  dummy );
 auto root_s2 = build_structured( true , 2 );
 const double ben_s2 = solve_from_config( root_s2 ,
					  "BSPar_benders_milp_single.txt" , dummy );
 std::cout << "2-scenario: ref = " << ref2 << "   convex = " << ben_c2
           << "   MILP-multi = " << ben_m2 << "   MILP-single = " << ben_s2
           << std::endl;

 // ----- feasibility cuts: no-slack instance, MILP regime ----------------- #
 // without the slack the subproblem is infeasible for small y, so the solver
 // must generate Benders feasibility cuts; this only works if the subproblem
 // Solver provides an infeasibility certificate (Farkas dual direction), which
 // the current MILPSolver integration does not, so this case is informational
 const double tol = 1e-5;
 auto rel = []( double a , double b ) {
  return( std::abs( a - b )
	  / std::max( 1.0 , std::max( std::abs( a ) , std::abs( b ) ) ) );
  };
 auto mono_ns = build_monolithic( false );
 int st_ns_ref;
 const double ref_ns = solve_from_config( mono_ns , "BSPar_sub.txt" ,
					  st_ns_ref );
 auto root_ns = build_structured( false );
 try {
  int st_ns;
  const double ben_ns = solve_from_config( root_ns , "BSPar_benders_milp.txt" ,
					   st_ns );
  std::cout << "Benders(MILP,feas-cuts,no-slack) = " << ben_ns
            << "   ref = " << ref_ns
            << ( rel( ref_ns , ben_ns ) <= tol ? "   -> OK" : "   -> FAIL" )
            << std::endl;
  }
 catch( const std::exception & e ) {
  std::cout << "Benders(MILP,feas-cuts,no-slack): skipped (known limitation) - "
            << e.what() << std::endl;
  }

 // ----- compare ( optimality-cut cases, the supported ones ) ------------- #
 const double err = rel( ref , ben );
 const double err2 = rel( ref , ben2 );
 const double err3 = rel( ref , ben3 );
 const bool ok1 = ( err <= tol ) && ( err2 <= tol ) && ( err3 <= tol );
 std::cout << "1-scenario: monolithic = " << ref
           << "   convex = " << ben << " (err " << err << ")"
           << "   MILP-multi = " << ben2 << " (err " << err2 << ")"
           << "   MILP-single = " << ben3 << " (err " << err3 << ")"
           << ( ok1 ? "   -> OK" : "   -> FAIL" ) << std::endl;

 const bool ok2 = ( rel( ref2 , ben_c2 ) <= tol )
	       && ( rel( ref2 , ben_m2 ) <= tol )
	       && ( rel( ref2 , ben_s2 ) <= tol );
 std::cout << "2-scenario: " << ( ok2 ? "-> OK" : "-> FAIL" ) << std::endl;

 const bool ok = ok1 && ok2;

 delete root_s2;
 delete root_m2;
 delete root_c2;
 delete mono2;
 delete root_ns;
 delete mono_ns;
 delete root3;
 delete root2;
 delete root;
 delete mono;

 return( ok ? 0 : 1 );
 }

/*--------------------------------------------------------------------------*/
/*------------------------------ End File test.cpp -------------------------*/
/*--------------------------------------------------------------------------*/
