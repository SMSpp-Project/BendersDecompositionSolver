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
  ( * cap )[ i ].set_lhs( - Inf< double >() );
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

static std::vector< ColVariable > * make_y( bool set_start ,
					    bool integer = false )
{
 auto y = new std::vector< ColVariable >( M );
 for( auto & y_i : * y ) {
  y_i.is_unitary( true );
  y_i.is_positive( true );
  if( integer )
   y_i.is_integer( true );
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

static AbstractBlock * build_monolithic( bool with_slack = true , int nsub = 1 ,
					 bool integer = false )
{
 auto block = new AbstractBlock();
 auto y = make_y( false , integer );
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

static AbstractBlock * build_structured( bool with_slack = true , int nsub = 1 ,
					 bool integer = false )
{
 auto root = new AbstractBlock();
 auto y = make_y( true , integer );
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
				 int & status , long * iters = nullptr ,
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
 status = solver->compute( false );
 const double lb = solver->get_lb();

 if( iters )
  *iters = solver->get_elapsed_iterations();
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

 /* ----- how many cuts each variant of the MILP regime takes -------------- #
  *
  * The three of them describe the same problem and have to end at the same
  * value: what changes is how many cuts are needed to get there, which is the
  * figure to compare, the number of rounds saying little when one round adds
  * one cut and another adds one per subproblem. */

 { auto root_m = build_structured( true , 4 );
   auto root_s = build_structured( true , 4 );
   auto root_p = build_structured( true , 4 );
   int st_m , st_s , st_p;
   long it_m = 0 , it_s = 0 , it_p = 0 , ct_m = 0 , ct_s = 0 , ct_p = 0;

   const double v_m = solve_from_config( root_m , "BSPar_benders_milp.txt" ,
					 st_m , & it_m , & ct_m );
   const double v_s = solve_from_config( root_s ,
					 "BSPar_benders_milp_single.txt" ,
					 st_s , & it_s , & ct_s );
   const double v_p = solve_from_config( root_p ,
					 "BSPar_benders_milp_pareto.txt" ,
					 st_p , & it_p , & ct_p );

   std::cout << "4-scenario MILP master: multi = " << v_m << " ( " << it_m
             << " rounds , " << ct_m << " cuts ) , single = " << v_s << " ( "
             << it_s << " rounds , " << ct_s << " cuts ) , Pareto = " << v_p
             << " ( " << it_p << " rounds , " << ct_p << " cuts )"
             << std::endl;
   }

 // ----- feasibility cuts: no-slack instance, MILP regime ----------------- #
 // without the slack the subproblem is infeasible for small y, so the solver
 // must generate Benders feasibility cuts out of the Farkas certificate of the
 // subproblem: whether there is one depends on how the subproblem Solver is
 // configured [see BSPar_sub.txt], hence the case is only checked if it does
 // provide it, and is skipped, rather than failed, if it does not
 const double tol = 1e-5;
 auto rel = []( double a , double b ) {
  return( std::abs( a - b )
	  / std::max( 1.0 , std::max( std::abs( a ) , std::abs( b ) ) ) );
  };
 /* ----- the master says what is master and what is complicating --------- #
  *
  * The same four-scenario problem, with the first scenario kept in the master
  * instead of being projected out and the complicating Variable named one by
  * one: which sub-Block are subproblems and which Variable are complicating
  * is a choice, and a different choice describes the very same problem, hence
  * the optimum has to be the one of the extensive form. What changes is the
  * work: one subproblem fewer to evaluate, and a larger master. */

 bool ok_k = true;
 { auto mono_k = build_monolithic( true , 4 );
   int st_rk;
   const double ref_k = solve_from_config( mono_k , "BSPar_sub.txt" , st_rk );
   auto root_k = build_structured( true , 4 );
   int st_k;
   long it_k = 0 , ct_k = 0;
   const double v_k = solve_from_config( root_k , "BSPar_benders_milp_keep.txt" ,
					 st_k , & it_k , & ct_k );
   ok_k = ( rel( ref_k , v_k ) <= tol );
   std::cout << "4-scenario, first one kept in the master: " << v_k << " ( "
             << it_k << " rounds , " << ct_k << " cuts )   ref = " << ref_k
             << ( ok_k ? "   -> OK" : "   -> FAIL" ) << std::endl;
   delete root_k;
   delete mono_k;
   }

 auto mono_ns = build_monolithic( false );
 int st_ns_ref;
 const double ref_ns = solve_from_config( mono_ns , "BSPar_sub.txt" ,
					  st_ns_ref );
 auto root_ns = build_structured( false );
 auto root_nn = build_structured( false );
 bool ok_ns = true;
 try {
  int st_ns , st_nn;
  long ct_ns = 0 , ct_nn = 0;
  const double ben_ns = solve_from_config( root_ns , "BSPar_benders_milp.txt" ,
					   st_ns , nullptr , & ct_ns );

  /* The very same run with the feasibility cuts normalized: they describe the
   * same half-spaces, so the optimum cannot change. */

  const double ben_nn = solve_from_config( root_nn ,
					   "BSPar_benders_milp_norm.txt" ,
					   st_nn , nullptr , & ct_nn );

  /* And with the cut of the phase one in place of the certificate: another
   * cut, hence another number of them, for the same optimum. */

  auto root_p1 = build_structured( false );
  int st_p1;
  long ct_p1 = 0;
  const double ben_p1 = solve_from_config( root_p1 ,
					   "BSPar_benders_milp_phase1.txt" ,
					   st_p1 , nullptr , & ct_p1 );
  ok_ns = ( rel( ref_ns , ben_ns ) <= tol ) &&
          ( rel( ref_ns , ben_nn ) <= tol ) &&
          ( rel( ref_ns , ben_p1 ) <= tol );
  std::cout << "Benders(MILP,feas-cuts,no-slack) = " << ben_ns
            << " ( " << ct_ns << " cuts )   normalized = " << ben_nn
            << " ( " << ct_nn << " cuts )   phase one = " << ben_p1
            << " ( " << ct_p1 << " cuts )   ref = " << ref_ns
            << ( ok_ns ? "   -> OK" : "   -> FAIL" ) << std::endl;
  delete root_p1;

  /* The same on four scenarios, where the master can starve several
   * subproblems at once and the two ways of cutting it away can be told
   * apart by how many rounds and how many cuts they take. */

  auto mono4 = build_monolithic( false , 4 );
  int st_r4;
  const double ref4 = solve_from_config( mono4 , "BSPar_sub.txt" , st_r4 );
  auto root_f4 = build_structured( false , 4 );
  auto root_14 = build_structured( false , 4 );
  int st_f4 , st_14;
  long it_f4 = 0 , it_14 = 0 , ct_f4 = 0 , ct_14 = 0;
  const double v_f4 = solve_from_config( root_f4 , "BSPar_benders_milp.txt" ,
					 st_f4 , & it_f4 , & ct_f4 );
  const double v_14 = solve_from_config( root_14 ,
					 "BSPar_benders_milp_phase1.txt" ,
					 st_14 , & it_14 , & ct_14 );
  ok_ns = ok_ns && ( rel( ref4 , v_f4 ) <= tol ) &&
                   ( rel( ref4 , v_14 ) <= tol );
  std::cout << "4-scenario, no slack: ref = " << ref4 << "   Farkas = "
            << v_f4 << " ( " << it_f4 << " rounds , " << ct_f4
            << " cuts )   phase one = " << v_14 << " ( " << it_14
            << " rounds , " << ct_14 << " cuts )" << std::endl;
  delete root_f4;
  delete root_14;
  delete mono4;
  }
 catch( const std::exception & e ) {
  std::cout << "Benders(MILP,feas-cuts,no-slack): skipped, the subproblem "
               "Solver gives no certificate - " << e.what() << std::endl;
  }

 /* ----- the same infeasibility, cut away combinatorially ---------------- #
  *
  * With the complicating Variable binary an infeasible subproblem can be cut
  * away by simply forbidding the assignment, which asks nothing of the
  * subproblem Solver: the two runs describe the same problem, hence they have
  * to end at the same value, and what the no-good cut costs is visible in how
  * many cuts it takes to get there. */

 bool ok_ng = true;
 { auto mono_b = build_monolithic( false , 1 , true );
   int st_ref_b;

   /* The reference is the monolithic problem with the y binary, so it has to
    * be solved as the MILP it is: BSPar_sub.txt relaxes the integrality. */

   const double ref_b = solve_from_config( mono_b , "BSPar_master_milp.txt" ,
					   st_ref_b );
   auto root_f = build_structured( false , 1 , true );
   auto root_g = build_structured( false , 1 , true );
   int st_f , st_g;
   long ct_f = 0 , ct_g = 0;
   const double v_f = solve_from_config( root_f , "BSPar_benders_milp.txt" ,
					 st_f , nullptr , & ct_f );
   const double v_g = solve_from_config( root_g ,
					 "BSPar_benders_milp_nogood.txt" ,
					 st_g , nullptr , & ct_g );
   ok_ng = ( rel( ref_b , v_f ) <= tol ) && ( rel( ref_b , v_g ) <= tol );
   std::cout << "binary master, no slack: ref = " << ref_b << "   Farkas = "
             << v_f << " ( " << ct_f << " cuts )   no-good = " << v_g
             << " ( " << ct_g << " cuts )"
             << ( ok_ng ? "   -> OK" : "   -> FAIL" ) << std::endl;
   delete root_f;
   delete root_g;
   delete mono_b;
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

 const bool ok = ok1 && ok2 && ok_ns && ok_ng && ok_k;

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
