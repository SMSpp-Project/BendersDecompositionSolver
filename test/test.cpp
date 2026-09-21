/*--------------------------------------------------------------------------*/
/*------------------------------ File test.cpp -----------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Validation test for BendersDecompositionSolver on a small two-stage linear
 * program that the test writes itself out of AbstractBlock, so that the
 * module knows of no Block but the ones it builds.
 *
 * Two equivalent models of the same problem are built:
 *
 * - a *monolithic* AbstractBlock (the first-stage Variable y and the
 *   second-stage Variable x in one Block) solved by a :MILPSolver, giving the
 *   reference optimum;
 *
 * - a *structured* AbstractBlock (a master Block with y, one nested sub-Block
 *   with x whose Constraint couple y) solved by
 *   BendersDecompositionSolver.
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
#include <random>

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

// a small instance: M first-stage Variable, N second-stage rows
static const int M = 3;
static const int N = 4;
static const double fixed_cost[ M ] = { 5 , 7 , 6 };
static const double capacity [ M ] = { 3 , 3 , 3 };
static const double demand   [ N ] = { 1 , 1 , 1 , 1 };
static const double cost[ M ][ N ] = { { 2 , 3 , 4 , 5 } ,
                                       { 4 , 1 , 2 , 3 } ,
                                       { 3 , 4 , 1 , 2 } };

// big-M cost of the slack of a row: it makes the subproblem always feasible
// (so only Benders optimality cuts are needed),
// while being large enough that the slack is never used at the optimum
static const double BigM = 1e2;

using array_type = boost::multi_array< ColVariable , 2 >;

/*--------------------------------------------------------------------------*/
/*--------------------------- BLOCK BUILDERS -------------------------------*/
/*--------------------------------------------------------------------------*/

// add to block the second stage of one scenario s ( the Variable x[ s ] and
// its slack, the rows they satisfy and the Constraint that couples them to
// the master Variable y ) and its cost terms to obj. With nsub scenarios this
// gives a 2-stage Benders structure with nsub subproblems coupled through y

static void add_second_stage( AbstractBlock * block ,
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
  add_second_stage( block , y , s , with_slack , f );

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
  add_second_stage( sub , y , s , with_slack , sf );
  auto sobj = new FRealObjective( sub , sf );
  sobj->set_sense( Objective::eMin );
  sub->set_objective( sobj );
  root->add_nested_Block( sub );
  }

 return( root );
 }

/*--------------------------------------------------------------------------*/

/* The same instance, with each scenario built as a tree: the second-stage
 * Variable that one first-stage Variable couples, and their cost, live in a
 * sub-Block of the scenario, which keeps the rows and the coupling
 * Constraint. This is the
 * shape a subproblem has whenever it is a model of its own rather than a bare
 * set of rows, and the Objective it is evaluated with is then the sum of
 * those of the Block it is made of; the separation problems that replicate it
 * have to reckon with all of them, which is what this instance checks. */

static AbstractBlock * build_nested( int nsub = 1 )
{
 auto root = new AbstractBlock();
 auto y = make_y( true , false );
 root->add_static_variable( * y , "y" );

 auto df = new LinearFunction();
 for( int i = 0 ; i < M ; ++i )
  df->add_variable( & ( * y )[ i ] , fixed_cost[ i ] );
 auto robj = new FRealObjective( root , df );
 robj->set_sense( Objective::eMin );
 root->set_objective( robj );

 for( int s = 0 ; s < nsub ; ++s ) {
  const std::string t = std::to_string( s );
  auto sub = new AbstractBlock( root );

  // one sub-Block per first-stage Variable, carrying the second-stage
  // Variable it couples and their cost
  std::vector< std::vector< ColVariable > * > x( M );

  for( int i = 0 ; i < M ; ++i ) {
   auto loc = new AbstractBlock( sub );
   x[ i ] = new std::vector< ColVariable >( N );
   for( auto & v : * x[ i ] )
    v.is_positive( true );
   loc->add_static_variable( * x[ i ] , "x" + t + "_" + std::to_string( i ) );

   auto lf = new LinearFunction();
   for( int j = 0 ; j < N ; ++j )
    lf->add_variable( & ( * x[ i ] )[ j ] , cost[ i ][ j ] );
   auto lobj = new FRealObjective( loc , lf );
   lobj->set_sense( Objective::eMin );
   loc->set_objective( lobj );

   sub->add_nested_Block( loc );
   }

  auto sl = new std::vector< ColVariable >( N );
  for( auto & v : * sl )
   v.is_positive( true );
  sub->add_static_variable( * sl , "s" + t );

  auto dem = new std::vector< FRowConstraint >( N );
  for( int j = 0 ; j < N ; ++j ) {
   auto f = new LinearFunction();
   for( int i = 0 ; i < M ; ++i )
    f->add_variable( & ( * x[ i ] )[ j ] , 1 );
   f->add_variable( & ( * sl )[ j ] , 1 );
   ( * dem )[ j ].set_function( f );
   ( * dem )[ j ].set_both( 1 );
   }
  sub->add_static_constraint( * dem , "demand" + t );

  auto cap = new std::vector< FRowConstraint >( M );
  for( int i = 0 ; i < M ; ++i ) {
   auto f = new LinearFunction();
   for( int j = 0 ; j < N ; ++j )
    f->add_variable( & ( * x[ i ] )[ j ] , demand[ j ] );
   f->add_variable( & ( * y )[ i ] , - capacity[ i ] );
   ( * cap )[ i ].set_function( f );
   ( * cap )[ i ].set_lhs( - Inf< double >() );
   ( * cap )[ i ].set_rhs( 0 );
   }
  sub->add_static_constraint( * cap , "capacity" + t );

  // the scenario itself only pays the slacks, the cost of the second stage being in
  // the Block it is made of
  auto sf = new LinearFunction();
  for( int j = 0 ; j < N ; ++j )
   sf->add_variable( & ( * sl )[ j ] , BigM );
  auto sobj = new FRealObjective( sub , sf );
  sobj->set_sense( Objective::eMin );
  sub->set_objective( sobj );

  root->add_nested_Block( sub );
  }

 return( root );
 }

/*--------------------------------------------------------------------------*/

// a sparse instance: each second-stage row is written over two first-stage
// Variable only, so that an infeasible subproblem has more than one way of
// being so, i.e., its dual has more than one extreme ray, and the coupling
// row i is written multiplied by row_scale[ i ] when scaled is true

static AbstractBlock * build_sparse( unsigned seed , bool scaled ,
				     int nsub = 4 )
{
 const int m = 6 , n = 12;
 std::mt19937 g( seed );
 std::uniform_real_distribution< double > u( 0 , 1 );

 std::vector< double > fc( m ) , cp( m ) , scale( m , 1 );
 std::vector< std::vector< double > > cs( m , std::vector< double >( n ) );
 std::vector< std::vector< int > > arc( m , std::vector< int >( n , 0 ) );
 for( int i = 0 ; i < m ; ++i ) {
  fc[ i ] = 5 + 10 * u( g );
  cp[ i ] = 3 + 4 * u( g );
  }
 for( int j = 0 ; j < n ; ++j ) {
  const int a = g() % m;
  int b = g() % m;
  if( b == a )
   b = ( a + 1 ) % m;
  arc[ a ][ j ] = arc[ b ][ j ] = 1;
  for( int i = 0 ; i < m ; ++i )
   cs[ i ][ j ] = 1 + 5 * u( g );
  }
 std::vector< std::vector< double > > dm( nsub , std::vector< double >( n ) );
 for( auto & d : dm )
  for( auto & dj : d )
   dj = 1 + 2 * u( g );
 if( scaled ) {
  std::uniform_real_distribution< double > e( -2 , 2 );
  for( auto & k : scale )
   k = std::pow( 10 , e( g ) );
  }

 auto root = new AbstractBlock();
 auto y = new std::vector< ColVariable >( m );
 for( auto & yi : * y ) {
  yi.is_unitary( true );
  yi.is_positive( true );
  }
 root->add_static_variable( * y , "y" );
 auto df = new LinearFunction();
 for( int i = 0 ; i < m ; ++i )
  df->add_variable( & ( * y )[ i ] , fc[ i ] );
 auto robj = new FRealObjective( root , df );
 robj->set_sense( Objective::eMin );
 root->set_objective( robj );

 for( int s = 0 ; s < nsub ; ++s ) {
  auto sub = new AbstractBlock( root );
  boost::array< array_type::index , 2 > shape = { m , n };
  auto x = new array_type( shape );
  for( auto q = x->data() ; q < x->data() + x->num_elements() ; ++q )
   q->is_positive( true );
  sub->add_static_variable( * x , "x" );

  auto dem = new std::vector< FRowConstraint >( n );
  for( int j = 0 ; j < n ; ++j ) {
   auto f = new LinearFunction();
   for( int i = 0 ; i < m ; ++i )
    if( arc[ i ][ j ] )
     f->add_variable( & ( * x )[ i ][ j ] , 1 );
   ( * dem )[ j ].set_function( f );
   ( * dem )[ j ].set_both( dm[ s ][ j ] );
   }
  sub->add_static_constraint( * dem , "demand" );

  auto cap = new std::vector< FRowConstraint >( m );
  for( int i = 0 ; i < m ; ++i ) {
   auto f = new LinearFunction();
   for( int j = 0 ; j < n ; ++j )
    if( arc[ i ][ j ] )
     f->add_variable( & ( * x )[ i ][ j ] , scale[ i ] );
   f->add_variable( & ( * y )[ i ] , - scale[ i ] * cp[ i ] );
   ( * cap )[ i ].set_function( f );
   ( * cap )[ i ].set_lhs( - Inf< double >() );
   ( * cap )[ i ].set_rhs( 0 );
   }
  sub->add_static_constraint( * cap , "capacity" );

  auto sf = new LinearFunction();
  for( int i = 0 ; i < m ; ++i )
   for( int j = 0 ; j < n ; ++j )
    if( arc[ i ][ j ] )
     sf->add_variable( & ( * x )[ i ][ j ] , cs[ i ][ j ] );
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
 /* Whatever goes wrong, from the configuration being refused onwards, the
  * Block is given back what the Solver had taken from it before the error
  * travels on, so that the caller can go on using the Block. */

 auto give_back = [ & ]() { bsc->clear(); bsc->apply( block ); delete bsc; };

 double lb;
 Solver * solver;
 try {
  bsc->apply( block );
  solver = block->get_registered_solvers().front();
  status = solver->compute( false );
  lb = solver->get_lb();
  }
 catch( ... ) { give_back(); throw; }

 if( iters )
  *iters = solver->get_elapsed_iterations();
 if( cuts )
  if( auto bds = dynamic_cast< BendersDecompositionSolver * >( solver ) )
   *cuts = bds->get_num_cuts();
 /* Reading the bound is all that was needed: the Solver is un-registered
  * and deleted by applying the cleared BlockSolverConfig, which is what
  * gives the Block back whatever the Solver had taken from it. */

 give_back();
 return( lb );
 }

/*--------------------------------------------------------------------------*/
/* The convex regime attaches a BundleSolver to the master, which keeps the
 * subproblems hanging from it as sub-Block; a bundle that wants each of its
 * sub-Block to be a bare function Block refuses it in set_Block(), and so
 * does one that does not carry the required parameters. In either case the
 * case is skipped rather than failed, and the reference value is returned so
 * that the comparisons downstream are satisfied. */

static double solve_by_bundle( AbstractBlock * block , const std::string & fn ,
			       int & status , bool & ran , double dflt )
{
 try {
  const double lb = solve_from_config( block , fn , status );
  ran = true;
  return( lb );
  }
 catch( const std::exception & e ) {
  std::cout << "Benders(convex): skipped, " << e.what() << std::endl;
  ran = false;
  return( dflt );
  }
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
 int st = 0;
 bool has_cvx;
 const double ben = solve_by_bundle( root , "BSPar_benders_convex.txt" , st ,
				     has_cvx , ref );
 if( has_cvx )
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
 // two second-stage scenarios coupled through y: BDS builds one
 // BendersBFunction per scenario, so multi-cut uses two epigraph Variable and
 // single-cut one; all must match the 2-scenario monolithic optimum
 auto mono2 = build_monolithic( true , 2 );
 int dummy;
 const double ref2 = solve_from_config( mono2 , "BSPar_sub.txt" , dummy );
 auto root_c2 = build_structured( true , 2 );
 bool has_cvx2;
 const double ben_c2 = solve_by_bundle( root_c2 , "BSPar_benders_convex.txt" ,
					dummy , has_cvx2 , ref2 );
 auto root_m2 = build_structured( true , 2 );
 const double ben_m2 = solve_from_config( root_m2 , "BSPar_benders_milp.txt" ,
					  dummy );
 auto root_s2 = build_structured( true , 2 );
 const double ben_s2 = solve_from_config( root_s2 ,
					  "BSPar_benders_milp_single.txt" , dummy );
 std::cout << "2-scenario: ref = " << ref2;
 if( has_cvx2 )
  std::cout << "   convex = " << ben_c2;
 std::cout << "   MILP-multi = " << ben_m2 << "   MILP-single = " << ben_s2
           << std::endl;

 /* ----- how many cuts each variant of the MILP regime takes -------------- #
  *
  * The three of them describe the same problem and have to end at the same
  * value: what changes is how many cuts are needed to get there, which is the
  * figure to compare, the number of rounds saying little when one round adds
  * one cut and another adds one per subproblem. */

 bool ok_t = true;
 { auto root_m = build_structured( true , 4 );
   auto root_s = build_structured( true , 4 );
   auto root_p = build_structured( true , 4 );
   int st_m , st_s , st_p;
   long it_m = 0 , it_s = 0 , it_p = 0 , ct_m = 0 , ct_s = 0 , ct_p = 0;

   const double v_m = solve_from_config( root_m , "BSPar_benders_milp.txt" ,
					 st_m , & it_m , & ct_m );

   /* The subproblems of a round evaluated by threads: the cuts are added in
    * the order of the subproblems whatever the order the threads finish in,
    * hence the master, and with it every round, is the one of the sequential
    * run, which is what is asserted. */

   auto root_t = build_structured( true , 4 );
   int st_t;
   long it_t = 0 , ct_t = 0;
   const double v_t = solve_from_config( root_t , "BSPar_benders_milp_par.txt" ,
					 st_t , & it_t , & ct_t );
   ok_t = ( v_t == v_m ) && ( it_t == it_m ) && ( ct_t == ct_m );
   std::cout << "4-scenario MILP master, 4 threads: " << v_t << " ( " << it_t
             << " rounds , " << ct_t << " cuts )"
             << ( ok_t ? "   -> OK" : "   -> FAIL" ) << std::endl;
   delete root_t;
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

  /* The same with the subproblems evaluated by threads, where some of them
   * are infeasible and their certificates are read inside the threads: the
   * run has to be the sequential one, round by round. */

  auto root_ft4 = build_structured( false , 4 );
  int st_ft4;
  long it_ft4 = 0 , ct_ft4 = 0;
  const double v_ft4 = solve_from_config( root_ft4 ,
					  "BSPar_benders_milp_par.txt" ,
					  st_ft4 , & it_ft4 , & ct_ft4 );
  const bool ok_ft = ( v_ft4 == v_f4 ) && ( it_ft4 == it_f4 ) &&
                     ( ct_ft4 == ct_f4 );
  std::cout << "4-scenario, no slack, 4 threads: " << v_ft4 << " ( "
            << it_ft4 << " rounds , " << ct_ft4 << " cuts )"
            << ( ok_ft ? "   -> OK" : "   -> FAIL" ) << std::endl;
  ok_ns = ok_ns && ok_ft;
  delete root_ft4;
  /* The unified cut on the same instance, which is where it has something to
   * do: the subproblems are infeasible at some x, and the one cut it
   * separates is a feasibility cut there and an optimality cut elsewhere,
   * without the two mechanisms the other two runs need. */

  auto root_u4 = build_structured( false , 4 );
  int st_u4;
  long it_u4 = 0 , ct_u4 = 0;
  const double v_u4 = solve_from_config( root_u4 ,
					 "BSPar_benders_milp_unified.txt" ,
					 st_u4 , & it_u4 , & ct_u4 );
  ok_ns = ok_ns && ( rel( ref4 , v_u4 ) <= tol );
  delete root_u4;

  /* The same cut normalized on its own coefficients rather than on the
   * slacks, which is the deepest one: here there is nothing to set, the
   * separation problem measuring how far the incumbent is from the epigraph
   * instead of how much the rows are violated. */

  auto root_d4 = build_structured( false , 4 );
  int st_d4;
  long it_d4 = 0 , ct_d4 = 0;
  const double v_d4 = solve_from_config( root_d4 ,
					 "BSPar_benders_milp_deepest.txt" ,
					 st_d4 , & it_d4 , & ct_d4 );
  ok_ns = ok_ns && ( rel( ref4 , v_d4 ) <= tol );
  delete root_d4;

  /* The same cut with the normalization written as the literature writes it,
   * one equation over the multipliers in place of a bound on each of them,
   * and with the weights of the static cut of Brandenberg and Stursberg. */

  auto root_s4 = build_structured( false , 4 );
  int st_s4;
  long it_s4 = 0 , ct_s4 = 0;
  const double v_s4 = solve_from_config( root_s4 ,
					 "BSPar_benders_milp_static.txt" ,
					 st_s4 , & it_s4 , & ct_s4 );
  ok_ns = ok_ns && ( rel( ref4 , v_s4 ) <= tol );
  delete root_s4;

  std::cout << "4-scenario, no slack: ref = " << ref4 << "   Farkas = "
            << v_f4 << " ( " << it_f4 << " rounds , " << ct_f4
            << " cuts )   phase one = " << v_14 << " ( " << it_14
            << " rounds , " << ct_14 << " cuts )   unified = " << v_u4
            << " ( " << it_u4 << " rounds , " << ct_u4 << " cuts )"
            << "   deepest = " << v_d4 << " ( " << it_d4 << " rounds , "
            << ct_d4 << " cuts )   static = " << v_s4 << " ( " << it_s4
            << " rounds , " << ct_s4 << " cuts )" << std::endl;
  delete root_f4;
  delete root_14;

  /* The phase one on the sparse instance, whose infeasible subproblems have
   * more than one extreme ray, hence where the normalization decides which
   * cut is taken: with the slacks costing one, writing the coupling rows
   * with factors between 1/100 and 100 changes the cuts, while with the
   * slacks of each row costing the inverse of its norm it does not, the
   * problem being the same. What is asserted is the optimum everywhere and
   * the invariance of the weighted run; the unit costs are printed only. */

  auto mono_sp = build_sparse( 4 , false );
  int st_sp;
  const double ref_sp = solve_from_config( mono_sp , "BSPar_sub.txt" , st_sp );
  delete mono_sp;

  struct { const char * cfg; bool scaled; double v; long it , ct; } sp[] = {
   { "BSPar_benders_milp_phase1.txt" , false , 0 , 0 , 0 } ,
   { "BSPar_benders_milp_phase1.txt" , true , 0 , 0 , 0 } ,
   { "BSPar_benders_milp_phase1_norm.txt" , false , 0 , 0 , 0 } ,
   { "BSPar_benders_milp_phase1_norm.txt" , true , 0 , 0 , 0 } };
  for( auto & r : sp ) {
   auto root = build_sparse( 4 , r.scaled );
   int st;
   r.v = solve_from_config( root , r.cfg , st , & r.it , & r.ct );
   delete root;
   }
  const bool ok_w = ( rel( ref_sp , sp[ 0 ].v ) <= tol ) &&
                    ( rel( ref_sp , sp[ 1 ].v ) <= tol ) &&
                    ( rel( ref_sp , sp[ 2 ].v ) <= tol ) &&
                    ( rel( ref_sp , sp[ 3 ].v ) <= tol ) &&
                    ( sp[ 2 ].it == sp[ 3 ].it ) &&
                    ( sp[ 2 ].ct == sp[ 3 ].ct );
  ok_ns = ok_ns && ok_w;
  std::cout << "sparse, phase one: ref = " << ref_sp << "   unit costs = "
            << sp[ 0 ].it << " rounds , rows scaled = " << sp[ 1 ].it
            << " rounds   row-norm costs = " << sp[ 2 ].it
            << " rounds , rows scaled = " << sp[ 3 ].it << " rounds"
            << ( ok_w ? "   -> OK" : "   -> FAIL" ) << std::endl;
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

 /* ----- the unified cut ------------------------------------------------- #
  *
  * One cut for feasibility and optimality, separated by the phase one with
  * the epigraph inequality among the rows that carry a slack: it has to end
  * at the same optimum as the ordinary cuts, on the flat instance and on the
  * one whose scenario is a tree, where the cost the epigraph inequality
  * bounds is the sum of those of three sub-Block. */

 bool ok_u = true;
 { auto root_u = build_structured( true , 2 );
   int st_u;
   long it_u = 0 , ct_u = 0;
   const double v_u = solve_from_config( root_u ,
					 "BSPar_benders_milp_unified.txt" ,
					 st_u , & it_u , & ct_u );
   delete root_u;

   auto mono_n = build_nested( 2 );
   int st_n;
   const double ref_n = solve_from_config( mono_n , "BSPar_sub.txt" , st_n );

   auto root_n = build_nested( 2 );
   int st_nu;
   long it_nu = 0 , ct_nu = 0;
   const double v_nu = solve_from_config( root_n ,
					  "BSPar_benders_milp_unified.txt" ,
					  st_nu , & it_nu , & ct_nu );
   delete root_n;

   auto root_nm = build_nested( 2 );
   int st_nm;
   const double v_nm = solve_from_config( root_nm , "BSPar_benders_milp.txt" ,
					  st_nm );
   delete root_nm;
   delete mono_n;

   ok_u = ( rel( ref2 , v_u ) <= tol ) && ( rel( ref_n , v_nu ) <= tol ) &&
          ( rel( ref_n , v_nm ) <= tol );

   std::cout << "unified: flat = " << v_u << " ( " << it_u << " rounds , "
             << ct_u << " cuts , ref " << ref2 << " )   nested = " << v_nu
             << " ( " << it_nu << " rounds , " << ct_nu << " cuts , ref "
             << ref_n << " , ordinary " << v_nm << " )"
             << ( ok_u ? "   -> OK" : "   -> FAIL" ) << std::endl;
   }

 /* ----- the Block given back ------------------------------------------- #
  *
  * With the reformulation undone at the end of every compute() the Block is
  * a problem of its own in between two of them, which is what lets another
  * Solver be attached to it: the same Block is therefore solved again, by a
  * :MILPSolver reading the whole tree, i.e., the extensive form, and the two
  * have to agree. Were the complicating Variable not given back to the
  * Constraint they were taken out of, the second solve would be of a problem
  * in which the subproblems are free of the master, hence of a relaxation. */

 bool ok_r = true;
 { auto root_r = build_structured( true , 2 );
   int st_r;
   long it_r = 0 , ct_r = 0;
   const double v_r = solve_from_config( root_r ,
					 "BSPar_benders_milp_restore.txt" ,
					 st_r , & it_r , & ct_r );

   int st_a;
   const double v_a = solve_from_config( root_r , "BSPar_sub.txt" , st_a );
   delete root_r;

   ok_r = ( rel( ref2 , v_r ) <= tol ) && ( rel( ref2 , v_a ) <= tol );

   std::cout << "Block given back: Benders = " << v_r << " ( " << it_r
             << " rounds , " << ct_r << " cuts )   the same Block then = "
             << v_a << "   ref = " << ref2
             << ( ok_r ? "   -> OK" : "   -> FAIL" ) << std::endl;

   /* The same, with the two Solver registered together rather than one after
    * the other, which is what a BlockSolverConfig of a cross-check does: each
    * of them is computed twice, in both orders, since what one does to the
    * Block the other has to see, or not see, whenever it is asked. */

   auto root_x = build_structured( true , 2 );
   auto cfg = Configuration::deserialize( "BSPar_benders_milp_cross.txt" );
   auto bsc = dynamic_cast< BlockSolverConfig * >( cfg );
   if( ! bsc ) {
    delete cfg;
    std::cerr << "BSPar_benders_milp_cross.txt is not a BlockSolverConfig"
              << std::endl;
    std::exit( 1 );
    }

   bsc->apply( root_x );
   auto & slvrs = root_x->get_registered_solvers();
   auto sb = slvrs.front();
   auto sm = slvrs.back();

   sb->compute( false );
   const double x_b1 = sb->get_lb();
   sm->compute( false );
   const double x_m1 = sm->get_lb();
   sb->compute( false );
   const double x_b2 = sb->get_lb();
   sm->compute( false );
   const double x_m2 = sm->get_lb();

   bsc->clear();
   bsc->apply( root_x );
   delete bsc;
   delete root_x;

   const bool ok_x = ( rel( ref2 , x_b1 ) <= tol ) &&
                     ( rel( ref2 , x_m1 ) <= tol ) &&
                     ( rel( ref2 , x_b2 ) <= tol ) &&
                     ( rel( ref2 , x_m2 ) <= tol );
   ok_r = ok_r && ok_x;

   std::cout << "the two together: Benders = " << x_b1 << " , " << x_b2
             << "   :MILPSolver = " << x_m1 << " , " << x_m2 << "   ref = "
             << ref2 << ( ok_x ? "   -> OK" : "   -> FAIL" ) << std::endl;
   }

 // ----- compare ( optimality-cut cases, the supported ones ) ------------- #
 const double err = rel( ref , ben );
 const double err2 = rel( ref , ben2 );
 const double err3 = rel( ref , ben3 );
 const bool ok1 = ( err <= tol ) && ( err2 <= tol ) && ( err3 <= tol );
 std::cout << "1-scenario: monolithic = " << ref;
 if( has_cvx )
  std::cout << "   convex = " << ben << " (err " << err << ")";
 std::cout << "   MILP-multi = " << ben2 << " (err " << err2 << ")"
           << "   MILP-single = " << ben3 << " (err " << err3 << ")"
           << ( ok1 ? "   -> OK" : "   -> FAIL" ) << std::endl;

 const bool ok2 = ( rel( ref2 , ben_c2 ) <= tol )
	       && ( rel( ref2 , ben_m2 ) <= tol )
	       && ( rel( ref2 , ben_s2 ) <= tol );
 std::cout << "2-scenario: " << ( ok2 ? "-> OK" : "-> FAIL" ) << std::endl;

 const bool ok = ok1 && ok2 && ok_ns && ok_ng && ok_k && ok_u && ok_r &&
                ok_t;

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
