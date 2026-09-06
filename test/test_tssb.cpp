/*--------------------------------------------------------------------------*/
/*--------------------------- File test_tssb.cpp ---------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * The two Benders forms of the same two-stage stochastic investment problem,
 * one against the other on the same data.
 *
 * The ad hoc one is an InvestmentBlock whose inner Block is the whole
 * TwoStageStochasticBlock: its InvestmentFunction is the value function, and
 * it FIXES the design in every scenario rather than mapping it into the
 * right-hand side, so it produces one aggregated linearization per iteration.
 *
 * The generic one is the structure this Solver asks for, i.e., the here-and-now
 * Variable in a single copy in the root, with no stochastic element in it, and
 * one LP sub-Block per scenario coupled to them only through its capacity
 * Constraint; it is built here out of the very same instance data, and the
 * BendersBFunction, hence the multi-cut, comes from
 * BendersDecompositionSolver.
 *
 * Being two formulations of one problem, the two optima must coincide, and
 * both must coincide with the monolithic solve of the extensive form, which is
 * the reference; the iterations and the times are reported.
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
#include <netcdf>

#include "AbstractBlock.h"
#include "BendersDecompositionSolver.h"
#include "BlockSolverConfig.h"
#include "ColVariable.h"
#include "FRealObjective.h"
#include "FRowConstraint.h"
#include "InvestmentBlock.h"
#include "InvestmentFunction.h"
#include "LinearFunction.h"
#include "Objective.h"
#include "OneVarConstraint.h"

using namespace SMSpp_di_unipi_it;
using namespace netCDF;

/*--------------------------------------------------------------------------*/
/*----------------------------- CONSTANTS ----------------------------------*/
/*--------------------------------------------------------------------------*/

// the relative gap at which two of these solves are the same number: the ad
// hoc form reads its objective back through a bundle whose required relative
// accuracy is 1e-4, hence nothing tighter can be asked of the comparison
static constexpr double tolerance = 1e-6;

static const std::string instance = "tssb_investment.nc4";

/*--------------------------------------------------------------------------*/
/*------------------------------- TYPES ------------------------------------*/
/*--------------------------------------------------------------------------*/

/// the instance, read once and used to build every form of it
/** Whatever the form, the data is the same, and it is read from the very file
 * the ad hoc form is deserialized from, so that the two cannot drift apart. */

struct Data {
 double cost;                       ///< cost of one unit of design
 double lb;                         ///< lower bound on the design
 double ub;                         ///< upper bound on the design
 std::vector< double > weight;      ///< one probability per scenario
 std::vector< std::vector< double > > demand;   ///< [ scenario ][ time ]
 std::vector< double > unit_max;    ///< max power, per unit of design
 std::vector< double > unit_cost;   ///< cost of the design-carrying unit
 std::vector< double > slack_max;   ///< max power of the slack
 std::vector< double > slack_cost;  ///< cost of the slack
 };

/*--------------------------------------------------------------------------*/
/*----------------------------- READING ------------------------------------*/
/*--------------------------------------------------------------------------*/

static std::vector< double > get_vec( const NcGroup & g ,
				      const std::string & name )
{
 auto v = g.getVar( name );
 std::size_t n = 1;
 for( int i = 0 ; i < v.getDimCount() ; ++i )
  n *= v.getDim( i ).getSize();
 std::vector< double > out( n );
 v.getVar( out.data() );
 return( out );
 }

/*--------------------------------------------------------------------------*/

static Data read( const std::string & file )
{
 Data d;
 NcFile f( file , NcFile::read );

 auto inv = f.getGroup( "Block_0" );
 d.cost = get_vec( inv , "Cost" ).front();
 d.lb = get_vec( inv , "LowerBound" ).front();
 d.ub = get_vec( inv , "UpperBound" ).front();

 auto tssb = inv.getGroup( "InnerBlock" );
 auto set = tssb.getGroup( "DiscreteScenarioSet" );
 d.weight = get_vec( set , "PoolWeights" );

 auto sc = set.getVar( "Scenarios" );
 const auto K = sc.getDim( 0 ).getSize() , T = sc.getDim( 1 ).getSize();
 auto flat = get_vec( set , "Scenarios" );
 d.demand.resize( K );
 for( std::size_t k = 0 ; k < K ; ++k )
  d.demand[ k ].assign( flat.begin() + k * T , flat.begin() + ( k + 1 ) * T );

 auto uc = tssb.getGroup( "StochasticBlock" ).getGroup( "Block" );
 auto u0 = uc.getGroup( "UnitBlock_0" ) , u1 = uc.getGroup( "UnitBlock_1" );
 d.unit_max = get_vec( u0 , "MaxPower" );
 d.unit_cost = get_vec( u0 , "ActivePowerCost" );
 d.slack_max = get_vec( u1 , "MaxPower" );
 d.slack_cost = get_vec( u1 , "ActivePowerCost" );

 return( d );
 }

/*--------------------------------------------------------------------------*/
/*--------------------------- BLOCK BUILDERS -------------------------------*/
/*--------------------------------------------------------------------------*/

/* BundleSolver does not take a general bound l <= x <= u on the Variable of
 * the Block it is attached to, so the design is shifted: the Variable is
 * x' = x - l in [ 0 , u - l ], and l is added back wherever x appears, which
 * leaves a constant cost * l in the Objective. The ad hoc form does the very
 * same thing, through the BlockConfig of the InvestmentBlock, so the two keep
 * describing one problem. */

static std::vector< ColVariable > * make_x( const Data & d )
{
 auto x = new std::vector< ColVariable >( 1 );
 ( * x )[ 0 ].is_positive( true );
 ( * x )[ 0 ].set_value( 0 );
 return( x );
 }

/*--------------------------------------------------------------------------*/

// the box on the shifted design, the only Constraint of the master

static void add_design_bound( AbstractBlock * block , const Data & d ,
			      std::vector< ColVariable > * x )
{
 auto bnd = new std::vector< BoxConstraint >( 1 );
 ( * bnd )[ 0 ].set_variable( & ( * x )[ 0 ] );
 ( * bnd )[ 0 ].set_lhs( 0 );
 ( * bnd )[ 0 ].set_rhs( d.ub - d.lb );
 block->add_static_constraint( * bnd , "design bound" );
 }

/*--------------------------------------------------------------------------*/

/* add to block the dispatch of one scenario: the power of the design-carrying
 * unit and of the slack, the demand Constraint and the capacity Constraint,
 * the latter being the only place where the design Variable appears, and the
 * costs of the scenario, weighted by its probability, to obj */

static void add_scenario( AbstractBlock * block , const Data & d ,
			  std::size_t k , std::vector< ColVariable > * x ,
			  LinearFunction * obj )
{
 const auto T = d.demand[ k ].size();

 auto g = new std::vector< ColVariable >( T );
 for( auto & v : * g ) v.is_positive( true );
 block->add_static_variable( * g , "g" );

 auto s = new std::vector< ColVariable >( T );
 for( auto & v : * s ) v.is_positive( true );
 block->add_static_variable( * s , "s" );

 auto dem = new std::vector< FRowConstraint >( T );
 for( std::size_t t = 0 ; t < T ; ++t ) {
  auto f = new LinearFunction();
  f->add_variable( & ( * g )[ t ] , 1 );
  f->add_variable( & ( * s )[ t ] , 1 );
  ( * dem )[ t ].set_function( f );
  ( * dem )[ t ].set_both( d.demand[ k ][ t ] );
  }
 block->add_static_constraint( * dem , "demand" );

 auto cap = new std::vector< FRowConstraint >( T );
 for( std::size_t t = 0 ; t < T ; ++t ) {
  auto f = new LinearFunction();
  f->add_variable( & ( * g )[ t ] , 1 );
  f->add_variable( & ( * x )[ 0 ] , - d.unit_max[ t ] );
  ( * cap )[ t ].set_function( f );
  ( * cap )[ t ].set_lhs( -Inf< double >() );
  ( * cap )[ t ].set_rhs( d.unit_max[ t ] * d.lb );
  }
 block->add_static_constraint( * cap , "capacity" );

 auto slk = new std::vector< FRowConstraint >( T );
 for( std::size_t t = 0 ; t < T ; ++t ) {
  auto f = new LinearFunction();
  f->add_variable( & ( * s )[ t ] , 1 );
  ( * slk )[ t ].set_function( f );
  ( * slk )[ t ].set_lhs( -Inf< double >() );
  ( * slk )[ t ].set_rhs( d.slack_max[ t ] );
  }
 block->add_static_constraint( * slk , "slack capacity" );

 for( std::size_t t = 0 ; t < T ; ++t ) {
  obj->add_variable( & ( * g )[ t ] , d.weight[ k ] * d.unit_cost[ t ] );
  obj->add_variable( & ( * s )[ t ] , d.weight[ k ] * d.slack_cost[ t ] );
  }
 }

/*--------------------------------------------------------------------------*/

// the monolithic model, every scenario in one Block: the reference optimum

static AbstractBlock * build_monolithic( const Data & d )
{
 auto block = new AbstractBlock();
 auto x = make_x( d );
 block->add_static_variable( * x , "x" );
 add_design_bound( block , d , x );

 auto f = new LinearFunction();
 f->add_variable( & ( * x )[ 0 ] , d.cost );
 f->set_constant_term( d.cost * d.lb );
 for( std::size_t k = 0 ; k < d.demand.size() ; ++k )
  add_scenario( block , d , k , x , f );

 auto obj = new FRealObjective( block , f );
 obj->set_sense( Objective::eMin );
 block->set_objective( obj );
 return( block );
 }

/*--------------------------------------------------------------------------*/

// the 2-level model: the design in the root, one LP sub-Block per scenario

static AbstractBlock * build_structured( const Data & d )
{
 auto root = new AbstractBlock();
 auto x = make_x( d );
 root->add_static_variable( * x , "x" );
 add_design_bound( root , d , x );

 auto df = new LinearFunction();
 df->add_variable( & ( * x )[ 0 ] , d.cost );
 df->set_constant_term( d.cost * d.lb );
 auto robj = new FRealObjective( root , df );
 robj->set_sense( Objective::eMin );
 root->set_objective( robj );

 for( std::size_t k = 0 ; k < d.demand.size() ; ++k ) {
  auto sub = new AbstractBlock( root );
  auto sf = new LinearFunction();
  add_scenario( sub , d , k , x , sf );
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

static double solve( Block * block , const std::string & cfg , double & secs ,
		     long & iters , int & status )
{
 auto c = Configuration::deserialize( cfg );
 auto bsc = dynamic_cast< BlockSolverConfig * >( c );
 if( ! bsc ) {
  std::cerr << cfg << " is not a BlockSolverConfig" << std::endl;
  std::exit( 1 );
  }
 bsc->apply( block );

 auto solver = block->get_registered_solvers().front();

 auto t0 = std::chrono::steady_clock::now();
 status = solver->compute();
 secs = std::chrono::duration< double >(
                              std::chrono::steady_clock::now() - t0 ).count();
 iters = solver->get_elapsed_iterations();
 const auto value = solver->get_var_value();

 bsc->clear();
 bsc->apply( block );
 delete bsc;
 return( value );
 }

/*--------------------------------------------------------------------------*/

/* the ad hoc form: the InvestmentBlock is deserialized from the file, its
 * bound Constraint are reformulated into the shifted box that a BundleSolver
 * takes, and the BlockSolverConfig of the inner Block travels to the
 * InvestmentFunction as the "extra" Configuration it expects */

static double solve_ad_hoc( double & secs , long & iters , int & status )
{
 auto root = Block::deserialize( instance );
 auto inv = dynamic_cast< InvestmentBlock * >( root );
 if( ! inv ) {
  std::cerr << instance << " does not hold an InvestmentBlock" << std::endl;
  std::exit( 1 );
  }

 auto config = new BlockConfig;
 config->f_static_constraints_Configuration =
                                        new SimpleConfiguration< int >( 1 );
 inv->set_BlockConfig( config );

 auto inner = dynamic_cast< BlockSolverConfig * >(
                       Configuration::deserialize( "TSSBInv_BSCfg.txt" ) );
 if( ! inner ) {
  std::cerr << "TSSBInv_BSCfg.txt is not a BlockSolverConfig" << std::endl;
  std::exit( 1 );
  }

 ComputeConfig cc;
 cc.f_extra_Configuration =
  new SimpleConfiguration< std::map< std::string , Configuration * > >
                                       ( { { "BlockSolverConfig" , inner } } );

 auto function = dynamic_cast< InvestmentFunction * >( inv->get_function() );
 if( ! function ) {
  std::cerr << "the InvestmentBlock has no InvestmentFunction" << std::endl;
  std::exit( 1 );
  }
 function->set_ComputeConfig( & cc );

 const auto value = solve( inv , "TSSBInv_BSPar.txt" , secs , iters , status );
 delete root;
 return( value );
 }

/*--------------------------------------------------------------------------*/
/*-------------------------------- MAIN ------------------------------------*/
/*--------------------------------------------------------------------------*/

int main( void )
{
 auto d = read( instance );
 std::cout << d.demand.size() << " scenarios, " << d.demand[ 0 ].size()
	   << " time steps, design cost " << d.cost << " in [ " << d.lb
	   << " , " << d.ub << " ]" << std::endl;

 double t_ref , t_bds , t_inv;
 long i_ref , i_bds , i_inv;
 int s_ref , s_bds , s_inv;

 auto mono = build_monolithic( d );
 const auto ref = solve( mono , "BSPar_sub.txt" , t_ref , i_ref , s_ref );

 auto bend = build_structured( d );
 const auto bds = solve( bend , "BSPar_benders_convex.txt" , t_bds , i_bds ,
			 s_bds );

 const auto inv = solve_ad_hoc( t_inv , i_inv , s_inv );

 std::cout.precision( 12 );
 std::cout << "\nmonolithic        " << ref << "  status " << s_ref
	   << "  in " << t_ref << " s"
	   << "\nBenders, generic  " << bds << "  status " << s_bds
	   << "  in " << t_bds << " s, " << i_bds << " iterations"
	   << "\nBenders, ad hoc   " << inv << "  status " << s_inv
	   << "  in " << t_inv << " s, " << i_inv << " iterations"
	   << std::endl;

 const auto scale = std::max( 1.0 , std::abs( ref ) );
 bool ok = true;

 if( std::abs( bds - ref ) > tolerance * scale ) {
  std::cout << "ERROR: the generic Benders form is off by "
	    << std::abs( bds - ref ) / scale << std::endl;
  ok = false;
  }

 if( std::abs( inv - ref ) > tolerance * scale ) {
  std::cout << "ERROR: the ad hoc Benders form is off by "
	    << std::abs( inv - ref ) / scale << std::endl;
  ok = false;
  }

 std::cout << ( ok ? "\nOK, the two forms agree with the extensive optimum"
	           : "\nthe forms do not describe the same problem" )
	   << std::endl;

 delete mono;
 delete bend;
 return( ok ? 0 : 1 );
 }

/*--------------------------------------------------------------------------*/
/*------------------------- End File test_tssb.cpp -------------------------*/
/*--------------------------------------------------------------------------*/
