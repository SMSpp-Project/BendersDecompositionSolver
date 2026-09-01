/*--------------------------------------------------------------------------*/
/*------------------ File BendersDecompositionSolver.cpp -------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Implementation of the BendersDecompositionSolver class.
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

#include "BendersDecompositionSolver.h"

#include "BlockSolverConfig.h"

#include "FRealObjective.h"

#include "FRowConstraint.h"

#include "DQuadFunction.h"

#include "LinearFunction.h"

#include "SMSTypedefs.h"

/*--------------------------------------------------------------------------*/
/*-------------------------- NAMESPACE & USING -----------------------------*/
/*--------------------------------------------------------------------------*/

using namespace SMSpp_di_unipi_it;

/*--------------------------------------------------------------------------*/
/*----------------------------- STATIC MEMBERS -----------------------------*/
/*--------------------------------------------------------------------------*/

// register BendersDecompositionSolver to the Solver factory

SMSpp_insert_in_factory_cpp_0( BendersDecompositionSolver );

/*--------------------------------------------------------------------------*/

// names of the BendersDecompositionSolver-specific parameters

static const std::vector< std::string > int_pars_BDSlv = {
 "int_BDSlv_iBCopy" ,
 "int_BDSlv_Regime" ,
 "int_BDSlv_CutType" ,
 "int_BDSlv_FeasCut" ,
 "int_BDSlv_MaxRounds"
 };

static const std::vector< std::string > str_pars_BDSlv = {
 "str_BDSlv_MSName" ,
 "str_Bsub_BSCfg" ,
 "str_Mstr_BSCfg"
 };

/*--------------------------------------------------------------------------*/
/*--------------- CONSTRUCTING AND DESTRUCTING -----------------------------*/
/*--------------------------------------------------------------------------*/

BendersDecompositionSolver::BendersDecompositionSolver( void )
 : CDASolver() {}

/*--------------------------------------------------------------------------*/

BendersDecompositionSolver::~BendersDecompositionSolver( void )
{
 /* The reformulation transforms (B) rather than copying it: the x terms are
  * stripped out of the Constraint of the subproblems, which are moved inside
  * the BendersBFunction, and the master is (B) itself. Undoing all this is
  * not attempted, the Block being *consumed* by the Solver; what is done here
  * is releasing what this Solver owns, and only that.
  *
  * Who owns the BendersBFunction depends on the regime: in the convex one
  * each of them is the Function of the Objective of a sub-Block of the
  * master, hence it goes with the master; in the MILP one nobody else has
  * them. In both cases deleting a BendersBFunction deletes the subproblem it
  * holds. The epigraph Variable and the cuts, instead, are referred to by the
  * master, which outlives this Solver, and are therefore left where they
  * are. */

 if( f_master_solver ) {
  if( f_master )
   f_master->unregister_Solver( f_master_solver );
  delete f_master_solver;
  f_master_solver = nullptr;
  }

 if( f_regime == eMILPMaster )
  for( auto bf : v_BF )
   delete bf;

 v_BF.clear();

 }  // end( BendersDecompositionSolver::~BendersDecompositionSolver )

/*--------------------------------------------------------------------------*/

void BendersDecompositionSolver::set_Block( Block * block )
{
 if( block == f_Block )  // nothing to do
  return;

 CDASolver::set_Block( block );

 if( block )
  reformulate();
 }

/*--------------------------------------------------------------------------*/
/*------------------------ SOLVING THE Block -------------------------------*/
/*--------------------------------------------------------------------------*/

int BendersDecompositionSolver::compute( bool changedvars )
{
 /* The reformulation is done here, rather than in set_Block(), because it
  * depends on the parameters, and those are set through a ComputeConfig
  * after the Solver has been registered to the Block. */

 reformulate();

 acquire_master_solver();

 f_solved = false;

 if( f_regime == eMILPMaster ) {
  const int status = solve_MILP_master();
  map_back_solution();
  return( status );
  }

 /* In the convex regime the master Solver is a bundle-type one, which drives
  * the cutting-plane loop by itself: the linearizations it asks the value
  * functions for *are* the Benders cuts. */

 const int status = f_master_solver->compute( changedvars );

 f_value = f_master_solver->get_lb();
 f_solved = f_master_solver->has_var_solution();

 map_back_solution();

 return( status );

 }  // end( BendersDecompositionSolver::compute )

/*--------------------------------------------------------------------------*/
/*------------------------ READING THE SOLUTION ----------------------------*/
/*--------------------------------------------------------------------------*/

Solver::OFValue BendersDecompositionSolver::get_lb( void )
{
 if( f_regime == eMILPMaster )
  return( f_value );

 return( f_master_solver ? f_master_solver->get_lb() : - Inf< OFValue >() );

 }  // end( BendersDecompositionSolver::get_lb )

/*--------------------------------------------------------------------------*/

Solver::OFValue BendersDecompositionSolver::get_ub( void )
{
 if( f_regime == eMILPMaster )
  return( f_solved ? f_value : Inf< OFValue >() );

 return( f_master_solver ? f_master_solver->get_ub() : Inf< OFValue >() );

 }  // end( BendersDecompositionSolver::get_ub )

/*--------------------------------------------------------------------------*/

bool BendersDecompositionSolver::has_var_solution( void )
{
 if( f_regime == eMILPMaster )
  return( f_solved );

 return( f_master_solver && f_master_solver->has_var_solution() );

 }  // end( BendersDecompositionSolver::has_var_solution )

/*--------------------------------------------------------------------------*/

void BendersDecompositionSolver::get_var_solution( Configuration * solc )
{
 /* The x are written into (B) by the master Solver, since the master *is*
  * (B); the y^k of each subproblem are written into the inner Block, which
  * holds the very Variable the original sub-Block had. */

 if( f_master_solver )
  f_master_solver->get_var_solution( solc );

 for( auto bf : v_BF )
  if( auto sub = bf->get_inner_block() )
   for( auto slvr : sub->get_registered_solvers() )
    if( slvr->has_var_solution() ) {
     slvr->get_var_solution( solc );
     break;
     }

 }  // end( BendersDecompositionSolver::get_var_solution )

/*--------------------------------------------------------------------------*/

bool BendersDecompositionSolver::has_dual_solution( void )
{
 return( ( f_regime == eConvexMaster ) && f_master_solver &&
         f_master_solver->has_dual_solution() );

 }  // end( BendersDecompositionSolver::has_dual_solution )

/*--------------------------------------------------------------------------*/

void BendersDecompositionSolver::get_dual_solution( Configuration * solc )
{
 if( ! has_dual_solution() )
  throw( std::logic_error( "BendersDecompositionSolver::get_dual_solution: "
                           "no dual solution is available" ) );

 f_master_solver->get_dual_solution( solc );

 }  // end( BendersDecompositionSolver::get_dual_solution )

/*--------------------------------------------------------------------------*/

long BendersDecompositionSolver::get_elapsed_iterations( void ) const
{
 return( f_master_solver ? f_master_solver->get_elapsed_iterations() : 0 );

 }  // end( BendersDecompositionSolver::get_elapsed_iterations )

/*--------------------------------------------------------------------------*/

double BendersDecompositionSolver::get_elapsed_time( void ) const
{
 return( f_master_solver ? f_master_solver->get_elapsed_time() : 0 );

 }  // end( BendersDecompositionSolver::get_elapsed_time )

/*--------------------------------------------------------------------------*/
/*------------------------- HANDLING PARAMETERS ----------------------------*/
/*--------------------------------------------------------------------------*/

Solver::idx_type BendersDecompositionSolver::get_num_int_par( void ) const
{
 // TODO: add the master Solver int parameters once f_master_solver exists
 return( intLastBDSlvPar );
 }

/*--------------------------------------------------------------------------*/

Solver::idx_type BendersDecompositionSolver::get_num_str_par( void ) const
{
 // TODO: add the master Solver str parameters once f_master_solver exists
 return( strLastBDSlvPar );
 }

/*--------------------------------------------------------------------------*/

int BendersDecompositionSolver::get_dflt_int_par( idx_type par ) const
{
 switch( par ) {
  case( int_BDSlv_iBCopy ):    return( 0 );
  case( int_BDSlv_Regime ):    return( eConvexMaster );
  case( int_BDSlv_CutType ):   return( eMultiCut );
  case( int_BDSlv_FeasCut ):   return( 1 );
  case( int_BDSlv_MaxRounds ): return( Inf< int >() );
  default:                     return( CDASolver::get_dflt_int_par( par ) );
  }
 }

/*--------------------------------------------------------------------------*/

const std::string & BendersDecompositionSolver::get_dflt_str_par(
						      idx_type par ) const
{
 static const std::string empty;
 if( ( par >= strLastParCDAS ) && ( par < strLastBDSlvPar ) )
  return( empty );
 return( CDASolver::get_dflt_str_par( par ) );
 }

/*--------------------------------------------------------------------------*/

Solver::idx_type BendersDecompositionSolver::int_par_str2idx(
					 const std::string & name ) const
{
 for( idx_type i = 0 ; i < int_pars_BDSlv.size() ; ++i )
  if( name == int_pars_BDSlv[ i ] )
   return( intLastParCDAS + i );
 return( CDASolver::int_par_str2idx( name ) );
 }

/*--------------------------------------------------------------------------*/

Solver::idx_type BendersDecompositionSolver::str_par_str2idx(
					 const std::string & name ) const
{
 for( idx_type i = 0 ; i < str_pars_BDSlv.size() ; ++i )
  if( name == str_pars_BDSlv[ i ] )
   return( strLastParCDAS + i );
 return( CDASolver::str_par_str2idx( name ) );
 }

/*--------------------------------------------------------------------------*/

const std::string & BendersDecompositionSolver::int_par_idx2str(
						      idx_type idx ) const
{
 if( ( idx >= intLastParCDAS ) && ( idx < intLastBDSlvPar ) )
  return( int_pars_BDSlv[ idx - intLastParCDAS ] );
 return( CDASolver::int_par_idx2str( idx ) );
 }

/*--------------------------------------------------------------------------*/

const std::string & BendersDecompositionSolver::str_par_idx2str(
						      idx_type idx ) const
{
 if( ( idx >= strLastParCDAS ) && ( idx < strLastBDSlvPar ) )
  return( str_pars_BDSlv[ idx - strLastParCDAS ] );
 return( CDASolver::str_par_idx2str( idx ) );
 }

/*--------------------------------------------------------------------------*/

void BendersDecompositionSolver::set_par( idx_type par , int value )
{
 switch( par ) {
  case( int_BDSlv_iBCopy ):    f_iBCopy = value;     return;
  case( int_BDSlv_Regime ):    f_regime = value;     return;
  case( int_BDSlv_CutType ):   f_cut_type = value;   return;
  case( int_BDSlv_FeasCut ):   f_feas_cut = value;   return;
  case( int_BDSlv_MaxRounds ): f_max_rounds = value; return;
  default:                     CDASolver::set_par( par , value );
  }
 }

/*--------------------------------------------------------------------------*/

void BendersDecompositionSolver::set_par( idx_type par ,
					  const std::string & value )
{
 switch( par ) {
  case( str_BDSlv_MSName ): f_MSName = value;     return;
  case( str_Bsub_BSCfg ):   f_Bsub_BSCfg = value; return;
  case( str_Mstr_BSCfg ):   f_Mstr_BSCfg = value; return;
  default:                  CDASolver::set_par( par , value );
  }
 }

/*--------------------------------------------------------------------------*/

int BendersDecompositionSolver::get_int_par( idx_type par ) const
{
 switch( par ) {
  case( int_BDSlv_iBCopy ):    return( f_iBCopy );
  case( int_BDSlv_Regime ):    return( f_regime );
  case( int_BDSlv_CutType ):   return( f_cut_type );
  case( int_BDSlv_FeasCut ):   return( f_feas_cut );
  case( int_BDSlv_MaxRounds ): return( f_max_rounds );
  default:                     return( CDASolver::get_int_par( par ) );
  }
 }

/*--------------------------------------------------------------------------*/

const std::string & BendersDecompositionSolver::get_str_par(
						      idx_type par ) const
{
 switch( par ) {
  case( str_BDSlv_MSName ): return( f_MSName );
  case( str_Bsub_BSCfg ):   return( f_Bsub_BSCfg );
  case( str_Mstr_BSCfg ):   return( f_Mstr_BSCfg );
  default:                  return( CDASolver::get_str_par( par ) );
  }
 }

/*--------------------------------------------------------------------------*/
/*------------------------- PROTECTED METHODS ------------------------------*/
/*--------------------------------------------------------------------------*/

void BendersDecompositionSolver::reformulate( void )
{
 static const std::string _prfx = "BendersDecompositionSolver::reformulate: ";

 if( f_reformulated )  // the reformulation is done already
  return;

 if( ! f_Block )
  throw( std::logic_error( _prfx + "no Block is set" ) );

 /* The complicating Variable x are the ColVariable of the root Block, in the
  * order in which it exposes them: the root is the master, and everything
  * that is not in a sub-Block is first-stage by definition. */

 v_x.clear();
 x_index.clear();

 auto take = [ this ]( ColVariable & var ) {
  x_index[ & var ] = v_x.size();
  v_x.push_back( & var );
  };

 for( const auto & el : f_Block->get_static_variables() )
  un_any_const_static( el , take , un_any_type< ColVariable >() );

 for( const auto & el : f_Block->get_dynamic_variables() )
  un_any_const_dynamic( el , take , un_any_type< ColVariable >() );

 if( v_x.empty() )
  throw( std::logic_error( _prfx + "no complicating Variable in the Block" ) );

 // one BendersBFunction per sub-Block - - - - - - - - - - - - - - - - - - - -

 const Index K = f_Block->get_number_nested_Blocks();
 if( ! K )
  throw( std::logic_error( _prfx + "the Block has no sub-Block, hence "
                           "nothing to project out" ) );

 v_BF.assign( K , nullptr );

 for( Index k = 0 ; k < K ; ++k )
  build_BendersBFunction( k );

 /* Each subproblem needs a Solver of its own, since evaluating the value
  * function means solving it: which one is a Configuration matter, so that
  * heterogeneous subproblems can be dealt with. */

 if( ! f_Bsub_BSCfg.empty() )
  for( auto bf : v_BF )
   apply_BSCfg( bf->get_inner_block() , f_Bsub_BSCfg );

 // assemble the master- - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( f_regime == eConvexMaster )
  build_convex_master();
 else
  build_MILP_master();

 f_reformulated = true;

 }  // end( BendersDecompositionSolver::reformulate )

/*--------------------------------------------------------------------------*/

void BendersDecompositionSolver::build_BendersBFunction( Index k )
{
 static const std::string _prfx =
                  "BendersDecompositionSolver::build_BendersBFunction: ";

 auto sub = f_Block->get_nested_Block( k );

 BendersBFunction::MultiVector A;
 BendersBFunction::RealVector b;
 BendersBFunction::ConstraintVector cns;
 BendersBFunction::ConstraintSideVector sides;

 /* The coupling of a subproblem to the master is whatever term in x appears
  * in its Constraint: each such Constraint gives one row of the affine
  * mapping, the x terms are stripped out of it and the mapping puts them
  * back into its side(s) at each evaluation of the value function.
  *
  * BendersBFunction *replaces* the side with A x + b, hence A is the
  * *opposite* of the coefficient the Constraint had, and b is the side it
  * has now: F^k x + E^k( y ) <= u becomes E^k( y ) <= u - F^k x. */

 auto scan = [ & ]( FRowConstraint & con ) {
  auto lf = dynamic_cast< LinearFunction * >( con.get_function() );
  if( ! lf )
   return;

  Subset nms;             // the positions of the x terms in the Function
  BendersBFunction::RealVector row;

  Index i = 0;
  for( auto & cp : lf->get_v_var() ) {
   auto xit = x_index.find( static_cast< const ColVariable * >( cp.first ) );
   if( xit != x_index.end() ) {
    if( row.empty() )
     row.assign( v_x.size() , 0 );
    row[ xit->second ] -= cp.second;
    nms.push_back( i );
    }
   ++i;
   }

  if( row.empty() )   // no x in this Constraint, it is not a coupling one
   return;

  const bool lhs = ( con.get_lhs() > - Inf< double >() );
  const bool rhs = ( con.get_rhs() < Inf< double >() );

  if( ! ( lhs || rhs ) )   // a Constraint with no side at all
   return;

  A.push_back( std::move( row ) );
  b.push_back( lhs ? con.get_lhs() : con.get_rhs() );
  cns.push_back( & con );
  sides.push_back( ( lhs && rhs ) ? BendersBFunction::eBoth
                                  : ( lhs ? BendersBFunction::eLHS
                                          : BendersBFunction::eRHS ) );

  /* The x terms leave the Constraint, which is what makes the subproblem a
   * problem in y alone. The Modification has to be issued, for otherwise the
   * Variable keep the Constraint among the "active" ones and whoever loads
   * the subproblem later finds them there. */

  lf->remove_variables( std::move( nms ) , true );
  };

 for( const auto & el : sub->get_static_constraints() )
  un_any_const_static( el , scan , un_any_type< FRowConstraint >() );

 for( const auto & el : sub->get_dynamic_constraints() )
  un_any_const_dynamic( el , scan , un_any_type< FRowConstraint >() );

 if( cns.empty() )
  throw( std::logic_error( _prfx + "sub-Block " + std::to_string( k ) +
                           " is not coupled to the master by any linear "
                           "Constraint" ) );

 v_BF[ k ] = new BendersBFunction(
  sub , BendersBFunction::VarVector( v_x ) , std::move( A ) , std::move( b ) ,
  std::move( cns ) , std::move( sides ) , nullptr );

 }  // end( BendersDecompositionSolver::build_BendersBFunction )

/*--------------------------------------------------------------------------*/

void BendersDecompositionSolver::build_convex_master( void )
{
 static const std::string _prfx =
                       "BendersDecompositionSolver::build_convex_master: ";

 /* The master is (B) itself: it already has the complicating Variable, the
  * first-stage Constraint and the linear part d( x ) of the Objective, which
  * is exactly what it means for it to be the master. What changes is what
  * hangs below it: the subproblems, which now live inside the
  * BendersBFunction, are replaced by one sub-Block per value function, each
  * with no Variable and no Constraint and an Objective that *is* the value
  * function. This is the structure a bundle-type Solver expects of a
  * sum-function: a linear term at the root plus one C05Function component
  * per sub-Block. */

 f_master = dynamic_cast< AbstractBlock * >( f_Block );
 if( ! f_master )
  throw( std::logic_error( _prfx + "the convex master has to be assembled "
                           "into the Block, which therefore has to be an "
                           "AbstractBlock" ) );

 auto & nested = f_master->access_nested_Blocks();
 nested.clear();

 const bool mx = f_master->get_objective() &&
                 ( f_master->get_objective()->get_sense() == Objective::eMax );

 for( auto bf : v_BF ) {
  auto wrap = new AbstractBlock( f_master );
  auto obj = new FRealObjective( wrap , bf );
  obj->set_sense( mx ? Objective::eMax : Objective::eMin , eNoMod );
  wrap->set_objective( obj , eNoMod );
  f_master->add_nested_Block( wrap );
  }

 }  // end( BendersDecompositionSolver::build_convex_master )

/*--------------------------------------------------------------------------*/

void BendersDecompositionSolver::build_MILP_master( void )
{
 static const std::string _prfx =
                         "BendersDecompositionSolver::build_MILP_master: ";

 /* The master is (B) itself, as in the convex regime, save that the value
  * functions are not handed to the master Solver as Objective: they are
  * inner-approximated by the Benders cuts, which are dynamic Constraint on
  * an epigraph Variable each, and it is this Solver that adds them. */

 f_master = dynamic_cast< AbstractBlock * >( f_Block );
 if( ! f_master )
  throw( std::logic_error( _prfx + "the MILP master has to be assembled into "
                           "the Block, which therefore has to be an "
                           "AbstractBlock" ) );

 auto obj = dynamic_cast< FRealObjective * >( f_master->get_objective() );
 if( ! obj )
  throw( std::logic_error( _prfx + "the Objective of the master is not a "
                           "FRealObjective" ) );

 if( obj->get_sense() != Objective::eMin )
  throw( std::logic_error( _prfx + "only a minimising master is supported, "
                           "the epigraph Variable of a maximising one taking "
                           "an initial bound that is problem-dependent" ) );

 /* The Objective of the master is not touched save for the epigraph
  * Variable, which enter it linearly: it can therefore be linear, as it is
  * in the combinatorial problems the method is classically applied to, or
  * quadratic, as it is whenever the master carries a regularisation term. */

 auto lf = dynamic_cast< LinearFunction * >( obj->get_function() );
 auto qf = dynamic_cast< DQuadFunction * >( obj->get_function() );

 if( ! ( lf || qf ) )
  throw( std::logic_error( _prfx + "the Objective of the master is neither "
                           "linear nor quadratic separable" ) );

 auto & nested = f_master->access_nested_Blocks();
 nested.clear();

 /* One epigraph Variable per subproblem, or a single one for all of them if
  * the cuts are aggregated. They are bounded below by zero, which is what
  * keeps the very first master bounded: this takes the value functions to be
  * non-negative, which is the case of the transportation-like subproblems
  * the method is meant for, and is checked at each round, since a cut is
  * generated whenever the epigraph Variable is below the value function. */

 const Index neta = ( f_cut_type == eSingleCut ) ? 1 : v_BF.size();

 v_eta = new std::vector< ColVariable >( neta );
 for( auto & eta : *v_eta ) {
  eta.is_positive( true );
  eta.set_value( 0 );
  }

 f_master->add_static_variable( *v_eta , "eta" );

 for( auto & eta : *v_eta )
  if( lf )
   lf->add_variable( & eta , 1 , eNoMod );
  else
   qf->add_variable( & eta , 1 , 0 , eNoMod );

 v_cuts = new std::list< FRowConstraint >;
 f_master->add_dynamic_constraint( *v_cuts , "cuts" );

 }  // end( BendersDecompositionSolver::build_MILP_master )

/*--------------------------------------------------------------------------*/

int BendersDecompositionSolver::solve_MILP_master( void )
{
 static const std::string _prfx =
                        "BendersDecompositionSolver::solve_MILP_master: ";

 const Index nx = v_x.size();
 const Index K = v_BF.size();
 const double tol = 1e-9;

 /* A cut is a Constraint on the master, i.e., eta - g x >= alpha for an
  * optimality one and - g x >= alpha for a feasibility one, the latter
  * having no epigraph Variable since it cuts away an x for which the
  * subproblem has no solution at all. */

 auto add_cut = [ & ]( ColVariable * eta ,
                       const std::vector< double > & g , double alpha ) {
  LinearFunction::v_coeff_pair cp;
  cp.reserve( nx + 1 );

  if( eta )
   cp.emplace_back( eta , 1 );

  for( Index i = 0 ; i < nx ; ++i )
   if( g[ i ] )
    cp.emplace_back( v_x[ i ] , - g[ i ] );

  std::list< FRowConstraint > nc( 1 );
  nc.front().set_function( new LinearFunction( std::move( cp ) ) );
  nc.front().set_lhs( alpha );
  nc.front().set_rhs( Inf< double >() );

  /* The Modification has to reach the master Solver, which is in the middle
   * of the loop this cut is generated by: it has to see the new Constraint
   * before it is asked to solve again. */

  f_master->add_dynamic_constraints( *v_cuts , nc , eModBlck );
  };

 int status = kOK;

 for( int round = 0 ; round < f_max_rounds ; ++round ) {

  status = f_master_solver->compute( round > 0 );

  if( ( status != kOK ) && ( status != kLowPrecision ) )
   return( status );

  f_master_solver->get_var_solution();

  f_value = f_master_solver->get_lb();

  // evaluate each value function at the master solution- - - - - - - - - - -

  std::vector< double > gs( nx , 0 );   // the aggregated linearization
  double as = 0;
  double fs = 0;
  bool all_feasible = true;
  Index added = 0;

  for( Index k = 0 ; k < K ; ++k ) {
   auto bf = v_BF[ k ];

   const int st = bf->compute();

   if( ( st != kOK ) && ( st != kInfeasible ) )
    return( st );

   const bool diagonal = ( st == kOK );

   if( ! bf->has_linearization( diagonal ) )
    if( ! bf->compute_new_linearization( diagonal ) ) {
     if( diagonal )
      throw( std::logic_error( _prfx + "no linearization of subproblem " +
                               std::to_string( k ) ) );

     throw( std::logic_error( _prfx + "subproblem " + std::to_string( k ) +
                              " is infeasible and its Solver gives no "
                              "unbounded dual direction, hence no feasibility "
                              "cut can be generated: the certificate only "
                              "exists if the infeasibility is proved by the "
                              "simplex, so ask the subproblem Solver for it "
                              "and switch its presolve off" ) );
     }

   std::vector< double > g( nx , 0 );
   bf->get_linearization_coefficients( g.data() , Range( 0 , nx ) );
   const double alpha = bf->get_linearization_constant();

   if( ! diagonal ) {   // a feasibility cut is never aggregated
    all_feasible = false;
    add_cut( nullptr , g , alpha );
    ++added;
    continue;
    }

   const double fk = bf->get_value();

   if( f_cut_type == eSingleCut ) {
    for( Index i = 0 ; i < nx ; ++i )
     gs[ i ] += g[ i ];
    as += alpha;
    fs += fk;
    continue;
    }

   // the cut is violated only if the epigraph Variable is below the value
   if( fk - (*v_eta)[ k ].get_value() > tol ) {
    add_cut( & (*v_eta)[ k ] , g , alpha );
    ++added;
    }
   }

  /* The aggregated cut is a valid one only if every subproblem has a value
   * to contribute to it, i.e., if none of them is infeasible. */

  if( ( f_cut_type == eSingleCut ) && all_feasible )
   if( fs - (*v_eta)[ 0 ].get_value() > tol ) {
    add_cut( & (*v_eta)[ 0 ] , gs , as );
    ++added;
    }

  if( ! added ) {   // no cut is violated: the master is the problem
   f_solved = true;
   return( status );
   }
  }

 return( kStopIter );

 }  // end( BendersDecompositionSolver::solve_MILP_master )

/*--------------------------------------------------------------------------*/

void BendersDecompositionSolver::map_back_solution( void )
{
 /* Nothing has to be moved around: the master *is* (B), hence the optimal x
  * are written by the master Solver into the very Variable of (B), and the
  * y^k are written by the Solver of each subproblem into the Variable of the
  * inner Block, which are the ones of the original sub-Block. */

 }  // end( BendersDecompositionSolver::map_back_solution )

/*--------------------------------------------------------------------------*/

void BendersDecompositionSolver::acquire_master_solver( void )
{
 static const std::string _prfx =
                    "BendersDecompositionSolver::acquire_master_solver: ";

 if( f_master_solver )  // it is there already
  return;

 if( f_Mstr_BSCfg.empty() )
  throw( std::logic_error( _prfx + "no BlockSolverConfig for the master" ) );

 auto cfg = Configuration::deserialize( f_Mstr_BSCfg );
 auto bsc = dynamic_cast< BlockSolverConfig * >( cfg );
 if( ! bsc ) {
  delete cfg;
  throw( std::invalid_argument( _prfx + f_Mstr_BSCfg +
                                " is not a BlockSolverConfig" ) );
  }

 if( bsc->get_SolverNames().empty() ) {
  delete bsc;
  throw( std::invalid_argument( _prfx + f_Mstr_BSCfg + " names no Solver" ) );
  }

 auto slvr = Solver::new_Solver( bsc->get_SolverName( 0 ) );
 if( ! slvr ) {
  const auto nm = bsc->get_SolverName( 0 );
  delete bsc;
  throw( std::invalid_argument( _prfx + nm + " is not in the Solver "
                                "factory" ) );
  }

 f_master_solver = dynamic_cast< CDASolver * >( slvr );
 if( ! f_master_solver ) {
  delete slvr;
  delete bsc;
  throw( std::invalid_argument( _prfx + "the master Solver is not a "
                                "CDASolver" ) );
  }

 if( auto cc = bsc->get_SolverConfig( 0 ) )
  f_master_solver->set_ComputeConfig( cc );

 delete bsc;

 f_master->register_Solver( f_master_solver );

 }  // end( BendersDecompositionSolver::acquire_master_solver )

/*--------------------------------------------------------------------------*/

void BendersDecompositionSolver::apply_BSCfg( Block * block ,
                                              const std::string & fn )
{
 auto cfg = Configuration::deserialize( fn );

 if( auto bsc = dynamic_cast< BlockSolverConfig * >( cfg ) ) {
  bsc->apply( block );
  bsc->clear();
  delete bsc;
  return;
  }

 /* A "meta-configuration" maps the classname() of a Block to the
  * BlockSolverConfig for it: it is dispatched over the whole sub-tree, so
  * that subproblems of different types get each the Solver that fits it. */

 auto meta = dynamic_cast< SimpleConfiguration<
                    std::map< std::string , Configuration * > > * >( cfg );
 if( ! meta ) {
  delete cfg;
  throw( std::invalid_argument( "BendersDecompositionSolver::apply_BSCfg: " +
                                fn + " is neither a BlockSolverConfig nor a "
                                "map of them" ) );
  }

 std::function< void( Block * ) > dispatch = [ & ]( Block * blk ) {
  auto it = meta->f_value.find( blk->classname() );
  if( it != meta->f_value.end() )
   if( auto bsc = dynamic_cast< BlockSolverConfig * >( it->second ) )
    bsc->apply( blk );

  for( auto sb : blk->get_nested_Blocks() )
   dispatch( sb );
  };

 dispatch( block );

 delete meta;

 }  // end( BendersDecompositionSolver::apply_BSCfg )

/*--------------------------------------------------------------------------*/

Solver::idx_type BendersDecompositionSolver::int_par_ms( idx_type par ) const
{
 return( par - intLastParCDAS + intLastBDSlvPar );
 }

/*--------------------------------------------------------------------------*/

Solver::idx_type BendersDecompositionSolver::str_par_ms( idx_type par ) const
{
 return( par - strLastParCDAS + strLastBDSlvPar );
 }

/*--------------------------------------------------------------------------*/
/*-------------- End File BendersDecompositionSolver.cpp -------------------*/
/*--------------------------------------------------------------------------*/
