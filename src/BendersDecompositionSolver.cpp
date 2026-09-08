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

#include <algorithm>

#include <cmath>

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
 "int_BDSlv_MaxRounds" ,
 "int_BDSlv_Pareto" ,
 "int_BDSlv_CutNorm"
 };

static const std::vector< std::string > dbl_pars_BDSlv = {
 "dbl_BDSlv_CoreMove"
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
 dismantle();

 }  // end( BendersDecompositionSolver::~BendersDecompositionSolver )

/*--------------------------------------------------------------------------*/

void BendersDecompositionSolver::set_Block( Block * block )
{
 if( block == f_Block )  // nothing to do
  return;

 dismantle();   // whatever was assembled around the previous (B) goes

 CDASolver::set_Block( block );

 if( block )
  reformulate();
 }

/*--------------------------------------------------------------------------*/

void BendersDecompositionSolver::dismantle( void )
{
 if( ! f_reformulated )   // there is nothing to dismantle
  return;

 if( f_master_solver ) {
  if( f_master )
   f_master->unregister_Solver( f_master_solver );
  delete f_master_solver;
  f_master_solver = nullptr;
  }

 /* Each subproblem is owned by (B), which still has it among its sub-Block:
  * the BendersBFunction has to let go of it, or it would be deleted twice.
  * Giving it back means restoring its father, too, a BendersBFunction being
  * a Block itself and having made itself the father when it took it in. */

 for( auto bf : v_BF )
  if( auto sub = bf ? bf->get_inner_block() : nullptr ) {
   bf->set_inner_block( nullptr , false );
   sub->set_f_Block( f_Block );
   }

 if( f_master == f_Block ) {
  /* The value-function sub-Block have been added to (B), and each of them
   * owns, through its Objective, the BendersBFunction it carries: deleting
   * them is deleting those, and (B) is left with the sub-Block it had. */

  auto & nested = f_master->access_nested_Blocks();
  for( auto wrap : v_wrap )
   nested.erase( std::remove( nested.begin() , nested.end() , wrap ) ,
                 nested.end() );

  for( auto wrap : v_wrap )
   delete wrap;

  v_wrap.clear();
  }
 else {
  /* The master is this Solver's own, and so are the BendersBFunction, which
   * nothing else refers to: deleting the master disposes of the epigraph
   * Variable, of the cuts and of the Objective they enter, and (B) has to be
   * un-grafted from it first, or it would be deleted with it. */

  for( auto bf : v_BF )
   delete bf;

  f_master->access_nested_Blocks().clear();
  if( f_Block )
   f_Block->set_f_Block( f_Block_father );
  delete f_master;
  v_eta = nullptr;
  v_cuts = nullptr;
  }

 v_BF.clear();
 f_ignored.clear();
 v_x.clear();
 x_index.clear();
 f_master = nullptr;
 f_Block_father = nullptr;
 f_reformulated = false;
 f_solved = false;

 }  // end( BendersDecompositionSolver::dismantle )

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
 /* The x are written into (B) by the master Solver, which has (B) either as
  * the Block it is attached to or as its only sub-Block; the y^k of each
  * subproblem are written into the inner Block, which holds the very
  * Variable the original sub-Block had. Note that in the MILP regime a
  * Configuration addressing a sub-Block by position refers to the master,
  * where (B) is the first sub-Block. */

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
 if( f_regime == eMILPMaster )
  return( f_rounds );

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

Solver::idx_type BendersDecompositionSolver::get_num_dbl_par( void ) const
{
 // TODO: add the master Solver dbl parameters once f_master_solver exists
 return( dblLastBDSlvPar );
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
  case( int_BDSlv_FeasCut ):   return( eFarkas );
  case( int_BDSlv_MaxRounds ): return( Inf< int >() );
  case( int_BDSlv_Pareto ):    return( eNoPareto );
  case( int_BDSlv_CutNorm ):   return( eNoNorm );
  default:                     return( CDASolver::get_dflt_int_par( par ) );
  }
 }

/*--------------------------------------------------------------------------*/

double BendersDecompositionSolver::get_dflt_dbl_par( idx_type par ) const
{
 switch( par ) {
  case( dbl_BDSlv_CoreMove ): return( 0.5 );
  default:                    return( CDASolver::get_dflt_dbl_par( par ) );
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

Solver::idx_type BendersDecompositionSolver::dbl_par_str2idx(
					 const std::string & name ) const
{
 for( idx_type i = 0 ; i < dbl_pars_BDSlv.size() ; ++i )
  if( name == dbl_pars_BDSlv[ i ] )
   return( dblLastParCDAS + i );
 return( CDASolver::dbl_par_str2idx( name ) );
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

const std::string & BendersDecompositionSolver::dbl_par_idx2str(
						      idx_type idx ) const
{
 if( ( idx >= dblLastParCDAS ) && ( idx < dblLastBDSlvPar ) )
  return( dbl_pars_BDSlv[ idx - dblLastParCDAS ] );
 return( CDASolver::dbl_par_idx2str( idx ) );
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
  case( int_BDSlv_Pareto ):    f_pareto = value;     return;
  case( int_BDSlv_CutNorm ):   f_cut_norm = value;   return;
  default:                     CDASolver::set_par( par , value );
  }
 }

/*--------------------------------------------------------------------------*/

void BendersDecompositionSolver::set_par( idx_type par , double value )
{
 switch( par ) {
  case( dbl_BDSlv_CoreMove ): f_core_move = value; return;
  default:                    CDASolver::set_par( par , value );
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
  case( int_BDSlv_Pareto ):    return( f_pareto );
  case( int_BDSlv_CutNorm ):   return( f_cut_norm );
  default:                     return( CDASolver::get_int_par( par ) );
  }
 }

/*--------------------------------------------------------------------------*/

double BendersDecompositionSolver::get_dbl_par( idx_type par ) const
{
 switch( par ) {
  case( dbl_BDSlv_CoreMove ): return( f_core_move );
  default:                    return( CDASolver::get_dbl_par( par ) );
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

 /* Once its Variable are projected out, a subproblem is no longer a part of
  * the master problem, but it is still a sub-Block of (B): it is the master
  * Solver that has to be told not to look at it. */

 f_ignored.clear();

 for( Index k = 0 ; k < K ; ++k ) {
  f_ignored.insert( f_Block->get_nested_Block( k ) );
  build_BendersBFunction( k );
  }

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
  * per sub-Block.
  *
  * The subproblems are not removed from (B): they are the sub-Block the
  * master Solver is told to ignore [see acquire_master_solver()], so that
  * (B) keeps owning them and the value-function sub-Block are simply added
  * next to them. */

 f_master = dynamic_cast< AbstractBlock * >( f_Block );
 if( ! f_master )
  throw( std::logic_error( _prfx + "the convex master has to be assembled "
                           "into the Block, which therefore has to be an "
                           "AbstractBlock" ) );

 const bool mx = f_master->get_objective() &&
                 ( f_master->get_objective()->get_sense() == Objective::eMax );

 v_wrap.clear();
 v_wrap.reserve( v_BF.size() );

 for( auto bf : v_BF ) {
  auto wrap = new AbstractBlock( f_master );
  auto obj = new FRealObjective( wrap , bf );
  obj->set_sense( mx ? Objective::eMax : Objective::eMin , eNoMod );
  wrap->set_objective( obj , eNoMod );
  f_master->add_nested_Block( wrap );
  v_wrap.push_back( wrap );
  }

 }  // end( BendersDecompositionSolver::build_convex_master )

/*--------------------------------------------------------------------------*/

void BendersDecompositionSolver::build_MILP_master( void )
{
 static const std::string _prfx =
                         "BendersDecompositionSolver::build_MILP_master: ";

 /* What the MILP regime adds to (B) are the epigraph Variable and the cuts,
  * neither of which belongs to (B): the master is therefore a Block of this
  * Solver's own, holding them, into which (B) is grafted as its only
  * sub-Block. The master Solver reads the x, the first-stage Constraint X
  * and the Objective d( x ) out of (B), which is left exactly as it is: its
  * Objective is not touched, hence it can be of whatever kind the master
  * Solver takes, linear as in the combinatorial problems the method is
  * classically applied to, or quadratic as whenever the master carries a
  * regularisation term. */

 if( f_Block->get_objective_sense() != Objective::eMin )
  throw( std::logic_error( _prfx + "only a minimising master is supported, "
                           "the epigraph Variable of a maximising one taking "
                           "an initial bound that is problem-dependent" ) );

 f_master = new AbstractBlock;
 f_Block_father = f_Block->get_f_Block();
 f_master->add_nested_Block( f_Block );

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

 /* The epigraph Variable enter the Objective linearly and with coefficient
  * one: since the Objective of the whole master is the sum of those of the
  * Block it is made of, this is an Objective of the master Block alone, and
  * the one of (B) stays what it was. */

 LinearFunction::v_coeff_pair cp;
 cp.reserve( neta );
 for( auto & eta : *v_eta )
  cp.emplace_back( & eta , 1 );

 auto obj = new FRealObjective( f_master ,
                                new LinearFunction( std::move( cp ) ) );
 obj->set_sense( Objective::eMin , eNoMod );
 f_master->set_objective( obj , eNoMod );

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
 /* How much a cut has to be violated to be worth adding, relative to the
  * value of the function it cuts: an absolute threshold would be below the
  * noise of the master on an instance whose values are large, and the loop
  * would keep separating cuts that say nothing. */

 const double tol = 1e-9;

 auto violated = [ tol ]( double f , double eta ) {
  return( f - eta > tol * std::max( 1.0 , std::abs( f ) ) );
  };

 /* A cut is a Constraint on the master, i.e., eta - g x >= alpha for an
  * optimality one and - g x >= alpha for a feasibility one, the latter
  * having no epigraph Variable since it cuts away an x for which the
  * subproblem has no solution at all. */

 auto add_cut = [ & ]( ColVariable * eta ,
                       const std::vector< double > & g , double alpha ,
                       bool scalable = true ) {
  /* A feasibility cut is a ray, hence its scale is arbitrary and putting all
   * of them on the same one is up to whoever generates them [see
   * feasibility_cut_norm_type]; a combinatorial one is not a ray, and its
   * coefficients mean what they say. */

  double scale = 1;
  if( scalable && ( ! eta ) && ( f_cut_norm != eNoNorm ) ) {
   double nrm = 0;
   for( Index i = 0 ; i < nx ; ++i )
    switch( f_cut_norm ) {
     case( eOneNorm ): nrm += std::abs( g[ i ] ); break;
     case( eTwoNorm ): nrm += g[ i ] * g[ i ]; break;
     default:          nrm = std::max( nrm , std::abs( g[ i ] ) );
     }

   if( f_cut_norm == eTwoNorm )
    nrm = std::sqrt( nrm );

   if( nrm > 0 )
    scale = 1 / nrm;
   }

  LinearFunction::v_coeff_pair cp;
  cp.reserve( nx + 1 );

  if( eta )
   cp.emplace_back( eta , 1 );

  for( Index i = 0 ; i < nx ; ++i )
   if( g[ i ] )
    cp.emplace_back( v_x[ i ] , - g[ i ] * scale );

  std::list< FRowConstraint > nc( 1 );
  nc.front().set_function( new LinearFunction( std::move( cp ) ) );
  nc.front().set_lhs( alpha * scale );
  nc.front().set_rhs( Inf< double >() );

  /* The Modification has to reach the master Solver, which is in the middle
   * of the loop this cut is generated by: it has to see the new Constraint
   * before it is asked to solve again. */

  f_master->add_dynamic_constraints( *v_cuts , nc , eModBlck );
  ++f_cuts;
  };

 /* Cutting away an x at which a subproblem has no solution: which cut that is
  * is a parameter [see feasibility_cut_type]. The no-good one forbids the
  * current assignment and nothing else, hence it is written out of the
  * incumbent rather than out of any certificate. */

 auto add_feasibility_cut = [ & ]( const std::vector< double > & g ,
                                   double alpha ) {
  if( f_feas_cut == eAlwaysFeasible )
   throw( std::logic_error( _prfx + "a subproblem is infeasible while "
                            "int_BDSlv_FeasCut says none can be" ) );

  if( f_feas_cut == eFarkas ) {
   add_cut( nullptr , g , alpha );
   return;
   }

  std::vector< double > ng( nx );
  double rhs = 1;
  for( Index i = 0 ; i < nx ; ++i ) {
   auto xi = v_x[ i ];
   if( ( ! xi->is_integer() ) || ( xi->get_lb() < 0 ) || ( xi->get_ub() > 1 ) )
    throw( std::logic_error( _prfx + "a combinatorial cut needs all the "
                             "complicating Variable to be binary" ) );

   if( xi->get_value() > 0.5 ) {  // the ones of the assignment
    ng[ i ] = 1;
    --rhs;
    }
   else
    ng[ i ] = -1;
   }

  add_cut( nullptr , ng , rhs , false );
  };

 /* Evaluating one value function at the point the x currently hold, and
  * reading the cut out of it: the linearization is the cut, a diagonal one
  * when the subproblem is feasible and a vertical one, i.e., the Farkas
  * certificate, when it is not. */

 auto get_cut = [ & ]( Index k , std::vector< double > & g , double & alpha ,
                       bool & diagonal ) {
  auto bf = v_BF[ k ];

  const int st = bf->compute();

  if( ( st != kOK ) && ( st != kInfeasible ) )
   return( st );

  diagonal = ( st == kOK );

  /* An infeasible subproblem that is cut away with something other than the
   * Farkas certificate does not need one, and asking for it would throw. */

  if( ( ! diagonal ) && ( f_feas_cut != eFarkas ) )
   return( int( kOK ) );

  if( ! bf->has_linearization( diagonal ) )
   if( ! bf->compute_new_linearization( diagonal ) ) {
    if( diagonal )
     throw( std::logic_error( _prfx + "no linearization of subproblem " +
                              std::to_string( k ) ) );

    throw( std::logic_error( _prfx + "subproblem " + std::to_string( k ) +
                             " is infeasible and its Solver gives no "
                             "unbounded dual direction, hence no feasibility "
                             "cut can be generated: switch the presolve of "
                             "the subproblem Solver off, it detecting the "
                             "infeasibility on the reduced problem and "
                             "leaving no certificate for the original one" ) );
    }

  bf->get_linearization_coefficients( g.data() , Range( 0 , nx ) );
  alpha = bf->get_linearization_constant();
  return( int( kOK ) );
  };

 /* The Pareto-optimal cuts of Papadakos: the very same evaluation, but at the
  * core point rather than at the incumbent, hence one more cut per subproblem
  * and per round. The cut is added whether or not it is violated, it being
  * generated at a point that has nothing to do with the incumbent, and the
  * core point is then moved towards the incumbent so that it keeps track of
  * where the master is going [see cut_strengthening_type]. */

 auto add_pareto_cuts = [ & ]( void ) {
  std::vector< double > x_inc( nx );
  for( Index i = 0 ; i < nx ; ++i ) {
   x_inc[ i ] = v_x[ i ]->get_value();
   v_x[ i ]->set_value( v_core[ i ] );
   }

  std::vector< double > g( nx );
  std::vector< double > gs( nx , 0 );
  double as = 0;
  bool all_feasible = true;

  for( Index k = 0 ; k < K ; ++k ) {
   double alpha;
   bool diagonal;
   if( get_cut( k , g , alpha , diagonal ) != kOK )
    break;

   /* The core point is not the incumbent, hence a no-good cut written out of
    * it would forbid an assignment that has nothing to do with the one being
    * looked at, and would not even be valid: only the certificate, which is a
    * cut of the value function itself, can be taken here. */

   if( ! diagonal ) {   // a feasibility cut is never aggregated
    all_feasible = false;
    if( f_feas_cut == eFarkas )
     add_cut( nullptr , g , alpha );
    continue;
    }

   if( f_cut_type == eSingleCut ) {
    for( Index i = 0 ; i < nx ; ++i )
     gs[ i ] += g[ i ];
    as += alpha;
    continue;
    }

   add_cut( & (*v_eta)[ k ] , g , alpha );
   }

  if( ( f_cut_type == eSingleCut ) && all_feasible )
   add_cut( & (*v_eta)[ 0 ] , gs , as );

  for( Index i = 0 ; i < nx ; ++i ) {
   v_core[ i ] += f_core_move * ( x_inc[ i ] - v_core[ i ] );
   v_x[ i ]->set_value( x_inc[ i ] );
   }
  };

 /* The core point starts wherever the x are when the loop does, which is
  * where whoever built the model left them: a point in the relative interior
  * of the master feasible set is what the theory asks for. */

 f_cuts = 0;
 f_rounds = 0;

 if( f_pareto == ePapadakos ) {
  v_core.resize( nx );
  for( Index i = 0 ; i < nx ; ++i )
   v_core[ i ] = v_x[ i ]->get_value();
  }

 int status = kOK;

 for( int round = 0 ; round < f_max_rounds ; ++round ) {

  ++f_rounds;

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
   std::vector< double > g( nx , 0 );
   double alpha;
   bool diagonal;

   const int st = get_cut( k , g , alpha , diagonal );
   if( st != kOK )
    return( st );

   if( ! diagonal ) {   // a feasibility cut is never aggregated
    all_feasible = false;
    add_feasibility_cut( g , alpha );
    ++added;
    continue;
    }

   const double fk = v_BF[ k ]->get_value();

   if( f_cut_type == eSingleCut ) {
    for( Index i = 0 ; i < nx ; ++i )
     gs[ i ] += g[ i ];
    as += alpha;
    fs += fk;
    continue;
    }

   // the cut is violated only if the epigraph Variable is below the value
   if( violated( fk , (*v_eta)[ k ].get_value() ) ) {
    add_cut( & (*v_eta)[ k ] , g , alpha );
    ++added;
    }
   }

  /* The aggregated cut is a valid one only if every subproblem has a value
   * to contribute to it, i.e., if none of them is infeasible. */

  if( ( f_cut_type == eSingleCut ) && all_feasible )
   if( violated( fs , (*v_eta)[ 0 ].get_value() ) ) {
    add_cut( & (*v_eta)[ 0 ] , gs , as );
    ++added;
    }

  if( ! added ) {   // no cut is violated: the master is the problem
   f_solved = true;
   return( status );
   }

  /* The Pareto-optimal cuts come after the ordinary ones, and only when
   * those exist: a round that adds none is the last one, and the value
   * functions have to be left evaluated at the incumbent, which is where
   * map_back_solution() reads the y^k of each subproblem from. */

  if( f_pareto == ePapadakos )
   add_pareto_cuts();
  }

 return( kStopIter );

 }  // end( BendersDecompositionSolver::solve_MILP_master )

/*--------------------------------------------------------------------------*/

void BendersDecompositionSolver::map_back_solution( void )
{
 /* Nothing has to be moved around: the x of the master *are* the Variable of
  * (B), whether the master is (B) itself or the Block it has been grafted
  * into, and the y^k are written by the Solver of each subproblem into the
  * Variable of the inner Block, which are the ones of the original
  * sub-Block. */

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

 /* The master Solver writes on the same stream as this one, without which
  * the cutting-plane loop is silent about the only thing that can go wrong
  * in it, i.e., the master problem; this happens *before* the ComputeConfig
  * is applied, so that a configuration naming a log of its own has the last
  * word, as a configuration always has. */

 if( f_log )
  f_master_solver->set_log( f_log );

 if( auto cc = bsc->get_SolverConfig( 0 ) )
  f_master_solver->set_ComputeConfig( cc );

 delete bsc;

 /* The exclusion list has to be installed *before* the Solver is attached to
  * the master, for it is at that moment that the Block tree is scanned and
  * the model loaded: the subproblems have to be invisible already. */

 if( ! f_ignored.empty() )
  f_master_solver->set_excluded_blocks( & f_ignored );

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
