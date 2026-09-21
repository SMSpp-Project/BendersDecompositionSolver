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

#include <functional>

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
 "int_BDSlv_CutNorm" ,
 "int_BDSlv_PhaseOneWeights" ,
 "int_BDSlv_Unified" ,
 "int_BDSlv_Restore"
 };

static const std::vector< std::string > dbl_pars_BDSlv = {
 "dbl_BDSlv_CoreMove" ,
 "dbl_BDSlv_EpiWeight"
 };

static const std::vector< std::string > vint_pars_BDSlv = {
 "vintMasterBlock" ,
 "vintMasterVars"
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

 /* Reformulating here is reformulating once and for all, which is what is
  * wanted unless (B) has to be a problem of its own in between two calls to
  * compute() [see block_handling_type]. */

 if( block && ( f_restore == eKeepReformulation ) )
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

 /* Projecting the y^k out has taken the x out of the Constraint of the
  * subproblem, and that is a change to (B), not to anything of this Solver's
  * own: they are put back here, with the coefficient they had, which the
  * mapping carries the opposite of, and with the sides they had, which the
  * mapping has been overwriting at every evaluation. (B) is thus the problem
  * it was before this Solver saw it, which is what lets anybody else look at
  * it [see block_handling_type]. */

 for( Index k = 0 ; k < v_BF.size() ; ++k ) {
  auto bf = v_BF[ k ];
  if( ( ! bf ) || ( k >= v_sides0.size() ) )
   continue;

  const auto & A = bf->get_A();
  const auto & cns = bf->get_constraints();
  const auto & sd0 = v_sides0[ k ];

  for( Index j = 0 ; ( j < cns.size() ) && ( j < sd0.size() ) ; ++j ) {
   auto cn = dynamic_cast< FRowConstraint * >( cns[ j ] );
   if( ! cn )
    continue;

   if( auto lf = dynamic_cast< LinearFunction * >( cn->get_function() ) ) {
    LinearFunction::v_coeff_pair cp;
    for( Index i = 0 ; ( i < A[ j ].size() ) && ( i < v_x.size() ) ; ++i )
     if( A[ j ][ i ] )
      cp.emplace_back( v_x[ i ] , - A[ j ][ i ] );

    if( ! cp.empty() ) {
     lf->add_variables( LinearFunction::v_coeff_pair( cp ) , eNoMod );

     for( auto & pr : cp )
      pr.first->add_active( cn );
     }
    }

   cn->set_lhs( sd0[ j ].first , eNoMod );
   cn->set_rhs( sd0[ j ].second , eNoMod );
   }
  }

 v_sides0.clear();

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
   * un-grafted from it first, or it would be deleted with it.
   *
   * The cuts are Constraint over the x, which are Variable of (B) and
   * survive the master: each of them is emptied here, which is what takes it
   * out of the list of the Constraint the x are active in, since a Variable
   * that kept a cut of a master that is no longer there would hand it to
   * whoever reads (B) next. */

  if( v_cuts ) {
   for( auto & cut : *v_cuts )
    cut.set_function( nullptr , eNoMod );

   v_cuts->clear();
   }

  for( auto bf : v_BF )
   delete bf;

  f_master->access_nested_Blocks().clear();
  if( f_Block )
   f_Block->set_f_Block( f_Block_father );
  delete f_master;
  v_eta = nullptr;
  v_cuts = nullptr;
  }

 /* The phase-one replicas are this Solver's own, and each of them is owned by
  * the BendersBFunction that took it in: deleting those disposes of both. */

 for( auto bf : v_BF1 )
  delete bf;

 v_BF1.clear();
 v_phase1.clear();

 v_BF.clear();
 v_sub.clear();
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
 /* The reformulation is done by set_Block(), with the parameters the Solver
  * has at that moment: a BlockSolverConfig sets them through its
  * ComputeConfig before it registers the Solver, while a parameter that
  * shapes the reformulation and is changed afterwards does not redo it.
  * Asking for it here does nothing when it is done, and complains when
  * there is no Block. */

 reformulate();

 acquire_master_solver();

 f_solved = false;

 int status;

 if( f_regime == eMILPMaster ) {
  status = solve_MILP_master();
  f_ub = f_solved ? f_value : Inf< OFValue >();
  }
 else {
  /* In the convex regime the master Solver is a bundle-type one, which drives
   * the cutting-plane loop by itself: the linearizations it asks the value
   * functions for *are* the Benders cuts. */

  status = f_master_solver->compute( changedvars );

  f_value = f_master_solver->get_lb();
  f_ub = f_master_solver->get_ub();
  f_solved = f_master_solver->has_var_solution();
  }

 map_back_solution();

 /* Everything this Solver has assembled around (B) is disposed of, (B) being
  * asked to be a problem of its own out of compute(); what has been read out
  * of the master is kept, so that the value and the solution can still be
  * asked for. */

 if( f_restore == eRestoreBlock )
  dismantle();

 return( status );

 }  // end( BendersDecompositionSolver::compute )

/*--------------------------------------------------------------------------*/
/*------------------------ READING THE SOLUTION ----------------------------*/
/*--------------------------------------------------------------------------*/

Solver::OFValue BendersDecompositionSolver::get_lb( void )
{
 if( ( f_regime == eMILPMaster ) || ( ! f_master_solver ) )
  return( f_value );

 return( f_master_solver->get_lb() );

 }  // end( BendersDecompositionSolver::get_lb )

/*--------------------------------------------------------------------------*/

Solver::OFValue BendersDecompositionSolver::get_ub( void )
{
 if( ( f_regime == eMILPMaster ) || ( ! f_master_solver ) )
  return( f_ub );

 return( f_master_solver->get_ub() );

 }  // end( BendersDecompositionSolver::get_ub )

/*--------------------------------------------------------------------------*/

bool BendersDecompositionSolver::has_var_solution( void )
{
 if( ( f_regime == eMILPMaster ) || ( ! f_master_solver ) )
  return( f_solved );

 return( f_master_solver->has_var_solution() );

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

Solver::idx_type BendersDecompositionSolver::get_num_vint_par( void ) const
{
 return( vintLastBDSlvPar );
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
  case( int_BDSlv_PhaseOneWeights ): return( eUnitWeights );
  case( int_BDSlv_Unified ):   return( eNoUnified );
  case( int_BDSlv_Restore ):   return( eKeepReformulation );
  default:                     return( CDASolver::get_dflt_int_par( par ) );
  }
 }

/*--------------------------------------------------------------------------*/

double BendersDecompositionSolver::get_dflt_dbl_par( idx_type par ) const
{
 switch( par ) {
  case( dbl_BDSlv_CoreMove ): return( 0.5 );
  case( dbl_BDSlv_EpiWeight ): return( 1 );
  default:                    return( CDASolver::get_dflt_dbl_par( par ) );
  }
 }

/*--------------------------------------------------------------------------*/

const std::vector< int > & BendersDecompositionSolver::get_dflt_vint_par(
						      idx_type par ) const
{
 // both empty: the root is the master, its sub-Block are the subproblems and
 // its ColVariable are the complicating ones
 static const std::vector< int > empty;
 if( ( par >= vintLastParCDAS ) && ( par < vintLastBDSlvPar ) )
  return( empty );
 return( CDASolver::get_dflt_vint_par( par ) );
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

Solver::idx_type BendersDecompositionSolver::vint_par_str2idx(
					 const std::string & name ) const
{
 for( idx_type i = 0 ; i < vint_pars_BDSlv.size() ; ++i )
  if( name == vint_pars_BDSlv[ i ] )
   return( vintLastParCDAS + i );
 return( CDASolver::vint_par_str2idx( name ) );
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

const std::string & BendersDecompositionSolver::vint_par_idx2str(
						      idx_type idx ) const
{
 if( ( idx >= vintLastParCDAS ) && ( idx < vintLastBDSlvPar ) )
  return( vint_pars_BDSlv[ idx - vintLastParCDAS ] );
 return( CDASolver::vint_par_idx2str( idx ) );
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
  case( int_BDSlv_PhaseOneWeights ): f_p1_weights = value; return;
  case( int_BDSlv_Unified ):   f_unified = value;    return;
  case( int_BDSlv_Restore ):   f_restore = value;    return;
  default:                     CDASolver::set_par( par , value );
  }
 }

/*--------------------------------------------------------------------------*/

void BendersDecompositionSolver::set_par( idx_type par , double value )
{
 switch( par ) {
  case( dbl_BDSlv_CoreMove ): f_core_move = value; return;
  case( dbl_BDSlv_EpiWeight ): f_epi_weight = value; return;
  default:                    CDASolver::set_par( par , value );
  }
 }

/*--------------------------------------------------------------------------*/

void BendersDecompositionSolver::set_par( idx_type par ,
					  std::vector< int > && value )
{
 switch( par ) {
  case( vintMasterBlock ): v_master_block = std::move( value );  return;
  case( vintMasterVars ):  v_master_vars = std::move( value );   return;
  default:                 CDASolver::set_par( par , std::move( value ) );
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
  case( int_BDSlv_PhaseOneWeights ): return( f_p1_weights );
  case( int_BDSlv_Unified ):   return( f_unified );
  case( int_BDSlv_Restore ):   return( f_restore );
  default:                     return( CDASolver::get_int_par( par ) );
  }
 }

/*--------------------------------------------------------------------------*/

double BendersDecompositionSolver::get_dbl_par( idx_type par ) const
{
 switch( par ) {
  case( dbl_BDSlv_CoreMove ):  return( f_core_move );
  case( dbl_BDSlv_EpiWeight ): return( f_epi_weight );
  default:                     return( CDASolver::get_dbl_par( par ) );
  }
 }

/*--------------------------------------------------------------------------*/

const std::vector< int > & BendersDecompositionSolver::get_vint_par(
						      idx_type par ) const
{
 switch( par ) {
  case( vintMasterBlock ): return( v_master_block );
  case( vintMasterVars ):  return( v_master_vars );
  default:                 return( CDASolver::get_vint_par( par ) );
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

 std::vector< ColVariable * > all;

 auto take = [ & all ]( ColVariable & var ) { all.push_back( & var ); };

 f_Block->for_each_variable_group( [ & take ]( const BaseGroup & group ) {
   group.for_each_as< ColVariable >( take ); } );

 /* Which of them are complicating is a choice, not a property of the Block
  * [see vintMasterVars]: saying nothing means all of them, which is the
  * convention this Solver was written with. */

 if( v_master_vars.empty() )
  v_x = std::move( all );
 else
  for( auto p : v_master_vars ) {
   if( ( p < 0 ) || ( std::size_t( p ) >= all.size() ) )
    throw( std::logic_error( _prfx + "vintMasterVars names the Variable " +
                             std::to_string( p ) + ", which the master does "
                             "not have" ) );
   v_x.push_back( all[ p ] );
   }

 for( Index i = 0 ; i < v_x.size() ; ++i )
  x_index[ v_x[ i ] ] = i;

 if( v_x.empty() )
  throw( std::logic_error( _prfx + "no complicating Variable in the Block" ) );

 // one BendersBFunction per subproblem - - - - - - - - - - - - - - - - - - -

 const Index K = f_Block->get_number_nested_Blocks();
 if( ! K )
  throw( std::logic_error( _prfx + "the Block has no sub-Block, hence "
                           "nothing to project out" ) );

 /* Which sub-Block are subproblems is the other half of the choice [see
  * vintMasterBlock]: the ones the master does not keep for itself. */

 v_sub.clear();
 for( Index k = 0 ; k < K ; ++k )
  if( std::find( v_master_block.begin() , v_master_block.end() , int( k ) )
      == v_master_block.end() )
   v_sub.push_back( k );

 if( v_sub.empty() )
  throw( std::logic_error( _prfx + "vintMasterBlock keeps every sub-Block in "
                           "the master, hence there is nothing to project "
                           "out" ) );

 const Index NS = v_sub.size();

 /* The unified cut is a cut in ( x , eta ), hence it needs the epigraph
  * Variable to exist while the replicas that separate it are built, i.e.,
  * before the master that holds them; and it needs one of them per
  * subproblem, the eta of a subproblem appearing in its own cuts alone. */

 if( f_unified != eNoUnified ) {
  if( f_regime != eMILPMaster )
   throw( std::invalid_argument( _prfx + "the unified cut needs the MILP "
                                 "regime, the convex one having no epigraph "
                                 "Variable to write it on" ) );

  if( f_cut_type != eMultiCut )
   throw( std::invalid_argument( _prfx + "the unified cut needs un-aggregated "
                                 "cuts, each of them carrying the epigraph "
                                 "Variable of its own subproblem" ) );

  v_eta = new std::vector< ColVariable >( NS );
  for( auto & eta : *v_eta ) {
   eta.is_positive( true );
   eta.set_value( 0 );
   }
  }

 v_BF.assign( NS , nullptr );
 v_sides0.assign( NS , {} );

 // the replica is the separation problem of the unified cut, too
 const bool p1 = ( f_feas_cut == ePhaseOne ) || ( f_unified != eNoUnified );

 v_BF1.assign( p1 ? NS : 0 , nullptr );
 v_phase1.assign( p1 ? NS : 0 , nullptr );

 /* Once its Variable are projected out, a subproblem is no longer a part of
  * the master problem, but it is still a sub-Block of (B): it is the master
  * Solver that has to be told not to look at it. */

 f_ignored.clear();

 for( Index k = 0 ; k < NS ; ++k ) {
  f_ignored.insert( f_Block->get_nested_Block( v_sub[ k ] ) );
  build_BendersBFunction( k );
  }

 /* Each subproblem needs a Solver of its own, since evaluating the value
  * function means solving it: which one is a Configuration matter, so that
  * heterogeneous subproblems can be dealt with. */

 if( ! f_Bsub_BSCfg.empty() ) {
  for( auto bf : v_BF )
   apply_BSCfg( bf->get_inner_block() , f_Bsub_BSCfg );

  for( auto bf : v_BF1 )
   if( bf )
    apply_BSCfg( bf->get_inner_block() , f_Bsub_BSCfg );
  }

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

 auto sub = f_Block->get_nested_Block( v_sub[ k ] );

 /* The sign of the linearization of the value function is read off the sense
  * of the Objective of the subproblem, hence a subproblem with no Objective
  * (or one whose sense is not set) silently gives cuts with the wrong sign;
  * it is rejected here rather than at the first evaluation. */

 if( sub->get_objective_sense() == Objective::eUndef )
  throw( std::invalid_argument( _prfx + "sub-Block " +
                                std::to_string( v_sub[ k ] ) +
                                " has no Objective, or its sense is not "
                                "set" ) );

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

 Subset cpl;        // the positions, in scanning order, of the coupling ones
 Index scanned = 0;

 // the sides each coupling Constraint has now, which is what it is given
 // back when the reformulation is undone [see dismantle()]
 auto & sd0 = v_sides0[ k ];
 sd0.clear();

 auto scan = [ & ]( FRowConstraint & con ) {
  const Index pos = scanned++;
  auto lf = dynamic_cast< LinearFunction * >( con.get_function() );
  if( ! lf )
   return;

  Subset nms;             // the positions of the x terms in the Function
  std::vector< ColVariable * > rmv;   // and the x themselves
  BendersBFunction::RealVector row;

  Index i = 0;
  for( auto & cp : lf->get_v_var() ) {
   auto xit = x_index.find( static_cast< const ColVariable * >( cp.first ) );
   if( xit != x_index.end() ) {
    if( row.empty() )
     row.assign( v_x.size() , 0 );
    row[ xit->second ] -= cp.second;
    nms.push_back( i );
    rmv.push_back( v_x[ xit->second ] );
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
  sd0.emplace_back( con.get_lhs() , con.get_rhs() );
  cns.push_back( & con );
  cpl.push_back( pos );
  sides.push_back( ( lhs && rhs ) ? BendersBFunction::eBoth
                                  : ( lhs ? BendersBFunction::eLHS
                                          : BendersBFunction::eRHS ) );

  /* The x terms leave the Constraint, which is what makes the subproblem a
   * problem in y alone. No Modification is issued: a Solver of the
   * subproblem does not exist yet, being attached to it once the
   * reformulation is done and reading it as it is then, while a Solver of
   * (B) that did exist would be told of a change that is none of its
   * business, this Solver undoing it before anybody is asked anything [see
   * block_handling_type]. What the Modification would do besides telling,
   * i.e., taking the Constraint out of the ones the Variable is active in,
   * is done here, since a Solver building its model by columns reads a
   * Variable's rows from that registration. */

  lf->remove_variables( std::move( nms ) , true , eNoMod );

  for( auto xi : rmv )
   xi->remove_active( & con );
  };

 sub->for_each_constraint_group( [ & scan ]( const BaseGroup & group ) {
   group.for_each_as< FRowConstraint >( scan ); } );

 if( cns.empty() )
  throw( std::logic_error( _prfx + "sub-Block " +
                           std::to_string( v_sub[ k ] ) +
                           " is not coupled to the master by any linear "
                           "Constraint" ) );

 /* The phase one is built before the mapping is handed over, it being built
  * on the same one. */

 if( ! v_BF1.empty() ) {
  std::vector< int > isides( sides.size() );
  for( Index i = 0 ; i < sides.size() ; ++i )
   isides[ i ] = int( sides[ i ] );
  build_phase_one( k , cpl , A , b , isides );
  }

 v_BF[ k ] = new BendersBFunction(
  sub , BendersBFunction::VarVector( v_x ) , std::move( A ) , std::move( b ) ,
  std::move( cns ) , std::move( sides ) , nullptr );

 }  // end( BendersDecompositionSolver::build_BendersBFunction )

/*--------------------------------------------------------------------------*/

void BendersDecompositionSolver::build_phase_one( Index k , const Subset & cpl ,
                       const std::vector< std::vector< double > > & A ,
                       const std::vector< double > & b ,
                       const std::vector< int > & sides )
{
 static const std::string _prfx =
                        "BendersDecompositionSolver::build_phase_one: ";

 auto sub = f_Block->get_nested_Block( v_sub[ k ] );

 /* The replica of the subproblem is the copy of its abstract representation
  * that any Block can be asked for [see AbstractBlock::mirror()], hence it
  * is the subproblem and not something assembled here out of the pieces of
  * it that this Solver happens to know how to read; the copy also says which
  * of its Constraint is the copy of which, which is how the coupling ones
  * are found below. What the phase one adds to it is one slack per side a
  * coupling Constraint can be violated on, and the Objective that minimizes
  * their sum in place of the one of the subproblem. */

 auto rep = new AbstractBlock();
 rep->mirror( sub );

 /* A Constraint the copy could not reproduce makes the phase one a
  * relaxation of the subproblem, which still cuts no feasible x but yields a
  * weaker cut: it is said rather than left to be found out from the
  * numbers. */

 for( auto & issue : rep->get_mirror_issues() )
  if( f_log )
   *f_log << "BendersDecompositionSolver: the phase one of sub-Block "
          << v_sub[ k ] << " is a relaxation, " << issue << std::endl;

 // the coupling Constraint of the subproblem, in the order cpl gives them
 std::vector< FRowConstraint * > orig( cpl.size() , nullptr );
 Index pos = 0 , nxt = 0;

 auto scan = [ & ]( FRowConstraint & con ) {
  if( ( nxt < cpl.size() ) && ( pos == cpl[ nxt ] ) )
   orig[ nxt++ ] = & con;
  ++pos;
  };

 sub->for_each_constraint_group( [ & scan ]( const BaseGroup & group ) {
   group.for_each_as< FRowConstraint >( scan ); } );

 /* The deepest cut bounds the coefficients of the cut rather than the
  * multipliers, which is a bound on an image of them: what that is, here, is
  * the displacement of ( x , eta ) that the rows are allowed instead of the
  * violation the slacks allow [see unified_cut_type]. Each column of the
  * mapping therefore gets a pair of nonnegative Variable, the positive and
  * the negative part of the displacement along it, and the sum of the pairs
  * is the Objective. */

 /* The normalization of the literature is one equation and not a box, and
  * what that is, here, is one slack for all the rows instead of one per row:
  * the column of that slack is the equation, the weight it enters each row
  * with is the coefficient of that row in it, and its cost is the Objective
  * [see phase_one_weight_type]. */

 const bool deepest = ( f_unified == eDeepest );
 const bool oneslack = ( ! deepest ) && ( f_p1_weights == eStaticBSWeights );
 const Index nx = v_x.size();

 auto sl = new std::vector< ColVariable >( deepest  ? 2 * ( nx + 1 ) :
                                           oneslack ? 1
                                                    : 2 * cpl.size() );
 for( auto & s : *sl )
  s.is_positive( true , eNoMod );

 rep->add_static_variable( *sl , deepest ? "d" : "s" );

 BendersBFunction::ConstraintVector cns;
 BendersBFunction::MultiVector A1;
 BendersBFunction::RealVector b1;
 BendersBFunction::ConstraintSideVector sides1;

 // the costs of the slacks, one per side as the slacks themselves
 std::vector< double > w( sl->size() , 1 );

 for( Index i = 0 ; i < cpl.size() ; ++i ) {
  if( ! orig[ i ] )
   continue;

  auto cp = dynamic_cast< FRowConstraint * >( rep->mirror_of( orig[ i ] ) );
  if( ! cp )
   continue;

  auto lf = dynamic_cast< LinearFunction * >( cp->get_function() );
  if( ! lf )
   continue;

  // the cost of the slacks of this row [see phase_one_weight_type]; the
  // deepest cut has no slacks, hence no cost of them to choose, and the
  // single-slack one has a weight per row in place of a cost per row
  if( ( ! deepest ) && ( f_p1_weights == eRowNormWeights ) ) {
   double n2 = 0;
   for( Index j = 0 ; j < lf->get_num_active_var() ; ++j )
    n2 += lf->get_coefficient( j ) * lf->get_coefficient( j );
   for( auto a : A[ i ] )
    n2 += a * a;
   if( n2 > 0 )
    w[ 2 * i ] = w[ 2 * i + 1 ] = 1 / std::sqrt( n2 );
   }

  const bool lhs = ( sides[ i ] != int( BendersBFunction::eRHS ) );
  const bool rhs = ( sides[ i ] != int( BendersBFunction::eLHS ) );

  /* The slack helps the side it is given to: + on a >=, - on a <=. It is
   * added issuing the Modification, although no Solver is there yet to get
   * it: the Constraint registers itself with a Variable coming in only when
   * it receives it, and a Solver building its model by columns reads a
   * Variable's rows from that registration, so a slack added without it
   * would sit in the row and be left out of the model.
   *
   * The displacement of the deepest cut goes in the same way and through the
   * row of the mapping: moving x by d moves the side of this row by A_i d,
   * which on the other side of it is - A_i d, and both sides move together,
   * a displacement being a displacement and not a violation.
   *
   * The single slack goes in with the weight of the row, the sum of the row
   * of the mapping, taken in absolute value: what the equation asks is that
   * the multipliers weigh one all together, so a negative weight would not
   * be an equation over a simplex. A row that weighs zero does not get it,
   * and is therefore a row that nothing can relax. */
  if( deepest )
   for( Index j = 0 ; j < A[ i ].size() ; ++j ) {
    if( A[ i ][ j ] == 0 )
     continue;
    lf->add_variable( & (*sl)[ 2 * j ]     , - A[ i ][ j ] );
    lf->add_variable( & (*sl)[ 2 * j + 1 ] ,   A[ i ][ j ] );
    }
  else
   if( oneslack ) {
    double om = 0;
    for( auto a : A[ i ] )
     om += a;
    om = std::abs( om );

    if( om > 0 )
     lf->add_variable( & (*sl)[ 0 ] , lhs ? om : - om );
    }
   else {
    if( lhs )
     lf->add_variable( & (*sl)[ 2 * i ] , 1 );
    if( rhs )
     lf->add_variable( & (*sl)[ 2 * i + 1 ] , -1 );
    }

  cp->set_lhs( lhs ? orig[ i ]->get_lhs() : - Inf< double >() , eNoMod );
  cp->set_rhs( rhs ? orig[ i ]->get_rhs() : Inf< double >() , eNoMod );

  cns.push_back( cp );
  A1.push_back( A[ i ] );
  b1.push_back( b[ i ] );
  sides1.push_back( BendersBFunction::ConstraintSide( sides[ i ] ) );
  }

 /* What the unified cut adds is one more row that can be violated, the
  * epigraph inequality c^T y <= eta: its slack is the multiplier pi_0 of the
  * literature, its cost is the normalization of that multiplier [see
  * unified_cut_type and dbl_BDSlv_EpiWeight], and eta reaches the row the
  * way x reaches the coupling ones, i.e., through the mapping, which is what
  * makes the linearization of the phase one a cut in ( x , eta ). */

 /* The Objective of a Block is the sum of its own and of those of the Block
  * it is made of, hence the cost of the subproblem is scattered over the tree
  * of the replica, and the phase one, which has to measure the violation of
  * the coupling Constraint and nothing else, has to get rid of all of them
  * and not of the one of the root alone. They are emptied here, and their
  * sum is collected while they are, the unified cut needing it to write the
  * epigraph inequality. */

 LinearFunction::v_coeff_pair ocp;
 double oct = 0;
 bool nonlinear = false;

 std::function< void( Block * ) > strip = [ & ]( Block * b ) {
  if( auto ob = dynamic_cast< FRealObjective * >( b->get_objective() ) ) {
   if( auto lf = dynamic_cast< LinearFunction * >( ob->get_function() ) ) {
    for( auto & vp : lf->get_v_var() )
     ocp.emplace_back( const_cast< ColVariable * >( vp.first ) , vp.second );
    oct += lf->get_constant_term();
    }
   else
    if( ob->get_function() )
     nonlinear = true;

   ob->set_function( new LinearFunction() , eNoMod );
   }

  for( Index i = 0 ; i < b->get_number_nested_Blocks() ; ++i )
   strip( b->get_nested_Block( i ) );
  };

 strip( rep );

 ColVariable * s0 = nullptr;

 if( f_unified != eNoUnified ) {
  if( sub->get_objective_sense() != Objective::eMin )
   throw( std::invalid_argument( _prfx + "the unified cut needs a minimising "
                                 "sub-Block, the epigraph of a maximising "
                                 "one lying on the other side" ) );

  if( nonlinear )
   throw( std::invalid_argument( _prfx + "the unified cut needs a linear "
                                 "Objective in sub-Block " +
                                 std::to_string( v_sub[ k ] ) ) );

  LinearFunction::v_coeff_pair cp( ocp );

  if( deepest ) {
   // eta reaches this row through the last column of the mapping, so its
   // displacement is the pair of that column
   cp.emplace_back( & (*sl)[ 2 * nx ]     , -1 );
   cp.emplace_back( & (*sl)[ 2 * nx + 1 ] ,  1 );
   }
  else
   if( oneslack )
    // the weight of this row in the equation is the omega_0 of the static
    // cut, which is one
    cp.emplace_back( & (*sl)[ 0 ] , -1 );
   else {
    auto s0v = new std::vector< ColVariable >( 1 );
    s0 = & s0v->front();
    s0->is_positive( true , eNoMod );
    rep->add_static_variable( *s0v , "s0" );

    cp.emplace_back( s0 , -1 );   // the slack helps the <= side
    }

  auto epiv = new std::vector< FRowConstraint >( 1 );
  auto epi = & epiv->front();
  epi->set_function( new LinearFunction( std::move( cp ) ) );
  epi->set_lhs( - Inf< double >() , eNoMod );
  epi->set_rhs( 0 , eNoMod );     // the mapping writes eta - c_0 here
  rep->add_static_constraint( *epiv , "epi" );

  for( auto & row : A1 )   // the coupling rows do not see eta
   row.push_back( 0 );

  BendersBFunction::RealVector erow( v_x.size() + 1 , 0 );
  erow.back() = 1;

  A1.push_back( std::move( erow ) );
  b1.push_back( - oct );
  cns.push_back( epi );
  sides1.push_back( BendersBFunction::eRHS );
  }

 /* The total weighted violation, which is what the phase one minimizes; for
  * the deepest cut the same sum is the l1 norm of the displacement, every
  * pair costing one, hence the distance of ( x , eta ) from the epigraph. */

 auto olf = new LinearFunction();
 for( Index j = 0 ; j < sl->size() ; ++j )
  olf->add_variable( & (*sl)[ j ] , w[ j ] , eNoMod );

 if( s0 )
  olf->add_variable( s0 , f_epi_weight , eNoMod );

 auto obj = new FRealObjective( rep , olf );
 obj->set_sense( Objective::eMin , eNoMod );

 delete rep->get_objective();   // the one the copy took from the subproblem
 rep->set_objective( obj , eNoMod );

 BendersBFunction::VarVector vx1( v_x );
 if( f_unified != eNoUnified )   // the unified cut is a cut in ( x , eta )
  vx1.push_back( & (*v_eta)[ k ] );

 v_phase1[ k ] = rep;
 v_BF1[ k ] = new BendersBFunction( rep , std::move( vx1 ) ,
                                    std::move( A1 ) , std::move( b1 ) ,
                                    std::move( cns ) , std::move( sides1 ) ,
                                    nullptr );

 }  // end( BendersDecompositionSolver::build_phase_one )

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

 if( ! v_eta ) {   // the unified cut has built them already
  v_eta = new std::vector< ColVariable >( neta );
  for( auto & eta : *v_eta ) {
   eta.is_positive( true );
   eta.set_value( 0 );
   }
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

 /* The unified cut is a cut in ( x , eta ): the separation problem is the
  * phase one with the epigraph inequality among the rows that carry a slack
  * [see unified_cut_type], its value is zero exactly at the points of the
  * epigraph of the value function, and its linearization there, asked to be
  * nonpositive, is the cut. Its scale is arbitrary, the cut being homogeneous
  * in the multipliers, hence it is scaled like a feasibility one. */

 auto add_unified_cut = [ & ]( ColVariable * eta ,
                               const std::vector< double > & g ,
                               double alpha ) {
  double scale = 1;
  if( f_cut_norm != eNoNorm ) {
   double nrm = 0;
   for( Index i = 0 ; i <= nx ; ++i )
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

  if( g[ nx ] )
   cp.emplace_back( eta , - g[ nx ] * scale );

  for( Index i = 0 ; i < nx ; ++i )
   if( g[ i ] )
    cp.emplace_back( v_x[ i ] , - g[ i ] * scale );

  std::list< FRowConstraint > nc( 1 );
  nc.front().set_function( new LinearFunction( std::move( cp ) ) );
  nc.front().set_lhs( alpha * scale );
  nc.front().set_rhs( Inf< double >() );

  f_master->add_dynamic_constraints( *v_cuts , nc , eModBlck );
  ++f_cuts;
  };

 /* One round of separation of the unified cuts: it returns how many have been
  * added, zero saying that the incumbent is in the epigraph of every value
  * function and the master is therefore the problem.
  *
  * With at_core the separation is done at the core point instead, which is
  * the Pareto-optimal cut of Papadakos read for a cut that carries the
  * epigraph Variable too [see cut_strengthening_type]: the separation
  * problem of a unified cut has, on the instances measured, an optimal face
  * with more than one vertex, so which supporting half-space comes out
  * depends on the algorithm that solves it, and separating at an interior
  * point is what chooses among them. The cut is added whether or not it is
  * violated, it being generated at a point that has nothing to do with the
  * incumbent. */

 auto separate_unified = [ & ]( int & status , bool at_core = false )
                                                               -> Index {
  Index added = 0;

  std::vector< double > x_inc;
  std::vector< double > eta_inc;

  if( at_core ) {
   x_inc.resize( nx );
   eta_inc.resize( K );
   for( Index i = 0 ; i < nx ; ++i ) {
    x_inc[ i ] = v_x[ i ]->get_value();
    v_x[ i ]->set_value( v_core[ i ] );
    }
   for( Index k = 0 ; k < K ; ++k ) {
    eta_inc[ k ] = (*v_eta)[ k ].get_value();
    (*v_eta)[ k ].set_value( v_core_eta[ k ] );
    }
   }

  for( Index k = 0 ; k < K ; ++k ) {
   auto bf = v_BF1[ k ];

   /* An empty separation problem is an answer and not a failure for the
    * normalizations that do not give every row something of its own to be
    * relaxed by: with the deepest one, what it says is that no displacement
    * of ( x , eta ) makes the subproblem consistent, and with the single
    * slack that the violated row is one the slack does not reach. Either
    * way, no subproblem is feasible for any x, hence (B) is not, and the
    * status travels up as it comes. */

   const int st = bf->compute();
   if( ( st != kOK ) && ( st != kLowPrecision ) ) {
    status = st;
    return( added );
    }

   if( ! bf->has_linearization( true ) )
    if( ! bf->compute_new_linearization( true ) )
     throw( std::logic_error( _prfx + "no linearization of the separation "
                              "problem of subproblem " +
                              std::to_string( k ) ) );

   const double viol = bf->get_value();
   if( ( ! at_core ) && ( viol <= 0 ) )   // the incumbent is in the epigraph
    continue;

   std::vector< double > g( nx + 1 , 0 );
   bf->get_linearization_coefficients( g.data() , Range( 0 , nx + 1 ) );

   /* How much the cut is violated is not the value of the separation
    * problem: that value is scaled by the multipliers, which the costs of
    * the slacks bound [see unified_cut_type], hence dividing it by the
    * coefficient the cut gives the epigraph Variable, or by the largest of
    * the others when that is zero, is what puts the violation back into the
    * units of the master and makes the test independent of those costs. */

   double cs = std::abs( g[ nx ] );
   double xs = std::abs( (*v_eta)[ k ].get_value() );

   if( cs == 0 )
    for( Index i = 0 ; i < nx ; ++i ) {
     cs = std::max( cs , std::abs( g[ i ] ) );
     xs = std::max( xs , std::abs( v_x[ i ]->get_value() ) );
     }

   if( ( ! at_core ) &&
       ( ( cs == 0 ) || ( viol <= tol * cs * std::max( 1.0 , xs ) ) ) )
    continue;

   add_unified_cut( & (*v_eta)[ k ] , g , bf->get_linearization_constant() );
   ++added;
   }

  /* The core point is then moved towards the incumbent, so that it keeps
   * track of where the master is going, and the incumbent is put back where
   * the master solver left it. */

  if( at_core ) {
   for( Index i = 0 ; i < nx ; ++i ) {
    v_core[ i ] += f_core_move * ( x_inc[ i ] - v_core[ i ] );
    v_x[ i ]->set_value( x_inc[ i ] );
    }
   for( Index k = 0 ; k < K ; ++k ) {
    v_core_eta[ k ] += f_core_move * ( eta_inc[ k ] - v_core_eta[ k ] );
    (*v_eta)[ k ].set_value( eta_inc[ k ] );
    }
   }

  return( added );
  };

 /* Cutting away an x at which a subproblem has no solution: which cut that is
  * is a parameter [see feasibility_cut_type]. The no-good one forbids the
  * current assignment and nothing else, hence it is written out of the
  * incumbent rather than out of any certificate. */

 auto add_feasibility_cut = [ & ]( Index kk , const std::vector< double > & g ,
                                   double alpha , bool check = true ) {
  if( f_feas_cut == eAlwaysFeasible )
   throw( std::logic_error( _prfx + "a subproblem is infeasible while "
                            "int_BDSlv_FeasCut says none can be" ) );

  if( f_feas_cut == eFarkas ) {
   add_cut( nullptr , g , alpha );
   return;
   }

  if( f_feas_cut == ePhaseOne ) {
   /* The cut is the linearization of the least violation, asked to be
    * nonpositive. The phase one is a relaxation of the subproblem whenever
    * something of it could not be replicated, so the cut may fail to cut the
    * incumbent: the certificate, which never fails to, is asked for then. */

   auto bf = v_BF1[ kk ];
   std::vector< double > g1( nx , 0 );
   double alpha1 = 0;

   if( bf && ( bf->compute() == kOK ) &&
       ( bf->has_linearization( true ) ||
         bf->compute_new_linearization( true ) ) ) {
    bf->get_linearization_coefficients( g1.data() , Range( 0 , nx ) );
    alpha1 = bf->get_linearization_constant();

    double viol = alpha1;
    for( Index i = 0 ; i < nx ; ++i )
     viol += g1[ i ] * v_x[ i ]->get_value();

    if( ( ! check ) || ( viol > tol ) ) {   // it does cut the incumbent
     add_cut( nullptr , g1 , alpha1 , false );
     return;
     }
    }

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

  /* An infeasible subproblem that is cut away by forbidding the assignment
   * needs no certificate, and asking for one would throw; the phase one does
   * ask for it, the certificate being what it falls back on. */

  if( ( ! diagonal ) && ( f_feas_cut == eCombinatorial ) )
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
    if( ( f_feas_cut == eFarkas ) || ( f_feas_cut == ePhaseOne ) )
     add_feasibility_cut( k , g , alpha , false );
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

  // the core point of a unified cut has an epigraph component too, the cut
  // being a cut in ( x , eta )
  if( f_unified != eNoUnified )
   v_core_eta.assign( v_eta->size() , 0 );
  }

 int status = kOK;

 for( int round = 0 ; round < f_max_rounds ; ++round ) {

  ++f_rounds;

  status = f_master_solver->compute( round > 0 );

  if( ( status != kOK ) && ( status != kLowPrecision ) )
   return( status );

  f_master_solver->get_var_solution();

  f_value = f_master_solver->get_lb();

  /* The unified cuts are separated on their own: one problem per subproblem,
   * telling feasibility and optimality apart by itself, hence neither the
   * feasibility cuts nor the aggregation have anything to do here. The
   * Pareto ones do: the separation problem of a unified cut has an optimal
   * face with more than one vertex, and a round at the core point is what
   * chooses among them. */

  if( f_unified != eNoUnified ) {
   int st = kOK;
   Index added = separate_unified( st );

   if( ( st != kOK ) && ( st != kLowPrecision ) )
    return( st );

   if( f_pareto == ePapadakos ) {
    added += separate_unified( st , true );

    if( ( st != kOK ) && ( st != kLowPrecision ) )
     return( st );
    }

   if( added )
    continue;

   /* The separation problems say that the incumbent is in the epigraph of
    * every value function, but they say it through multipliers that the
    * costs of the slacks bound, hence a cut can be missed when those costs
    * are far from the scale of the model. What closes the loop is therefore
    * the value functions themselves, evaluated here, which is also what
    * leaves the y^k where map_back_solution() reads them from: if one of
    * them is above its epigraph Variable after all, its ordinary cut is
    * added and the loop goes on. */

   Index late = 0;

   for( Index k = 0 ; k < K ; ++k ) {
    std::vector< double > g( nx , 0 );
    double alpha;
    bool diagonal;

    const int st = get_cut( k , g , alpha , diagonal );
    if( st != kOK )
     return( st );

    if( ! diagonal ) {
     add_feasibility_cut( k , g , alpha );
     ++late;
     continue;
     }

    if( violated( v_BF[ k ]->get_value() , (*v_eta)[ k ].get_value() ) ) {
     add_cut( & (*v_eta)[ k ] , g , alpha );
     ++late;
     }
    }

   if( late )
    continue;

   f_solved = true;
   return( status );
   }

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
    add_feasibility_cut( k , g , alpha );
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
