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
 // TODO: detach and delete f_master_solver, f_master and the v_BF; revert
 //       any eviction of the original sub-Block of (B)
 }

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
 // TODO: drive the master solver according to f_regime;
 //       - eConvexMaster: delegate to f_master_solver->compute()
 //       - eMILPMaster:   run solve_MILP_master()
 throw( std::logic_error(
  "BendersDecompositionSolver::compute: not implemented yet" ) );
 }

/*--------------------------------------------------------------------------*/
/*------------------------ READING THE SOLUTION ----------------------------*/
/*--------------------------------------------------------------------------*/

Solver::OFValue BendersDecompositionSolver::get_lb( void )
{
 throw( std::logic_error(
  "BendersDecompositionSolver::get_lb: not implemented yet" ) );
 }

/*--------------------------------------------------------------------------*/

Solver::OFValue BendersDecompositionSolver::get_ub( void )
{
 throw( std::logic_error(
  "BendersDecompositionSolver::get_ub: not implemented yet" ) );
 }

/*--------------------------------------------------------------------------*/

bool BendersDecompositionSolver::has_var_solution( void )
{
 throw( std::logic_error(
  "BendersDecompositionSolver::has_var_solution: not implemented yet" ) );
 }

/*--------------------------------------------------------------------------*/

void BendersDecompositionSolver::get_var_solution( Configuration * solc )
{
 // TODO: read x from f_master, recover y^k from each subproblem, then
 //       map_back_solution()
 throw( std::logic_error(
  "BendersDecompositionSolver::get_var_solution: not implemented yet" ) );
 }

/*--------------------------------------------------------------------------*/

bool BendersDecompositionSolver::has_dual_solution( void )
{
 throw( std::logic_error(
  "BendersDecompositionSolver::has_dual_solution: not implemented yet" ) );
 }

/*--------------------------------------------------------------------------*/

void BendersDecompositionSolver::get_dual_solution( Configuration * solc )
{
 throw( std::logic_error(
  "BendersDecompositionSolver::get_dual_solution: not implemented yet" ) );
 }

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
 // TODO: scan f_Block, detect the complicating Variable x (v_x), build one
 //       BendersBFunction per sub-Block (build_BendersBFunction), then
 //       assemble the master (build_convex_master / build_MILP_master)
 throw( std::logic_error(
  "BendersDecompositionSolver::reformulate: not implemented yet" ) );
 }

/*--------------------------------------------------------------------------*/

void BendersDecompositionSolver::build_BendersBFunction( Index k )
{
 throw( std::logic_error(
  "BendersDecompositionSolver::build_BendersBFunction: not implemented yet" )
  );
 }

/*--------------------------------------------------------------------------*/

void BendersDecompositionSolver::build_convex_master( void )
{
 throw( std::logic_error(
  "BendersDecompositionSolver::build_convex_master: not implemented yet" ) );
 }

/*--------------------------------------------------------------------------*/

void BendersDecompositionSolver::build_MILP_master( void )
{
 throw( std::logic_error(
  "BendersDecompositionSolver::build_MILP_master: not implemented yet" ) );
 }

/*--------------------------------------------------------------------------*/

int BendersDecompositionSolver::solve_MILP_master( void )
{
 throw( std::logic_error(
  "BendersDecompositionSolver::solve_MILP_master: not implemented yet" ) );
 }

/*--------------------------------------------------------------------------*/

void BendersDecompositionSolver::map_back_solution( void )
{
 throw( std::logic_error(
  "BendersDecompositionSolver::map_back_solution: not implemented yet" ) );
 }

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
