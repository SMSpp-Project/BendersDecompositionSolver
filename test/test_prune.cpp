/*--------------------------------------------------------------------------*/
/*--------------------------- File test_prune.cpp --------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Unit test for the cut/row pruning that used to live in the standalone
 * CutProcessing class and now lives, properly, on the classes that own the
 * data: PolyhedralFunction::remove_parallel_rows() (geometric) and the static
 * PolyhedralFunctionBlock::remove_redundant_rows() (LP-based, the epigraph).
 *
 * A convex PolyhedralFunction is built with a known parallel (dominated) row
 * and a known inactive row, and the pruning is checked to remove exactly them.
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

#include <iostream>

#include "BlockSolverConfig.h"
#include "ColVariable.h"
#include "PolyhedralFunction.h"
#include "PolyhedralFunctionBlock.h"

using namespace SMSpp_di_unipi_it;

int main( void )
{
 // f( x ) = max { rows } over R^2, convex:
 //  0: x0      1: -x0     2: x1      3: -x1   ( grow in all directions )
 //  4: x0 - 3            ( parallel to row 0, dominated -> remove_parallel )
 //  5: -9               ( constant below the others -> remove_redundant )
 std::vector< ColVariable > x( 2 );
 PolyhedralFunction::VarVector xv{ & x[ 0 ] , & x[ 1 ] };
 PolyhedralFunction::MultiVector A{ { 1 , 0 } , { -1 , 0 } , { 0 , 1 } ,
				    { 0 , -1 } , { 1 , 0 } , { 0 , 0 } };
 PolyhedralFunction::RealVector b{ 0 , 0 , 0 , 0 , -3 , -9 };

 PolyhedralFunction f( std::move( xv ) , std::move( A ) , std::move( b ) ,
		       - Inf< PolyhedralFunction::FunctionValue >() , true );

 const auto n0 = f.get_nrows();
 std::cout << "rows: " << n0;

 // geometric: drops the parallel dominated row ( row 4 )
 f.remove_parallel_rows();
 const auto n1 = f.get_nrows();
 std::cout << " -> after parallel: " << n1;

 // LP-based: drops the inactive row ( row 5 )
 auto c = Configuration::deserialize( "BSPar_sub.txt" );
 auto cfg = dynamic_cast< BlockSolverConfig * >( c );
 if( ! cfg ) { std::cerr << "BSPar_sub.txt not a BlockSolverConfig\n";
               return( 1 ); }
 PolyhedralFunctionBlock pfb( nullptr , & f );
 pfb.remove_redundant_rows( cfg );
 const auto n2 = f.get_nrows();
 std::cout << " -> after redundant: " << n2 << std::endl;
 cfg->clear();
 delete cfg;

 const bool ok = ( n0 == 6 ) && ( n1 == 5 ) && ( n2 == 4 );
 std::cout << ( ok ? "-> OK ( parallel + inactive rows removed )"
		   : "-> FAIL" ) << std::endl;
 return( ok ? 0 : 1 );
 }

/*--------------------------------------------------------------------------*/
/*------------------------ End File test_prune.cpp -------------------------*/
/*--------------------------------------------------------------------------*/
