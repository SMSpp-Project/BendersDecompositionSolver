/*--------------------------------------------------------------------------*/
/*------------------- File BendersDecompositionSolver.h --------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the BendersDecompositionSolver class, which implements the
 * Solver interface, in particular in its CDASolver version, automating the
 * Benders decomposition of a block-structured problem: the "complicating"
 * (first-stage) Variable are kept in a master problem, while each sub-Block
 * is projected out into its value function, represented by a BendersBFunction
 * that produces the corresponding Benders (optimality and feasibility) cuts.
 *
 * The user is assumed to be familiar with the algorithm: refer to
 *
 *  W. van Ackooij, A. Frangioni, W. de Oliveira "Inexact Stabilized Benders'
 *  Decomposition Approaches, with Application to Chance-Constrained Problems
 *  with Finite Support" Computational Optimization and Applications 65(3),
 *  637 - 669, 2016
 *
 * available at
 *
 * \link
 *  http://www.di.unipi.it/~frangio/abstracts.html#COAP15
 * \endlink
 *
 * for the stabilized (and possibly inexact) variant of the method, in which
 * the master problem is solved by a bundle-type approach, and to
 *
 *  D. Baena, J. Castro, A. Frangioni "Stabilized Benders Methods for
 *  Large-scale Combinatorial Optimization, with Application to Data Privacy"
 *  Management Science 66(7), 3051 - 3068, 2020
 *
 * available at
 *
 * \link
 *  http://www.di.unipi.it/~frangio/abstracts.html#ManSci18
 * \endlink
 *
 * for the case in which the master problem is a (stabilized) combinatorial
 * problem, solved as a MILP with Benders cuts added as lazy constraints.
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
/*----------------------------- DEFINITIONS --------------------------------*/
/*--------------------------------------------------------------------------*/

#ifndef __BendersDecompositionSolver
 #define __BendersDecompositionSolver
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "AbstractBlock.h"

#include "BendersBFunction.h"

#include "CDASolver.h"

#include "FRowConstraint.h"

#include "LinearFunction.h"

#include "UpdateSolver.h"

#include <list>

#include <map>

/*--------------------------------------------------------------------------*/
/*--------------------------- NAMESPACE ------------------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Structured Modeling System++ (SMS++)
namespace SMSpp_di_unipi_it
{
/*--------------------------------------------------------------------------*/
/*-------------------------- CLASS BendersDecompositionSolver --------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- GENERAL NOTES --------------------------------*/
/*--------------------------------------------------------------------------*/
/// CDASolver automating the Benders decomposition of a Block
/** BendersDecompositionSolver is a CDASolver that, attached to a suitably
 * structured Block (B), automatically constructs and solves its Benders
 * decomposition: the "complicating" (first-stage) Variable of (B) are kept in
 * a master problem, while the Variable of each sub-Block are projected out,
 * the resulting value function being represented by a BendersBFunction that
 * produces the corresponding Benders cuts.
 *
 * The Block (B) is assumed to have the following structure:
 *
 *  (O)   min / max { d( x ) + \sum_{k \in K} c^k( y^k ) :
 *                    x \in X , ( x , y^k ) \in Z^k , k \in K }
 *
 * where:
 *
 * - x are the "complicating" (first-stage) Variable; they are the
 *   ColVariable of the root (B) itself, which plays the role of the *master*;
 *
 * - d( x ), together with the (first-stage) Constraint X over x, is the
 *   Objective and the Constraint of the root (B);
 *
 * - each sub-Block (B^k), k \in K, is a *subproblem*: it has its own Variable
 *   y^k, its own Objective c^k( y^k ) and its own (second-stage) Constraint;
 *   it is coupled to the master only through the appearance of (some of) the
 *   x in (some of) its Constraint, i.e., Z^k is described by Constraint of
 *   the form l^k <= F^k x + E^k( y^k ) <= u^k .
 *
 * By projecting out y^k, (O) is reformulated as
 *
 *  (O')  min / max { d( x ) + \sum_{k \in K} v^k( x ) : x \in X }
 *
 * where the value function of the k-th subproblem
 *
 *  v^k( x ) = min / max { c^k( y^k ) :
 *                         ( l^k - F^k x ) <= E^k( y^k ) <= ( u^k - F^k x ) }
 *
 * is exactly what a BendersBFunction represents: the active Variable of the
 * BendersBFunction are the x, and the affine mapping x -> F^k x is injected
 * into the right-hand side of the coupling Constraint of (B^k). Each v^k() is
 * a convex (for minimization) C05Function whose linearizations are the
 * Benders optimality cuts (when the subproblem is bounded) and feasibility
 * cuts (Farkas certificates, when the subproblem is infeasible); these are
 * produced by BendersBFunction through the standard C05Function interface.
 *
 * BendersDecompositionSolver therefore does NOT re-implement the cut
 * generation, which lives in BendersBFunction. Its job is the *automation*
 * layer:
 *
 * 1. scan the root (B) to detect the complicating Variable x and, in each
 *    sub-Block (B^k), the Constraint that are linear in x; strip the F^k x
 *    term from those Constraint (recording F^k) and build one BendersBFunction
 *    v^k( x ) per sub-Block;
 *
 * 2. assemble an internal master Block representing (O') and solve it;
 *
 * 3. map the optimal x back to (B), together with the y^k recovered from the
 *    last solution of each subproblem.
 *
 * Two regimes are supported for how the master (O') is solved, selected by
 * the int_BDSlv_Regime parameter:
 *
 * - convex regime (the x are continuous, or a continuous relaxation is being
 *   solved): the master Objective is d( x ) + \sum_k v^k( x ), with the v^k()
 *   entering as C05Function; the master is handed to a CDASolver of the bundle
 *   family (e.g., BundleSolver), which drives the cutting-plane / stabilized
 *   cutting-plane loop pulling linearizations from the BendersBFunction. This
 *   is the setting of the COAP15 reference;
 *
 * - MILP regime (the x are integer): the master is a MILP with one epigraph
 *   Variable per subproblem (multi-cut) or a single aggregated one
 *   (single-cut), and the Benders cuts are added as dynamic Constraint; the
 *   master is solved by a MILPSolver, with BendersDecompositionSolver driving
 *   the outer cut loop (solve master, evaluate each v^k at the incumbent x,
 *   add the violated cuts, re-solve) until no violated cut remains. This is
 *   the setting of the ManSci18 reference.
 *
 * In both regimes the BendersBFunction objects are the same; only the master
 * assembly and the loop driver differ. The "inner" solver of the master is
 * created via the Solver factory and configured via Configuration, so that
 * "how the master is solved" is not hard-wired (see str_BDSlv_MSName).
 *
 * Current limitations (each is meant to become a pluggable extension point):
 *
 * - the coupling of each subproblem to x must be *affine* (linear terms in x
 *   inside FRowConstraint), as required by the BendersBFunction mapping;
 *
 * - the subproblems must be *convex* (typically LP), so that the LP-dual
 *   based Benders cuts are valid; integer subproblems (combinatorial /
 *   logic-based Benders) are not supported yet.
 *
 * Note that both kinds of cut are duals of the subproblem, hence the Solver
 * that is given to it has to be asked for a *vertex* solution: an
 * interior-point one has optimal but non-basic duals, which give a valid yet
 * weaker optimality cut, and above all it proves infeasibility without
 * producing the unbounded dual direction, i.e., the Farkas certificate, that
 * the feasibility cut is. The same happens if the infeasibility is detected
 * by the presolve rather than by the simplex, whence the presolve of the
 * subproblem Solver has to be switched off if feasibility cuts are wanted at
 * all; if the certificate is missing, exception is thrown rather than
 * silently converging to a wrong optimum. Note that a subproblem that is
 * always feasible, e.g. because the constraints that the master can make
 * unsatisfiable carry a slack with a large cost, needs no feasibility cut in
 * the first place, and is therefore the robust choice whenever the model
 * allows it. */

class BendersDecompositionSolver : public CDASolver
{

/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 public:

/*--------------------------------------------------------------------------*/
/*---------------------------- PUBLIC TYPES --------------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Public types
 *  @{ */

 // "import" basic types from Block
 using Index = Block::Index;
 using c_Index = Block::c_Index;
 using Subset = Block::Subset;
 using Range = Block::Range;

 /// regime used to solve the master problem (O')
 enum master_regime_type {
  eConvexMaster = 0 ,  ///< convex master, solved by a bundle-type CDASolver
  eMILPMaster   = 1    ///< MILP master, solved with lazy Benders cuts
  };

 /// how the Benders cuts are aggregated in the master
 enum cut_aggregation_type {
  eMultiCut  = 0 ,  ///< one epigraph Variable / cut family per subproblem
  eSingleCut = 1    ///< one aggregated epigraph Variable for all subproblems
  };

/** @} ---------------------------------------------------------------------*/
/*------------- CONSTRUCTING AND DESTRUCTING BendersDecompositionSolver ----*/
/*--------------------------------------------------------------------------*/
/** @name Constructing and destructing BendersDecompositionSolver
 *  @{ */

 /// constructor: does nothing besides initialising fields to defaults
 BendersDecompositionSolver( void );

/*--------------------------------------------------------------------------*/

 /// destructor: tears down the internal master Block and BendersBFunction
 ~BendersDecompositionSolver( void ) override;

/*--------------------------------------------------------------------------*/

 /// attach the Solver to the Block (B), triggering the reformulation
 /** Besides the standard Solver::set_Block() bookkeeping, this scans (B) as
  * described in the general notes, detects the complicating Variable x and
  * the coupling Constraint of each sub-Block, and builds the internal master
  * Block together with one BendersBFunction per sub-Block. */
 void set_Block( Block * block ) override;

/** @} ---------------------------------------------------------------------*/
/*-------------------- METHODS FOR SOLVING THE Block -----------------------*/
/*--------------------------------------------------------------------------*/
/** @name Solving the Block
 *  @{ */

 /// solve the Benders decomposition of (B)
 /** Drives the master solver according to the selected regime, returning the
  * usual Solver status codes. In the MILP regime this runs the outer Benders
  * cut loop; in the convex regime it delegates to the bundle-type master
  * solver, which performs the cutting-plane loop internally. */
 int compute( bool changedvars = true ) override;

/** @} ---------------------------------------------------------------------*/
/*---------------------- METHODS FOR READING RESULTS -----------------------*/
/*--------------------------------------------------------------------------*/
/** @name Reading the solution
 *  @{ */

 OFValue get_lb( void ) override;

 OFValue get_ub( void ) override;

 bool has_var_solution( void ) override;

 /// write the optimal x (and the recovered y^k) back into (B)
 void get_var_solution( Configuration * solc = nullptr ) override;

 bool has_dual_solution( void ) override;

 void get_dual_solution( Configuration * solc = nullptr ) override;

/** @} ---------------------------------------------------------------------*/
/*-------------------- METHODS FOR HANDLING PARAMETERS ---------------------*/
/*--------------------------------------------------------------------------*/
/** @name Handling parameters
 *
 * All the parameters of the master Solver are exposed as native parameters of
 * BendersDecompositionSolver, the indices being shifted past the
 * BendersDecompositionSolver-specific ones; the *_par_ms() helpers translate
 * between the two index spaces.
 *  @{ */

 /// public enum of the int parameters specific to BendersDecompositionSolver
 enum int_par_type_BDSlv {
  int_BDSlv_iBCopy = intLastParCDAS ,
  ///< copy (R3_Block) vs evict the sub-Block into the BendersBFunction

  int_BDSlv_Regime ,
  ///< master regime, a master_regime_type value

  int_BDSlv_CutType ,
  ///< cut aggregation, a cut_aggregation_type value

  int_BDSlv_FeasCut ,
  ///< 0 = always-feasible subproblem (slack), 1 = Farkas feasibility cuts

  int_BDSlv_MaxRounds ,
  ///< cap on the number of cut rounds in the MILP regime

  intLastBDSlvPar  ///< first allowed parameter value for derived classes
  };

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 /// public enum of the string parameters specific to BendersDecompositionSolver
 enum str_par_type_BDSlv {
  str_BDSlv_MSName = strLastParCDAS ,
  ///< classname of the master Solver ("BundleSolver", a "*MILPSolver", ...)

  str_Bsub_BSCfg ,
  ///< filename of the default BlockSolverConfig for the subproblems

  str_Mstr_BSCfg ,
  ///< filename of the BlockSolverConfig for the master Block

  strLastBDSlvPar  ///< first allowed parameter value for derived classes
  };

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 [[nodiscard]] idx_type get_num_int_par( void ) const override;

 [[nodiscard]] idx_type get_num_str_par( void ) const override;

 [[nodiscard]] int get_dflt_int_par( idx_type par ) const override;

 [[nodiscard]] const std::string & get_dflt_str_par( idx_type par )
  const override;

 [[nodiscard]] idx_type int_par_str2idx( const std::string & name )
  const override;

 [[nodiscard]] idx_type str_par_str2idx( const std::string & name )
  const override;

 [[nodiscard]] const std::string & int_par_idx2str( idx_type idx )
  const override;

 [[nodiscard]] const std::string & str_par_idx2str( idx_type idx )
  const override;

 void set_par( idx_type par , int value ) override;

 void set_par( idx_type par , const std::string & value ) override;

 [[nodiscard]] int get_int_par( idx_type par ) const override;

 [[nodiscard]] const std::string & get_str_par( idx_type par ) const override;

/** @} ---------------------------------------------------------------------*/
/*-------------------- PROTECTED PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 protected:

/*--------------------------------------------------------------------------*/
/*----------------------- PROTECTED METHODS --------------------------------*/
/*--------------------------------------------------------------------------*/

 /// scan (B), build the master Block and the per-subproblem BendersBFunction
 void reformulate( void );

 /// detect, in sub-Block sbi, the Constraint that are linear in x and turn
 /// them into the affine mapping (A_k, b_k) of the k-th BendersBFunction
 void build_BendersBFunction( Index k );

 /// assemble the convex master Objective d( x ) + \sum_k v^k( x )
 void build_convex_master( void );

 /// assemble the MILP master with epigraph Variable and the cut family
 void build_MILP_master( void );

 /// run the outer Benders cut loop driving the MILP master
 int solve_MILP_master( void );

 /// translate the master x and the subproblem y^k back into (B)
 void map_back_solution( void );

 /// create the master Solver out of str_Mstr_BSCfg and register it to (B)
 /** Creates the Solver named by the BlockSolverConfig in str_Mstr_BSCfg,
  * gives it the corresponding ComputeConfig and registers it to the master
  * Block. The Solver is registered *additively*, i.e., by hand rather than by
  * applying the BlockSolverConfig: applying it would replace the Solver
  * registered to (B), which is this BendersDecompositionSolver, and destroy
  * it in the middle of its own compute(). */

 void acquire_master_solver( void );

 /// apply to \p block the BlockSolverConfig, or meta-configuration, in \p fn
 /** Applies to \p block the Configuration in the file \p fn, which is either
  * a BlockSolverConfig, applied as it is, or a "meta-configuration", i.e., a
  * SimpleConfiguration< std::map< std::string , Configuration * > > mapping
  * the classname() of a Block to the BlockSolverConfig for it, which is
  * dispatched by classname() over the whole sub-tree. The latter is what
  * makes heterogeneous subproblems configurable. */

 void apply_BSCfg( Block * block , const std::string & fn );

 /// translate a master-Solver int parameter index into the BDSlv index space
 [[nodiscard]] idx_type int_par_ms( idx_type par ) const;

 /// translate a master-Solver str parameter index into the BDSlv index space
 [[nodiscard]] idx_type str_par_ms( idx_type par ) const;

/*--------------------------------------------------------------------------*/
/*----------------------- PROTECTED FIELDS ---------------------------------*/
/*--------------------------------------------------------------------------*/

 /// the internal master Block representing (O') (an AbstractBlock)
 AbstractBlock * f_master{};

 /// the master Solver, created via factory and attached to f_master
 CDASolver * f_master_solver{};

 /// one BendersBFunction v^k( x ) per sub-Block of (B)
 std::vector< BendersBFunction * > v_BF;

 /// the complicating (first-stage) Variable x, in master order
 std::vector< ColVariable * > v_x;

 /// the position of each complicating Variable in v_x
 std::map< const ColVariable * , Index > x_index;

 /// the epigraph Variable of the MILP master, one per cut family
 std::vector< ColVariable > * v_eta{};

 /// the Benders cuts of the MILP master
 std::list< FRowConstraint > * v_cuts{};

 /// true once the reformulation has been done
 bool f_reformulated = false;

 /// the value of the master at the last compute()
 OFValue f_value = 0;

 /// true if a solution of the master is available
 bool f_solved = false;

 // ----- parameters -------------------------------------------------------

 int f_iBCopy = 0;       ///< int_BDSlv_iBCopy

 int f_regime = eConvexMaster;  ///< int_BDSlv_Regime

 int f_cut_type = eMultiCut;    ///< int_BDSlv_CutType

 int f_feas_cut = 1;     ///< int_BDSlv_FeasCut

 int f_max_rounds = Inf< int >();  ///< int_BDSlv_MaxRounds

 std::string f_MSName;   ///< str_BDSlv_MSName

 std::string f_Bsub_BSCfg;  ///< str_Bsub_BSCfg

 std::string f_Mstr_BSCfg;  ///< str_Mstr_BSCfg

/*--------------------------------------------------------------------------*/
/*------------------------- PRIVATE PART OF THE CLASS ----------------------*/
/*--------------------------------------------------------------------------*/

 private:

/*--------------------------------------------------------------------------*/
/*---------------------------- PRIVATE FIELDS ------------------------------*/
/*--------------------------------------------------------------------------*/

 SMSpp_insert_in_factory_h;

/*--------------------------------------------------------------------------*/

 };  // end( class BendersDecompositionSolver )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

}  // end( namespace SMSpp_di_unipi_it )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#endif  /* BendersDecompositionSolver.h included */

/*--------------------------------------------------------------------------*/
/*------------- End File BendersDecompositionSolver.h ----------------------*/
/*--------------------------------------------------------------------------*/
