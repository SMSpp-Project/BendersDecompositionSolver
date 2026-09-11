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
 * problem, solved as a MILP with Benders cuts added as lazy constraints. The
 * Pareto-optimal cuts are those of
 *
 *  T.L. Magnanti, R.T. Wong "Accelerating Benders Decomposition: Algorithmic
 *  Enhancement and Model Selection Criteria" Operations Research 29(3),
 *  464 - 484, 1981
 *
 * in the one-step form of
 *
 *  N. Papadakos "Practical Enhancements to the Magnanti-Wong Method"
 *  Operations Research Letters 36(4), 444 - 449, 2008
 *
 * the combinatorial cuts are those of
 *
 *  G. Codato, M. Fischetti "Combinatorial Benders' Cuts for Mixed-Integer
 *  Linear Programming" Operations Research 54(4), 756 - 766, 2006
 *
 * and the reason why a feasibility cut is normalized, i.e., that which of the
 * many cuts a given infeasibility offers is selected is a matter of the
 * normalization imposed on the certificate, is discussed in
 *
 *  M. Fischetti, D. Salvagnin, A. Zanette "A Note on the Selection of
 *  Benders' Cuts" Mathematical Programming 124(1-2), 175 - 182, 2010
 *
 * a broad survey of all this being
 *
 *  R. Rahmaniani, T.G. Crainic, M. Gendreau, W. Rei "The Benders
 *  Decomposition Algorithm: A Literature Review" European Journal of
 *  Operational Research 259(3), 801 - 817, 2017
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

#include <unordered_set>

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
 * 1. take the complicating Variable x and the subproblems as the two
 *    vector-of-int parameters say [see vintMasterBlock and vintMasterVars],
 *    the same Block admitting many Benders reformulations and which one is
 *    wanted being a choice of whoever poses the problem rather than something
 *    that can be read off it; find, in each subproblem (B^k), the Constraint
 *    that are linear in x, strip the F^k x term from them (recording F^k) and
 *    build one BendersBFunction v^k( x ) per subproblem;
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
 * Where the master is assembled differs, too, because the two regimes need
 * to add different things to (B), and adding something to a Block, which is
 * the jealous guardian of its own contents, means embedding it into a larger
 * one:
 *
 * - in the MILP regime what is added are the epigraph Variable and the cuts,
 *   which do not belong to (B) at all: the master is therefore a Block of
 *   this Solver's own, holding them, into which (B) is *grafted* as its only
 *   sub-Block. The master Solver sees the x, the X and the d( x ) through
 *   (B), and (B) is left exactly as it was, whence *any* Block can be the
 *   root here;
 *
 * - in the convex regime what is added are the sub-Block carrying the value
 *   functions, and they have to be sub-Block of the Block that holds the x,
 *   for that is the structure a bundle-type Solver reads: the master is (B)
 *   itself, which therefore has to be an AbstractBlock, the only Block that
 *   lets a Solver add a sub-Block to it.
 *
 * In neither regime are the subproblems removed from (B): they now live
 * inside the BendersBFunction, and are declared *excluded* to the master
 * Solver [see Solver::set_excluded_blocks()], which is what "the master does
 * not see them" means without (B) having to be mutilated to say it.
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
 * Both kinds of cut are duals of the subproblem, but they ask different things
 * of the Solver that produces them. The optimality cut is a linearization of
 * the value function, and *any* optimal dual solution gives a valid one,
 * basic or not: an interior-point method, whose solution is "central" rather
 * than a vertex, is not only allowed but can be expected to give stronger
 * cuts. The feasibility cut, instead, *is* the unbounded dual direction, i.e.,
 * the Farkas certificate, which is a ray of the dual polyhedron and therefore
 * has nothing to do with the algorithm: a Solver worth its name produces one
 * from an interior point as well. What does destroy it is the presolve, which
 * detects the infeasibility on the reduced problem and returns no certificate
 * for the original one, whence the presolve of the subproblem Solver has to
 * be switched off if feasibility cuts are wanted at all. If the certificate is
 * missing, exception is thrown rather than silently converging to a wrong
 * optimum. Note that a subproblem that is always feasible, e.g. because the
 * constraints that the master can make unsatisfiable carry a slack with a
 * large cost, needs no feasibility cut in the first place, and is therefore
 * the robust choice whenever the model allows it, nothing being then asked of
 * the subproblem Solver beyond solving it. */

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

 /// what an infeasible subproblem is cut away with
 /** A subproblem that has no solution at the current x says that that x is
  * not feasible for the original problem either, and the master has to be
  * told so with a cut that no optimality one can give, the value function
  * being infinite there. There are three ways out:
  *
  * - the subproblem is never infeasible, which is the case whenever the
  *   constraints the master can make unsatisfiable carry a slack with a large
  *   enough cost: nothing is asked of the subproblem Solver then, and this is
  *   the robust choice whenever the model allows it;
  *
  * - the Farkas certificate of the infeasibility, i.e., the unbounded ray of
  *   the dual, which is the classical feasibility cut and the tightest of the
  *   three, but has to be produced by the subproblem Solver;
  *
  * - a combinatorial, or no-good, cut, which merely forbids the current
  *   assignment: it asks nothing of the subproblem Solver, but it is only
  *   available when the complicating Variable are all binary, and it is much
  *   weaker, cutting away one point at a time [Codato and Fischetti];
  *
  * - the cut of a *phase one*, i.e., of the problem that minimizes the total
  *   violation of the coupling Constraint by giving each of them a slack of
  *   unit cost. That problem is feasible whatever x, its value is zero
  *   exactly where the subproblem has a solution, and it is a value function
  *   like any other: its linearization at the incumbent, asked to be
  *   nonpositive, is a feasibility cut. Which cut a given infeasibility
  *   yields is decided by the normalization the certificate is subject to,
  *   and here the multipliers are bounded by the unit costs, which is the
  *   normalization Fischetti, Salvagnin and Zanette recommend: the cut is
  *   therefore not an arbitrary ray but the one a bounded separation problem
  *   selects. The subproblem is *replicated*, so nothing of the user's is
  *   touched, and the replica is the copy of its abstract representation
  *   that the Block makes of itself [see AbstractBlock::mirror()] rather
  *   than something assembled here: whatever that copy cannot reproduce is
  *   left out and said, which makes the phase one a relaxation, and a
  *   relaxation still cuts no feasible x; when its cut does not cut the
  *   incumbent the Farkas certificate is asked for instead. */

 enum feasibility_cut_type {
  eAlwaysFeasible = 0 ,  ///< the subproblem cannot be infeasible
  eFarkas         = 1 ,  ///< the Farkas certificate of the infeasibility
  eCombinatorial  = 2 ,  ///< a no-good cut on the binary complicating Variable
  ePhaseOne       = 3    ///< the cut of the minimum-violation problem
  };

 /// which of the many cuts a degenerate subproblem offers is taken
 /** A subproblem is very often dual degenerate, and then its optimal duals
  * are a face rather than a point: every one of them gives a valid cut, but
  * some of those cuts dominate the others, and which one the Solver happens
  * to return is arbitrary. A cut that is *not* dominated by any other is
  * called Pareto-optimal.
  *
  * [
  *   v^k( x ) \geq lpha + g^	op x \quad \mbox{dominates} \quad
  *   v^k( x ) \geq lpha' + g'^	op x
  * ]
  * when it is above it everywhere on the master feasible set, so which cut
  * is Pareto-optimal depends on that set and not on the incumbent alone.
  * Magnanti and Wong obtain one by maximizing the cut at a *core point*, a
  * point in the relative interior of the master feasible set, subject to the
  * duals staying optimal for the incumbent; Papadakos observed that the extra
  * constraint can be dropped, the cut generated by solving the subproblem at
  * the core point alone being Pareto-optimal already, which turns the second
  * problem into another ordinary evaluation of the value function. */

 enum cut_strengthening_type {
  eNoPareto  = 0 ,  ///< the cut the incumbent gives, whichever it is
  ePapadakos = 1    ///< one more cut per round, generated at the core point
  };

 /// how a feasibility cut is scaled before it enters the master
 /** A feasibility cut is the Farkas certificate of an infeasible subproblem,
  * i.e., a ray of its dual polyhedron: a ray is defined up to a positive
  * multiplier, so how large the coefficients of the cut come out is an
  * accident of how the Solver normalized the certificate. The half-space the
  * cut describes does not change with the scaling, but everything the master
  * measures on it does, from the violation that decides whether it is
  * separated to the numerics of the row itself, whence dividing it by a norm
  * of its coefficients puts all the feasibility cuts on the same footing.
  * Which cut a given infeasibility yields is in fact decided by the
  * normalization imposed on the certificate, as Fischetti, Salvagnin and
  * Zanette show by writing the separation as a problem of its own.
  * Optimality cuts are not scaled, and cannot be: the epigraph Variable in
  * them has coefficient one, which fixes their scale. */

 enum feasibility_cut_norm_type {
  eNoNorm  = 0 ,  ///< the cut as the certificate comes out of the Solver
  eOneNorm = 1 ,  ///< divided by the 1-norm of its coefficients
  eTwoNorm = 2 ,  ///< divided by the 2-norm of its coefficients
  eInfNorm = 3    ///< divided by the largest of its coefficients
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

/*--------------------------------------------------------------------------*/
 /// the number of iterations of the cutting-plane loop
 /** The loop is run by the Solver of the master, be it the bundle of the
  * convex regime or this Solver itself in the MILP one: it is therefore the
  * count of the master Solver in the former case and the number of rounds of
  * the cutting-plane loop in the latter, and it is zero if nothing has been
  * solved yet. */

 long get_elapsed_iterations( void ) const override;

/*--------------------------------------------------------------------------*/
 /// the number of Benders cuts generated in the MILP regime
 /** The number of cuts the cutting-plane loop of the MILP regime has added to
  * the master so far, which is what tells a multi-cut run from a single-cut
  * one: they may well do the same number of iterations while adding a very
  * different number of cuts. In the convex regime the cuts are the
  * linearizations the master Solver asks for, and they are counted by it. */

 [[nodiscard]] long get_num_cuts( void ) const { return( f_cuts ); }

/*--------------------------------------------------------------------------*/
 /// the time spent in the cutting-plane loop, master and subproblems

 double get_elapsed_time( void ) const override;

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
  ///< how an infeasible subproblem is cut away, a feasibility_cut_type value

  int_BDSlv_MaxRounds ,
  ///< cap on the number of cut rounds in the MILP regime

  int_BDSlv_Pareto ,
  ///< cut strengthening, a cut_strengthening_type value

  int_BDSlv_CutNorm ,
  ///< scaling of the feasibility cuts, a feasibility_cut_norm_type value

  intLastBDSlvPar  ///< first allowed parameter value for derived classes
  };

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 /// public enum of the double parameters specific to BendersDecompositionSolver
 enum dbl_par_type_BDSlv {
  dbl_BDSlv_CoreMove = dblLastParCDAS ,
  ///< how far the core point moves towards the incumbent, in [ 0 , 1 ]

  dblLastBDSlvPar  ///< first allowed parameter value for derived classes
  };

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 /// public enum of the vector-of-int parameters specific to this Solver
 /** Which sub-Block are the master and which Variable are the complicating
  * ones is not something that can be read off a Block: the same Block admits
  * many Benders reformulations, and which one is wanted is a choice of
  * whoever poses the problem. These two say it, in the simplest way that
  * covers the models whose subproblems are sub-Block of the master: the
  * *positions* of the sub-Block that belong to the master, the subproblems
  * being the ones that are left, and the positions of the complicating
  * Variable among those the master exposes.
  *
  * Both empty, which is the default, is the convention this Solver had
  * before they existed: the root is the master, every sub-Block of it is a
  * subproblem, and every ColVariable of the root is complicating. */

 enum vint_par_type_BDSlv {
  vintMasterBlock = vintLastParCDAS ,
  ///< positions of the sub-Block that are part of the master

  vintMasterVars ,
  ///< positions, among those of the master, of the complicating Variable

  vintLastBDSlvPar  ///< first allowed parameter value for derived classes
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

 [[nodiscard]] idx_type get_num_dbl_par( void ) const override;

 [[nodiscard]] idx_type get_num_vint_par( void ) const override;

 [[nodiscard]] idx_type get_num_str_par( void ) const override;

 [[nodiscard]] int get_dflt_int_par( idx_type par ) const override;

 [[nodiscard]] double get_dflt_dbl_par( idx_type par ) const override;

 [[nodiscard]] const std::vector< int > & get_dflt_vint_par( idx_type par )
  const override;

 [[nodiscard]] const std::string & get_dflt_str_par( idx_type par )
  const override;

 [[nodiscard]] idx_type dbl_par_str2idx( const std::string & name )
  const override;

 [[nodiscard]] idx_type int_par_str2idx( const std::string & name )
  const override;

 [[nodiscard]] idx_type vint_par_str2idx( const std::string & name )
  const override;

 [[nodiscard]] idx_type str_par_str2idx( const std::string & name )
  const override;

 [[nodiscard]] const std::string & dbl_par_idx2str( idx_type idx )
  const override;

 [[nodiscard]] const std::string & int_par_idx2str( idx_type idx )
  const override;

 [[nodiscard]] const std::string & vint_par_idx2str( idx_type idx )
  const override;

 [[nodiscard]] const std::string & str_par_idx2str( idx_type idx )
  const override;

 void set_par( idx_type par , int value ) override;

 void set_par( idx_type par , double value ) override;

 void set_par( idx_type par , const std::string & value ) override;

 void set_par( idx_type par , std::vector< int > && value ) override;

 [[nodiscard]] int get_int_par( idx_type par ) const override;

 [[nodiscard]] double get_dbl_par( idx_type par ) const override;

 [[nodiscard]] const std::string & get_str_par( idx_type par ) const override;

 [[nodiscard]] const std::vector< int > & get_vint_par( idx_type par )
  const override;

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

/*--------------------------------------------------------------------------*/
 /// build the phase-one value function of subproblem \p k
 /** Replicates subproblem \p k into an AbstractBlock of this Solver's own,
  *  gives each of its coupling Constraint, the ones at the positions \p cpl
  *  in the order the Constraint of the subproblem are scanned, a slack of
  *  unit cost, and makes the sum of those slacks the Objective: the value of
  *  the problem is then the least total violation of the coupling, which is
  *  zero exactly where the subproblem has a solution. The affine mapping is
  *  the one of the subproblem, \p A and \p b, on the same \p sides.
  *
  *  Only ColVariable and linear FRowConstraint are replicated: anything else
  *  is left out, which turns the replica into a relaxation and the cut into a
  *  weaker, but still valid, one [see feasibility_cut_type]. */

 void build_phase_one( Index k , const Subset & cpl ,
                       const std::vector< std::vector< double > > & A ,
                       const std::vector< double > & b ,
                       const std::vector< int > & sides );

 /// assemble the convex master Objective d( x ) + \sum_k v^k( x )
 void build_convex_master( void );

 /// assemble the MILP master with epigraph Variable and the cut family
 void build_MILP_master( void );

 /// run the outer Benders cut loop driving the MILP master
 int solve_MILP_master( void );

 /// give (B) back what belongs to it and release what this Solver owns
 /** Undoes what reformulate() has assembled: the master Solver is released,
  * each subproblem is given back to (B), which is its owner, and the master
  * Block is dismantled, i.e., deleted if this Solver has built it and
  * stripped of the value-function sub-Block if it is (B) itself.
  *
  * What is *not* undone is the reformulation of the subproblems: the x terms
  * that have been stripped out of their Constraint stay inside the
  * BendersBFunction, i.e., (B) is consumed by this Solver as far as the
  * coupling Constraint are concerned. */

 void dismantle( void );

 /// translate the master x and the subproblem y^k back into (B)
 void map_back_solution( void );

 /// create the master Solver out of str_Mstr_BSCfg and register it
 /** Creates the Solver named by the BlockSolverConfig in str_Mstr_BSCfg,
  * gives it the corresponding ComputeConfig, tells it which sub-Block it has
  * to ignore and registers it to the master Block. The Solver is registered
  * *additively*, i.e., by hand rather than by applying the
  * BlockSolverConfig: applying it would replace the Solver registered to
  * (B), which is this BendersDecompositionSolver, and destroy it in the
  * middle of its own compute(). */

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

 /// the master Block representing (O'): (B) itself, or a Block of this own
 AbstractBlock * f_master{};

 /// the sub-Block of (B) the master Solver has to ignore
 /** The subproblems now live inside the BendersBFunction, but they are still
  * sub-Block of (B): rather than being removed from it, which only an
  * AbstractBlock would allow and which would mutilate (B), they are declared
  * excluded to the master Solver [see Solver::set_excluded_blocks()]. */

 std::unordered_set< Block * > f_ignored;

 /// the father (B) had before being grafted into the master, if it was
 Block * f_Block_father{};

 /// the master Solver, created via factory and attached to f_master
 CDASolver * f_master_solver{};

 /// one BendersBFunction v^k( x ) per sub-Block of (B)
 std::vector< BendersBFunction * > v_BF;

 /// the sub-Block carrying the value functions in the convex master
 std::vector< AbstractBlock * > v_wrap;

 /// the complicating (first-stage) Variable x, in master order
 std::vector< ColVariable * > v_x;

 /// the position of each complicating Variable in v_x
 std::map< const ColVariable * , Index > x_index;

 /// which sub-Block of (B) are the subproblems, in the order they are dealt
 /// with: everything that vintMasterBlock does not keep in the master

 std::vector< Index > v_sub;

 /// the epigraph Variable of the MILP master, one per cut family
 std::vector< ColVariable > * v_eta{};

 /// the Benders cuts of the MILP master
 std::list< FRowConstraint > * v_cuts{};

 /// the phase-one value function of each subproblem, if it is being used
 /** One per subproblem, each on a replica of it that carries a slack of unit
  * cost on every coupling Constraint [see feasibility_cut_type]; empty unless
  * int_BDSlv_FeasCut says #ePhaseOne. */

 std::vector< BendersBFunction * > v_BF1;

 /// the replicas the phase-one value functions are built on
 std::vector< Block * > v_phase1;

 /// true once the reformulation has been done
 bool f_reformulated = false;

 /// the value of the master at the last compute()
 OFValue f_value = 0;

 /// true if a solution of the master is available
 bool f_solved = false;

 /// the core point the Pareto-optimal cuts are generated at, in master order
 std::vector< double > v_core;

 /// the number of cuts added to the MILP master
 long f_cuts = 0;

 /// the number of rounds of the cutting-plane loop of the MILP regime
 long f_rounds = 0;

 // ----- parameters -------------------------------------------------------

 int f_iBCopy = 0;       ///< int_BDSlv_iBCopy

 int f_regime = eConvexMaster;  ///< int_BDSlv_Regime

 int f_cut_type = eMultiCut;    ///< int_BDSlv_CutType

 int f_feas_cut = eFarkas;  ///< int_BDSlv_FeasCut

 int f_max_rounds = Inf< int >();  ///< int_BDSlv_MaxRounds

 int f_pareto = eNoPareto;  ///< int_BDSlv_Pareto

 double f_core_move = 0.5;  ///< dbl_BDSlv_CoreMove

 int f_cut_norm = eNoNorm;  ///< int_BDSlv_CutNorm

 std::string f_MSName;   ///< str_BDSlv_MSName

 std::string f_Bsub_BSCfg;  ///< str_Bsub_BSCfg

 std::string f_Mstr_BSCfg;  ///< str_Mstr_BSCfg

 std::vector< int > v_master_block;  ///< vintMasterBlock

 std::vector< int > v_master_vars;   ///< vintMasterVars

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
