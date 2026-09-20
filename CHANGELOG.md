# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- `int_BDSlv_PhaseOneWeights`, the cost of the slacks of the phase one, which
  is the normalization that selects its cut: every slack costing one, as
  before and by default, or the slacks of each coupling row costing the
  inverse of the 2-norm of its coefficients, which makes the cut independent
  of how the rows are scaled

- `int_BDSlv_Unified`, which separates one cut for both feasibility and
  optimality, as Cornuejols and Lemarechal, Fischetti, Salvagnin and Zanette
  and Brandenberg and Stursberg do: the epigraph inequality joins the coupling
  Constraint that carry a slack in the phase one, the epigraph Variable enters
  its value function the way the complicating ones do, and the linearization
  in ( x , eta ), asked to be nonpositive, is the cut. Together with it,
  `dbl_BDSlv_EpiWeight`, the cost of the slack of the epigraph inequality,
  whose ratio to the costs of the coupling ones is the normalization of the
  separation problem. Since that normalization bounds the multipliers, a
  separation problem whose costs are far from the scale of the model can miss
  a cut, hence a round that separates none is not believed on its word: the
  value functions are evaluated, which is what leaves the solution of each
  subproblem where it is read from anyway, and the loop ends only if none of
  them is above its epigraph Variable. The ratio cannot be pushed arbitrarily
  far, though: the costs are coefficients of the separation problem, an
  extreme one makes it ill-conditioned, and multipliers that are only
  approximately feasible for its dual give a cut that is only approximately
  valid, which shows up as a master value above the optimum

- `int_BDSlv_Restore`, which gives the Block back as it was at the end of each
  compute(): the complicating Variable return to the Constraint they were
  taken out of, with the coefficient and the sides they had, and the master
  the Block was grafted into is disposed of. This is what lets another Solver
  be attached to the same Block, i.e., what a cross-check in a
  BlockSolverConfig is made of; what it costs is the reformulation, paid at
  every compute() rather than once, and the cuts, which are pieces of the
  master and go with it

### Changed

- the tests are no longer here, all three of them having to be posed on a
  problem, i.e., on a Block this module does not carry: the run of every way
  this Solver has of writing a cut, against the monolithic optimum of a small
  facility location instance, is in `tests/CapacitatedFacilityLocation`
  beside the comparison with the ad hoc Benders of that Block; the comparison
  of the two dual decompositions of a support vector machine, which is a test
  of what `SVMBlock::set_structure()` deals out, is in `tests/SVMBlock`; and
  the unit test of the pruning of the rows of a `PolyhedralFunction`, which
  names nothing of Benders, is in `tests/PolyhedralFunctionBlock`. A test of
  this Solver against itself, needing no problem, would belong here

- the two comparisons that are not about this Solver alone, i.e., the three
  forms of a two-stage stochastic investment problem and the two Benders
  decompositions of a facility location instance, move to the test suite of
  the Block each of them is posed on, where the batches of a Block live;
  what stays here is what tests the Solver itself

- taking the complicating Variable out of the Constraint of a subproblem, and
  putting them back, issues no Modification: the Solver of the subproblem is
  attached after the reformulation and reads it as it is then, while a Solver
  of the Block that was there before would be told of a change that is undone
  before anything is asked of it, and would rebuild its model out of the two
  halves of it. What the Modification did besides telling, i.e., the
  registration of the Constraint among the ones the Variable is active in, is
  now done explicitly

### Fixed

- the phase one minimized the violation of the coupling Constraint *plus the
  cost of the sub-Block of the subproblem*: the Objective of a Block is the
  sum of its own and of those of the Block it is made of, while only the one
  of the root of the replica was replaced by the total violation. On a
  subproblem that is itself a tree, which is the common case, the phase one
  was therefore measuring something else, and its cut could be invalid. The
  replica is now walked, every Objective in it is emptied, and their sum is
  what the unified cut writes its epigraph inequality with

- a sub-Block with no Objective, or whose Objective has no sense, is rejected
  instead of yielding cuts with the wrong sign: the sign of the linearization
  of a value function is read off the sense of the Objective of the
  subproblem

- the phase one never gave its cut: its slacks were added to the coupling
  rows without issuing the Modification, and a Constraint registers itself
  with a Variable coming in only when it receives it, so that a Solver
  building its model by columns, as the :MILPSolver do, left the slacks out
  of the rows. The phase one was then infeasible at every x, and every
  infeasible subproblem was silently cut away by the Farkas certificate
  instead. The test now checks the phase one on an instance whose infeasible
  subproblems have more than one extreme ray, where its cut differs from the
  certificate

## [0.1.1] - 2026-09-13

### Fixed

- the master needs the least numerical care of GUROBI, not none of it, and
  the residual is declared

## [0.1.0] - 2026-09-12

### Added

- initial scaffold of the BendersDecompositionSolver module: CDASolver
  interface, parameter handling and the skeleton of the reformulation and of
  the master cut loop

### Changed

- the replica of the subproblem that the phase one of the feasibility cut
  works on is the copy of its abstract representation that the Block makes of
  itself [see `AbstractBlock::mirror()`], rather than one assembled here out
  of the pieces this Solver knows how to read: the copy keeps the shape of
  the groups, the fixings and the Constraint that are not linear rows, and
  says what it could not reproduce instead of dropping it in silence

- the version of the module is the git tag of its repository, or the
  VERSION.txt of a release tarball, and the shared library carries it: its
  SONAME is major.minor while the major is 0, and it is installed with an
  RPATH relative to itself, so that an installed tree keeps working wherever
  it is moved

[Unreleased]: https://gitlab.com/smspp/bendersdecompositionsolver/-/compare/0.1.1...develop
[0.1.1]: https://gitlab.com/smspp/bendersdecompositionsolver/-/compare/0.1.0...0.1.1
[0.1.0]: https://gitlab.com/smspp/bendersdecompositionsolver/-/tags/0.1.0
