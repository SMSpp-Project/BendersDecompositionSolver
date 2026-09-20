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
  separation problem

### Changed

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
