# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

### Changed

### Fixed

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

[Unreleased]: https://gitlab.com/smspp/bendersdecompositionsolver/-/compare/0.1.0...develop
[0.1.0]: https://gitlab.com/smspp/bendersdecompositionsolver/-/tags/0.1.0
