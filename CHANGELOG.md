# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

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

### Fixed
