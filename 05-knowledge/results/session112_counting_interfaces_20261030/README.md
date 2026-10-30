# Session 112: loss-aware counting interfaces

30 October 2026. **ACTIVE research session; second verified checkpoint.**

This session uses Rudich–Rousseau's *Deterministic Factoring in Time N^{o(1)}*
as inspiration for continuing the repository's symmetry, enumeration, modular
clock and arithmetic work. [Source scope](paper/PAPER_SCOPE.md) separates the
paper's theorem, conjectured stronger rate, inspected formal source and
external registry-reported checks. We have not independently built its Lean
proof or reproduced the general counting engine. No new factoring complexity
result is claimed in this packet.

## First proved interface

[Modular observers](arithmetic/MODULAR_OBSERVERS.md) identify exactly what is
retained by observations of a scaled count sD: the count modulo
lcm(m_i)/gcd(lcm(m_i),s). The note derives an elementary cutoff avoiding the
prime support of the scale, an exact error/erasure correction criterion, and
prefix-sensitive observation budgets. Burnside's fixed-point numerator gives
a direct bridge to stabilizer-sensitive symmetry counts. Equal group orders
or equal prime supports do not identify the groups or their actions.

The [checker](arithmetic/verify_modular.py) supplies exhaustive bounded
kernel and coding tests, cutoff checks and a nonfree Burnside control.
[Its output](arithmetic/modular_checks.json) records those finite universes.
Run from this directory:

    python3 arithmetic/verify_modular.py

The elementary proofs, rather than finite checks alone, support the stated
all-integer interfaces. Cheap modular observations do not imply that computing
the observed counts is cheap.

## Signed constraints, selectors and canonical counts

[Norm signs](norms/NORM_SIGNS.md) give a necessary and sufficient binary
rectangle condition for Euclidean norm composition in a normalized signed-XOR
basis. Four named probes recover the promised-valid eight-square table; a
five-equation contradiction excludes dimension sixteen in this model. An
explicit coordinate map transports the old divisor-sum Hamming bank to these
16 valid tables, with squareclass 30 selecting the composing orientations.
This retains multiplication signs rather than only unsigned Fano support.

[Random selectors](symmetry/PROOFS.md), conditional on the inherited universal
ramified-group classification, separate residue coverage from the remaining
reflection phase. At 625 vertices exactly 1/125 of random minimal transversals
retain that phase; one pin move repairs every bad transversal. Exact sampling
probabilities and expected stopping times retain the depth-one exception.
These are colored-carrier base counts, not a classification of every single
tournament's automorphism group.

[Canonical ranking](counting/README.md) turns cumulative saturated cuts into
an exact rank/unrank oracle for tournament orbits. Ordinary prefix orbit
counts are not additive. Saturation and exact union indicators repair the
semantics, but may have enormous representation cost. All cuts through five
vertices and selected affine-rank implementations are finite controls; no
fast unrestricted enumeration algorithm follows.

The three lanes have independently replayed normal and optimized checkers:
[norm checks](norms/checks.json), [selector checks](symmetry/results.json), and
[counting checks](counting/results.json). Their proofs and inherited
classification assumptions remain explicit in the corresponding notes.

## Continuing lanes

- Counting/search: preserve unique auxiliary completion, positivity and the
  cost of actual residue evaluation.
- Symmetry/enumeration: integrate shared-residue bank actions, universal
  ramified vertex kernels and the loss of marked bases.
- Norms/signs: use the exact syndrome transport to test further coherent
  arithmetic/algebraic interfaces, retaining signs and validity promises.
- Arithmetic clocks: pursue loss-aware modular and valuation interfaces while
  retaining the integer clock and source data.

These are active research directions, not announced theorems. The remaining
unrestricted tournament classification and Collatz convergence questions stay
open.
