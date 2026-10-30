# Session 112: loss-aware counting interfaces

30 October 2026. **ACTIVE research session; first verified checkpoint.**

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

## Continuing lanes

- Counting/search: preserve unique auxiliary completion, positivity and the
  cost of actual residue evaluation.
- Symmetry/enumeration: integrate shared-residue bank actions, universal
  ramified vertex kernels and the loss of marked bases.
- Norms/signs: investigate the signed eight-square/octonion mechanism against
  existing Paley/Fano carriers, without discarding orientation signs.
- Arithmetic clocks: pursue loss-aware modular and valuation interfaces while
  retaining the integer clock and source data.

These are active research directions, not announced theorems. The remaining
unrestricted tournament classification and Collatz convergence questions stay
open.
