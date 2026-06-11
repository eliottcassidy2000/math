# HYP-2411 - Eta, Gleason, and OCF are three cancellation gates

**Status:** OPEN synthesis / proof program.
**Source:** codex-2026-06-11-P1.
**Companions:** HYP-2409, HYP-2410, THM-466, THM-478, T780.

## Statement

Three recent threads share one proof shape:

1. **Euler eta gate:** the infinite sparse denominator
   `prod(1-q^n)` has pentagonal support and a reciprocal with subexponential
   partition growth. Random signs on the same support appear to produce positive
   reciprocal Lyapunov growth.
2. **Gleason Type II gate:** the invariant ring in `A,B` forces the low-weight
   prefix of a length-72 Type II weight enumerator to vanish, leaving a unique
   nonnegative scalar partition function.
3. **OCF / digit gate:** the tournament Hamiltonian-path count is an odd-cycle
   partition function whose first 2-adic layers sit in low-degree / Reed-Muller
   strata (THM-466, THM-478) before higher digits become pseudorandom.

Conjectural common principle: a hard counting object becomes proof-facing only
after one identifies the cancellation gate that keeps the reciprocal/partition
function from generic exponential behavior. The gate may be an infinite product
(eta), a finite invariant ring (Gleason), or an odd-cycle independence polynomial
(OCF).

## What the Quotient Preserves and Destroys

The quotient preserves:

- which low layers vanish;
- whether growth is subexponential or exponential;
- which algebraic symmetry is responsible for cancellation.

The quotient destroys:

- root locations for eta-like denominators;
- codeword support and neighborhood data for `[72,36,16]`;
- individual odd-cycle geometry for tournament OCF.

That destruction is exactly where proof work remains.

## Proof Route

Formalize a "cancellation gate" as a filtered algebra plus a partition function
whose generic signed perturbation has positive Lyapunov exponent, while the
structured object has zero or low-depth growth because the filtered algebra
forces initial layers to cancel.

Candidate test cases:

- eta vs random pentagonal signs;
- Gleason degree-72 gate vs random Type II-looking signed support;
- `H mod 4` Reed-Muller flat recurrences vs random Boolean functions on the
  tiling cube.

## Challenged Assumption

Do not assume the proof vertices are the original objects. The useful vertices
may be sign laws, weight layers, Gleason monomials, odd-cycle collections, flats,
automorphism-prime cases, or proof obligations. In this session the sign-law
quotient was most productive, but it is explicitly lossy.
