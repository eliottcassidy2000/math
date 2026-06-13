# HYP-2452 - Factor-capture witnesses prune integral convolution lifts over coefficient tilings

**Status:** OPEN synthesis; exact degree-4/5 convolution scout verified.
**Source:** codex-2026-06-12.
**Companions:** HYP-2451, HYP-2450, HYP-2449, HYP-2448, HYP-2447,
HYP-2445, HYP-2444, HYP-2430.
**Computation:** `04-computation/convolution_factor_capture_tiling_codex.py`;
stored output `05-knowledge/results/convolution_factor_capture_tiling_codex.out`.
**External anchors:** Singh arXiv:2411.18366; Iravanian arXiv:2410.15880;
Abu Salem-Gao-Lauder ISSAC 2004 polytope factorization.

## Statement

HYP-2451 establishes the residue/valuation split-survivor carrier for hidden
convolution lifts.  This addendum keeps that carrier and adds the complementary
integer side: exact degree-4/5 integral lift search plus factor-capture
witnesses from values `f(m)`.

HYP-2450 projects fixed-path tournament tilings to coefficient profiles.  This
hypothesis asserts that the next proof-bearing object is the hidden integral
convolution lift behind the coefficient profile.

For a target polynomial

```text
f(x) = a_0 + a_1 x + ... + a_n x^n
```

an integral two-factor lift is a pair of nonconstant integer coefficient rows

```text
g(x) = b_0 + ... + b_r x^r,
h(x) = c_0 + ... + c_s x^s,
r+s=n,
```

such that every coefficient is a diagonal sum

```text
a_k = sum_{i+j=k} b_i c_j.
```

Thus reducibility is an integer tiling problem.  The target coefficient vector
is a boundary-total shadow; a factorization is a hidden 2D grid whose diagonals
sum to that shadow.  Irreducibility means no nontrivial integral lift exists.

This is the coefficient analogue of the repo's scalar/support split:

```text
coefficient profile       -> scalar boundary totals
factor coefficient grid   -> hidden support/incidence lift
irreducibility            -> no nontrivial integral lift
```

## Exact Finite Evidence

The pilot implements a hand-derived convolution-lift search that is exact for
primitive polynomials of degree at most 5.  In that range, any reducible
polynomial has a linear or quadratic factor, so boundary divisors plus the
inward diagonal equations are complete.

Stored scans:

```text
degree 4, |a_i| <= 3, primitive leading-positive rows:
  rows = 3856
  sympy reducible = 792
  convolution lift yes = 792
  mismatches = 0
  fixed-divisor blocked = 768

degree 5, |a_i| <= 2, primitive leading-positive rows:
  rows = 2016
  sympy reducible = 488
  convolution lift yes = 488
  mismatches = 0
  fixed-divisor blocked = 480
```

The degree-4 example

```text
x^4 - 2x^3 - 3x^2 - 3x - 3 = (x+1)(x^3 - 3x^2 - 3)
```

is displayed as a literal lift.  The degree-5 example

```text
x^5 + x^4 - 2x^3 - 2x^2 - 2x - 2
  = (x+1)(x^4 - 2x^2 - 2)
```

confirms the same boundary-to-interior recursion in odd degree.

## Factor-Capture Hypertournament

At a witness integer `m`, the value

```text
f(m) = product p_i^{e_i}
```

breaks into prime tokens.  Any factorization `f=gh` forces those tokens to be
allocated among the values `g(m)` and `h(m)`.  This is more naturally a
hypertournament/allocation game than a graph: vertices can be token blocks,
factor slots, or proof obligations, and an edge direction compares the number
of unresolved allocations after choosing a witness.

The script records the cheap score

```text
Omega(f(m)) = sum_i e_i
```

and the two-factor allocation count

```text
product_i (e_i+1) - 2.
```

Examples:

```text
x^4 - x^3 - x^2 - x - 1:
  m=3 gives f(m)=41, Omega=1, allocation_count=0, irreducible degree [4].

2x^4 - 2x^3 - x^2 - x - 1:
  m=3 gives f(m)=95=5*19, Omega=2, allocation_count=2,
  actual factor degrees [2,2].

1+x+...+x^8:
  m=1 gives f(m)=9, Omega=2, allocation_count=1,
  factor degrees [2,6].

1+x+...+x^10:
  m=1 gives f(m)=11, Omega=1, allocation_count=0,
  irreducible degree [10].
```

Singh's value-factor philosophy supplies the anchor: sufficiently informative
integer values and root-location data can bound the number of irreducible
factors.  The repo use is computational and conservative: low `Omega(f(m))`
is a pruning score for possible factor slots, not a proof by itself unless the
paper's hypotheses are checked.

## Residue And Sign-Cube Tournaments

Residue tournaments use vertices `0,...,p-1` and orient

```text
r -> s  iff  v_p(f(r)) < v_p(f(s)),
```

with residue order as the fixed Hamiltonian tie path.  The fixed-divisor
example `x^2+x+2` has `all_bad_mod_p=True` at `p=2`; the other sample rows do
not.  The small residue tournaments in this pilot are transitive because the
valuation score is scalar plus a fixed tie path.  That is a useful warning:
the first residue score detects local obstruction, but richer edge labels are
needed to expose cycles.

Sign-cube tournaments fix coefficient magnitudes and orient sign chambers by
`f(B)` at a base `B`.  As expected, the base-value tournament is transitive,
but chamber statistics are informative:

```text
magnitudes (1,1,1,1,1,1), base 5:
  64 chambers, 24 irreducible, 30 with a prime hit <=30.

magnitudes (1,1,0,1,0,1), base 5:
  16 distinct chambers, 4 irreducible, 7 with a prime hit <=30.

magnitudes (1,2,1,2,1,1), base 7:
  64 chambers, 36 irreducible, 33 with a prime hit <=30.
```

This sharpens HYP-2450's magnetization slice problem: some slices look
irreducibility-rich, but the quotient needs convolution-lift obstruction data
before it can become a theorem.

## Newton-Slope Extension

The univariate convolution grid is the one-dimensional shadow of a Newton
polytope factorization program.  Abu Salem-Gao-Lauder factorization via
polytopes suggests the higher-dimensional upgrade:

```text
Newton polytope decomposition
  -> boundary factorisation
  -> inward layer lifting through polynomial equations
  -> obstruction when no compatible interior lift exists.
```

That is exactly the same grammar as the coefficient convolution lift, but with
vertices now allowed to be boundary edges, slopes, layers, or compatibility
obligations rather than coefficients.

## Transfer To LRC14 And The 72-Code Gate

For LRC14, the analogue of `a_k=sum b_i c_j` is a resource ledger whose scalar
row says a denominator is blocked, while the hidden lift says which runners,
owners, divisor fibers, and carries realize the block.  HYP-2443/HYP-2444
already show that raw scalar denominators are too coarse; the convolution
view says to solve the hidden support lift, not merely count blocked q's.

For a self-dual `[72,36,16]` code, the weight enumerator is the coefficient
profile.  It can be scalar-feasible and still fail to be realized by binary
support/design/matroid incidence.  The code problem is therefore a
convolution-lift problem in a broader incidence algebra:

```text
weight enumerator coefficients   -> boundary totals
supports of weight-16 words       -> hidden incidence lift
Type II constraints               -> diagonal compatibility equations
nonexistence route                -> no integral/binary support lift
```

This matches HYP-2430 and HYP-2450: scalar modular positivity is not the hard
gate at length 72; support realization is.

## Assumption Challenge

Alternate vertex sets considered: factor coefficients, convolution grid cells,
diagonal sums, witness integers `m`, prime tokens of `f(m)`, residue classes,
sign chambers, Newton boundary edges, real-factor recombination subsets, LRC
denominator resources, and `[72,36,16]` support/design obligations.

Chosen quotient: hidden 2D coefficient convolution lifts.

Preserved: exact reducibility witnesses through degree 5, diagonal-sum tiling
geometry, local residue winners, and witness-value factor complexity.

Destroyed: root geometry, Galois action, large-degree factor search, and
asymptotic prime-value information.  These are side channels to reattach.

Challenged assumption: tournament vertices need not be polynomials, arcs, or
primes.  In this session the best vertices are proof obligations and hidden
factor tiles; tiebreaks are arbitrary fixed Hamiltonian paths.

## Next Moves

1. Add a bounded ILP/SAT engine for degree `>=6` convolution lifts, using the
   exact degree-4/5 solver as the regression oracle.
2. Attach Newton-slope boundary layers to sparse and multivariate examples:
   boundary factorisation first, inward quadratic constraints second.
3. Use factor-capture witnesses as search pruning: low `Omega(f(m))` gives few
   token allocations and few possible irreducible factor slots.
4. Transfer the scalar/fiber split to `[72,36,16]`: weight enumerator
   coefficients are boundary totals, while support/design incidence is the
   hidden lift that must exist.
