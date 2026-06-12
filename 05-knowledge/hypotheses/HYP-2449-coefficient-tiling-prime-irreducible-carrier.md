# HYP-2449 - Coefficient tiling is a real tournament carrier, but irreducibility needs marked row addresses

**Status:** OPEN synthesis; finite coefficient-tiling scout verified.
**Source:** codex-2026-06-12.
**Companions:** HYP-2447, HYP-2448, OPEN-Q-070, HYP-2446,
HYP-2443, THM-474.
**Computation:** `04-computation/coefficient_tiling_prime_irreducible_codex.py`;
stored output `05-knowledge/results/coefficient_tiling_prime_irreducible_codex.out`.

## Statement

The user's coefficient-tiling picture is literally present in the fixed-path
tournament tiling cube.

For a degree `d` polynomial

```text
f(x)=a_0+a_1 x+...+a_d x^d,
```

the skip-row decomposition of a linearly ordered tournament has row sizes

```text
d, d-1, ..., 2, 1.
```

Thus for degree `5`, a row-constant sign tournament can place

```text
a5 on the apex row of size 1,
a4 on the row of size 2,
a3 on the row of size 3,
a2 on the row of size 4,
a1 on the row of size 5.
```

There are two natural versions.

```text
top_no_constant:
  degree d -> tournament on d+1 vertices;
  a_s orients all edges of skip s, for s=1..d;
  a_0 remains outside as a scalar/evaluation anchor.

constant_spine:
  degree d -> tournament on d+2 vertices;
  a_0 orients the adjacent Hamiltonian-path row;
  a_d is the apex row.
```

The second version is the stronger one.  It turns the idea around: the
constant coefficient is not external to the tournament.  It is the spine, the
evaluation-at-zero gate, and the first fixed-divisor witness.

However, the bare unmarked tournament is far too coarse.  The proof carrier is
not just the sign tournament, but

```text
marked skip-row tiling
+ coefficient magnitudes
+ local residue / valuation addresses
+ evaluation-depth certificates
+ recombination survivor data.
```

In short: coefficient signs give a real tiling cube; irreducibility and prime
production live in the addresses that survive projection.

## Evidence

The script exhaustively checks degree-4 polynomials with `|a_i| in {1,2,3}`
and positive leading sign:

```text
rows = 3888
irreducible = 3096
reducible = 792
fixed_divisor > 1 = 800
no prime hit on n=0..40 = 552
```

The unmarked coefficient tournaments leak every tested arithmetic predicate:

```text
top_no_constant canonical key:
  groups = 6, mixed irreducibility buckets = 6

constant_spine canonical key:
  groups = 9, mixed irreducibility buckets = 9

score sequence only:
  groups = 4, mixed irreducibility buckets = 4
```

Even the marked sign vector is not enough:

```text
marked_signs:
  groups = 16, mixed irreducibility buckets = 16
```

The same coefficient-sign tournament can contain both irreducible and
reducible rows:

```text
signs = (-,-,-,-,+)

x^4 - x^3 - x^2 - x - 1
  irreducible, fixed divisor 1, prime_hits=11

2*x^4 - 2*x^3 - x^2 - x - 1
  factor degrees (2,2), fixed divisor 1, prime_hits=0
```

The same signs can also hide fixed-divisor obstruction:

```text
x^4 - x^3 - x^2 - x - 1
  fixed divisor 1

x^4 - x^3 - x^2 - x - 2
  fixed divisor 2, local zero prime (2)
```

Adding the fixed-divisor residue address removes that particular leakage in
the scout:

```text
marked_signs + local_zero_primes:
  groups = 48, mixed fixed-divisor buckets = 0
```

This mirrors HYP-2446/HYP-2448: scalar or sign shadows are not proof objects
until the local obstruction channel is retained.

## Cohn Warning

Cohn digit polynomials show why sign-only coefficient tournaments cannot be
the final answer.

Both of these binary/base-3 repunit rows have all centered digit signs
positive, hence a transitive sign tournament with `H=1`:

```text
9841 in base 3 -> 1+x+...+x^8
  Omega(9841)=2, factor degrees (2,6), reducible.

2047 in base 2 -> 1+x+...+x^10
  Omega(2047)=2, factor degrees (10), irreducible.
```

The sign tournament is identical in spirit, but the arithmetic differs.  Cohn
does not live in signs alone; it lives in the place-value address.

## Procedural Carrier Tournament

The script compares candidate carriers:

```text
row_sign_plus_local_residues
newton_polygon_valuation_tiling
constant_spine_tiling
cohn_digit_carry_address
marked_skip_row_tiling
singh_value_depth_state
trace_subset_recombination
coefficient_vertex_tournament
raw_unmarked_skip_tournament
```

The majority tournament is transitive in this first scoring:

```text
score_hist = {0:1, 1:1, 4:3, 5:1, 7:1, 8:2}
directed_3cycles = 0
Hamiltonian paths = 84
leaders = row_sign_plus_local_residues,
          newton_polygon_valuation_tiling
```

This says the strongest next carrier is not merely the user's sign tiling, but
the sign tiling with local obstruction and valuation data attached.

## New Directions

1. **Gauss lemma as switching.**  Content is a global gauge; passing to the
   primitive part is the quotient where row-gcd debt has been removed.
2. **Eisenstein as a valuation tournament.**  A prime `p` orients lower
   coefficient rows into a p-adic sink, with the leading row as escape and the
   constant row as the `p^2` guard.
3. **Newton polygon tiling.**  Replace sign rows by p-adic valuation rows and
   orient edges by slope comparisons.  Irreducibility criteria become lower
   hull rigidity statements.
4. **Cohn as exponential row weighting.**  Base `b` gives skip rows a
   place-value address; primality says no nontrivial row subset recombines at
   that address.
5. **LRC14 fixed-divisor analogue.**  A denominator/resource family is
   dangerous only if some runner/resource row vanishes on every local residue,
   the same way `x^2+x+2` is always even.
6. **A000568/switching-class transfer.**  THM-474 already says the fixed-path
   tiling cube is the switching-class coordinate system.  Polynomial
   coefficient signs should be studied as marked switching classes, not as
   unmarked tournament isomorphism classes.

## Assumption Challenge

Alternate vertex sets considered: coefficient tiles, skip rows, coefficients,
monomials, polynomial sign patterns, coefficient magnitudes, local residue
classes, p-adic valuation rows, digit positions, evaluation times, factor
subsets, roots, Newton polygon edges, LRC runners, LRC denominator resources,
switching classes, and proof obligations.

Chosen first carrier: marked skip rows plus residue/valuation/evaluation
addresses.

Preserved: coefficient row geometry, fixed-path tiling structure, local
obstruction detection, and transfer to Cohn/Singh/recombination/LRC ledgers.

Destroyed by the bare quotient: coefficient magnitudes, local residue
obstructions, place-value carries, trace subset ambiguity, and actual
asymptotic prime production.

Challenged assumption: if a tournament appears in this bridge, its vertices
need not be polynomials or primes.  The best finite tournament may live on
coefficient rows, valuation slopes, or certificate states.

## Next Moves

1. Build exact p-adic Newton-row tournaments for Eisenstein, Dumas, and
   Perron-style irreducibility families.
2. Extend the coefficient sweep to degree `5` and `6` with cached factor
   profiles, then measure edge flips as local residue primes are added.
3. Turn Cohn's base-`b` digit address into a weighted skip-row tournament and
   compare prime `N` versus composite low-Omega `N`.
4. Translate the fixed-divisor row detector to LRC14 Q27 resource rows:
   a blocked denominator should be a local zero row, not a scalar count.
