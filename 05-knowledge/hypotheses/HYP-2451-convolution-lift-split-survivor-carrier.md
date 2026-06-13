# HYP-2451 - Irreducibility is absence of a convolution lift, and split survivors are the right carrier

**Status:** OPEN synthesis; exact local scout verified.
**Source:** codex-2026-06-12.
**Companions:** HYP-2450, HYP-2449, HYP-2448, HYP-2447,
OPEN-Q-071, OPEN-Q-070, THM-474.
**Computation:** `04-computation/convolution_lift_irreducibility_carrier_codex.py`;
stored output `05-knowledge/results/convolution_lift_irreducibility_carrier_codex.out`.

## Statement

The coefficient-row tournament from HYP-2449 should be lifted one level:

```text
f(x)=a_0+a_1x+...+a_dx^d
```

is reducible exactly when the coefficient row is the diagonal shadow of a
nontrivial hidden coefficient grid

```text
a_k = sum_{i+j=k} b_i c_j.
```

So the useful carrier is not only a sign chamber, and not only a marked
coefficient row.  It is a split-survivor state:

```text
degree split r x (d-r)
+ residue-field convolution survivors
+ p-adic valuation/Newton survivors
+ fixed-divisor residues
+ evaluation-depth / recombination certificates.
```

The arbitrary fixed Hamiltonian path appears as a deterministic order for
checking split obligations and resolving ties.  It is not the source of the
structure.  The source is the hidden diagonal-sum lift.

## Evidence

The computation directly compares symbolic factorization over `F_p` with an
exact brute-force diagonal-lift enumerator in small examples:

```text
x^4 - x^3 - x^2 - x - 1 mod 2:
  factor_degrees=(4,)
  symbolic_survivors=none
  brute_survivors=none

2*x^4 - 2*x^3 - x^2 - x - 1 mod 3:
  factor_degrees=(1,1,2)
  symbolic_survivors=(1,2)
  brute_survivors=(1,2)
```

Thus "no split survives mod p" is literally "no local diagonal grid lift
exists."

For a primitive integer polynomial, a clean prime `p` with no nontrivial
factorization of `f mod p` proves irreducibility over `Z`.  In the same
degree-4 coefficient scout used by HYP-2449 (`|a_i| in {1,2,3}`, positive
leading sign):

```text
rows = 3888
irreducible = 3096
reducible = 792
least mod-p convolution blocker among p <= 31:
  3058 / 3096 irreducibles certified = 98.77%
  false positives = 0
```

The least-blocker histogram is:

```text
2: 784
3: 508
5: 744
7: 472
11: 264
13: 120
17: 76
19: 58
23: 20
29: 12
None: 38
```

This dramatically refines the sign-only quotient:

```text
key = signs:
  groups = 16
  mixed irreducibility buckets = 16

key = signs + least_modp_convolution_blocker:
  groups = 162
  mixed irreducibility buckets = 8
```

All reducible rows remain uncertified, as they must.

## Residue Versus Valuation

The important caution is that residue blockers are strong but not complete.
Some irreducible polynomials look reducible after reduction mod the same prime
that proves them by valuation.

For example:

```text
x^5 + 6*x + 3, p=3
  mod-p split survivors = (1,2)
  Newton lower hull = [(0,1),(5,0)]
  edge_gcd = 1
  certificate = True
```

and

```text
x^4 + 10*x^2 + 5, p=5
  mod-p split survivors = (1,2)
  Newton lower hull = [(0,1),(4,0)]
  edge_gcd = 1
  certificate = True
```

So the local carrier has two complementary faces:

```text
residue face:
  reduce the diagonal grid mod p and ask which split rectangles survive

valuation face:
  keep p-adic height data and ask whether the Newton hull permits a split
```

The residue face empties split sets directly.  The valuation face detects
rigidity when the residue shadow has collapsed.

## Tournament Analysis

Vertices in the printed tournament are the small local certificate channels
`2,3,5,7,11,13,17,19,23,29,31`.

Pairwise observable:

```text
p -> q iff p is the least convolution blocker for more irreducible scout rows
than q.
```

Switch/gauge:

```text
pass to the primitive part by Gauss lemma; ignore degree-drop primes.
```

Tie Hamiltonian path:

```text
2 -> 3 -> 5 -> ... -> 31.
```

Fingerprint:

```text
least_blocker_counts={2:784, 3:508, 5:744, 7:472, 11:264,
                      13:120, 17:76, 19:58, 23:20, 29:12, 31:0}
score_hist={0:1,1:1,2:1,3:1,4:1,5:1,6:1,7:1,8:1,9:1,10:1}
directed_3cycles=0
Hamiltonian_paths=1
ranking=[2,5,3,7,11,13,17,19,23,29,31]
```

The transitivity is not a proof of anything deep; it says the finite scout has
a clear low-prime certificate order.

## Why This Matters

HYP-2449 showed that the coefficient sign tiling is real but too lossy.
HYP-2451 says what should replace it:

```text
which hidden factor rectangles survive?
```

This is a better invariant than sign, score sequence, or unmarked tournament
isomorphism class.  It also absorbs several of the suggested directions:

1. **Factor-capture hypertournament.**  At a witness value `m`, prime tokens
   of `f(m)` are allocated among possible factors.  This is downstream of the
   same split-survivor state.
2. **Convolution/tiling lift.**  This is the exact hidden object:
   coefficients are diagonal totals of a 2D grid.
3. **Newton-slope tournament.**  The valuation face is the p-adic lower-hull
   version of the same lift question.
4. **Residue tournaments.**  Residue primes rank by how quickly they empty
   the split-survivor set.

## Transfer To LRC14

The LRC analogue should not ask only whether a denominator `q` is blocked.
It should ask which local lift obligations survive:

```text
denominator/resource split
+ residue fiber
+ owner/carry address
+ survivor set after local gates.
```

A dangerous LRC14 row should be one whose local survivor set never empties
under the available gates, just as a stubborn polynomial has no small
mod-p convolution blocker and needs valuation/Singh/recombination channels.

This reframes the next HYP-2443/HYP-2444 target: replace scalar blocker counts
by split-survivor ledgers over the Q27 resource fibers.

## Next Moves

1. Extend the degree-4 scout by adding p-adic Newton certificates to the 38
   irreducibles with no blocker up to `31`.
2. Build degree-5 and degree-6 split-survivor signatures, but cache local
   factorization aggressively.
3. Add Singh-depth and Cohn-depth channels only for rows that survive all
   cheap residue/valuation tests.
4. Translate `survivors_p(f)` to LRC14: for each denominator family, store the
   surviving resource splits rather than the scalar "blocked/unblocked" bit.
