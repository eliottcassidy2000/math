---
id: THM-962
title: Scale-seven Hamming-six product-character and square-sum obstruction
status: PROVED STRUCTURAL + FINITE-EXACT + LEAN TERMINAL BRICK GREEN — all 108,673,488 primitive proper AP-centred c=7 leave-one-out order/unit contexts fail necessary common-sheet coverage; capacity forces all six effective orders to seven, a row-product character kills 922 of 924 supports, and a mod-seven square-sum contradiction kills the two quadratic-coset exceptions; the final square-sum obstruction is kernel-checked in Lean with no sorry/native_decide
source: codex-2026-07-17-S63 scale-seven product/square audit
depends_on: [THM-765, THM-810, THM-823, THM-859, THM-860]
related: [THM-862, THM-957, THM-958, THM-960, THM-963, THM-969, THM-970, THM-974, HYP-6820]
verification:
  - 04-computation/lrc13_scale_seven_hamming_six_product_obstruction_codex_S63.py
  - 05-knowledge/results/lrc13_scale_seven_hamming_six_product_obstruction_codex_S63.out
  - 04-computation/lean/TournamentH7/TournamentH7/LRCScaleSevenSquareSum.lean
---

# THM-962 — scale seven is killed by a product character and a square sum

Let `R subset F_13^*` have six elements, put `P=[12] minus R`, and consider

```text
A=7P union {w_r:r in R},
w_r=7r (mod 13),              w_r>0,              w_r!=7r.       (1)
```

If `M(A)<=1/13`, THM-860 supplies common-sheet coverage at every replacement
owner and the hereditary leave-one-out law for

```text
D_r=7/gcd(7,w_r),             e_r=(D_r w_r/7) mod D_r.           (2)
```

Then no such primitive proper packet exists.  Equivalently, the complete
AP-centred common-scale-seven Hamming-six sheet bank is empty, so every packet
in (1) has `M(A)>1/13` and no metric component recursion is born.

## 1. Capacity forces the all-order-seven stratum

Since seven is prime, the effective states are

```text
(D,e)=(1,0),(7,1),(7,2),(7,3),(7,4),(7,5),(7,6).                 (3)
```

The leave-one-out lcm law says that at least two coordinates have order seven.
If exactly `k` coordinates do, choose an order-seven owner.  Literal CRT
reduction gives:

- an order-one provider fills all seven sheets at its own owner and no sheet
  at a different owner;
- an order-seven owner fills exactly two sheets at itself;
- every other order-seven provider fills exactly one sheet at that owner.

Thus this owner sees at most

```text
2+(k-1)=k+1                                                     (4)
```

sheets.  Coverage forces `k=6`.  The exact state-word and labelled-context
strata are

| number `k` of order-seven colours | state words | labelled contexts |
|---:|---:|---:|
| 2 | 540 | 498,960 |
| 3 | 4,320 | 3,991,680 |
| 4 | 19,440 | 17,962,560 |
| 5 | 46,656 | 43,110,144 |
| 6 | 46,656 | 43,110,144 |

Hence (4) eliminates `65,563,344` of the total

```text
924*(7^6-1-6*6)=108,673,488                                   (5)
```

contexts and leaves `924*6^6=43,110,144` all-order-seven contexts.

## 2. The exact local object is a multiplicative Latin row

Fix an owner `o`.  The six self masks share one fixed sheet.  Use the other
sheet in the self mask of unit `e` to label the six nonfixed sheets by
`e in F_7^*`.  This is an owner-dependent sheet gauge, not a circle symmetry.
In that gauge an off-owner provider `r` of unit `e_r` fills the singleton

```text
lambda(r/o)e_r in F_7^*,                                      (6)
```

where, for the representative `a in {1,...,12}`,

```text
lambda(a)=ceil(a/2)^(-1) mod 7,

ratio a       {1,2} {3,4} {5,6} {7,8} {9,10} {11,12}
lambda(a)        1     4     5     2      3       6.           (7)
```

The self provider contributes its fixed sheet and the symbol `e_o`, which is
also (6) because `lambda(1)=1`.  Since the remaining five masks are
singletons outside the fixed sheet, owner `o` is covered exactly when

```text
{lambda(r/o)e_r:r in R}=F_7^*.                                (8)
```

Thus the faithful sheet object at scale seven is the `6 by 6` character
matrix

```text
Lambda_R=(lambda(r/o))_(o,r in R),                             (9)
```

with column scaling by the unit word.  Common-sheet coverage asks that every
scaled row be a permutation of `F_7^*`.

## 3. A row-product character leaves only two supports

Put

```text
E=product_(r in R)e_r,
G_o=product_(r in R)lambda(r/o).                               (10)
```

If (8) holds, Wilson's theorem in `F_7` gives

```text
E G_o=product_(x in F_7^*)x=-1.                               (11)
```

The left factor `E` is independent of `o`, so all six `G_o` must agree.  An
exact scan of the `binom(12,6)=924` supports gives the maximum multiplicity
distribution of the values `(G_o:o in R)`:

```text
maximum multiplicity     2    3    4    5    6
number of supports      438  328  132   24    2.               (12)
```

The two constant-product supports are exactly the quadratic residues and
their nonsquare coset:

```text
Q ={1,3,4,9,10,12},
NQ={2,5,6,7,8,11}.                                            (13)
```

Therefore (11) kills `922*6^6=43,016,832` unit contexts without inspecting
their units.  In particular all 64 signed-doubling six-cycle supports that
organized scales four through six die at this product stage.  Scale seven is
not another refinement of their affine nerve: its hard supports have changed.

The classification in (12)--(13) is finite-exact.  The verifier constructs
the signature directly from (7) on every support and freezes the complete
924-row payload.  No floating-point or randomized step occurs.

## 4. The two quadratic cosets fail by square sums

It remains to exclude (13).  The two cases are multiplicative translates, so
write either support as `alpha Q` and cyclically order it by

```text
q_i=4^i mod 13=(1,4,3,12,9,10),             i in Z/6Z.         (14)
```

Along this cycle the squared character row is

```text
(lambda(q_0)^2,...,lambda(q_5)^2)=(1,2,2,1,2,2) in F_7.       (15)
```

Put `z_i=e_(alpha q_i)^2` and `S=sum_i z_i`.  If owner `alpha q_i`
is covered, its six symbols form `F_7^*`; hence their squared sum is

```text
sum_(x in F_7^*)x^2=0.                                       (16)
```

Using (15), equation (16) is exactly

```text
2S-z_i-z_(i+3)=0,                         i in Z/6Z.            (17)
```

Summing the six equations gives `10S=3S=0`, so `S=0`.  Substitution
in (17) gives

```text
z_(i+3)=-z_i.                                                 (18)
```

But every `z_i` is a nonzero square modulo seven.  The squares are
`{1,2,4}`, while their negatives are `{6,5,3}`.  They are disjoint, so (18)
is impossible.  This excludes both supports in (13) and proves the theorem.

As an independent literal check, the verifier replays all
`2*6^6=93,312` unit words on the raw seven-sheet masks.  Each exceptional
support has the exact owner-satisfaction profile

```text
owners satisfied       0       2      4
unit words          44,712   1,728    216,                         (19)
```

and no other satisfaction count occurs.

## 5. Tournament audit and challenged vertices

Runner vertices do not retain (10) or (16), but there is a natural pairwise
telemetry quotient.  For `a<b` take the directed pair observable

```text
kappa(a,b)=lambda(b/a)/lambda(a/b) in F_7^*.                  (20)
```

Write `kappa=3^d`.  The binary switch orients `a -> b` for
`d in {1,2}`, orients `b -> a` for `d in {4,5}`, and declares the
self-inverse cases `d in {0,3}` tied.  Ties are oriented forward along the
increasing-label Hamiltonian path.  Across all 924 supports the completion
has

```text
joint fingerprints             245
directed triangles             0,...,8
SCC profiles                   (6), (5,1), (4,1,1), (3,3),
                               (3,1,1,1), (1,1,1,1,1,1)
tie-edge counts                7,...,12,15
Hamiltonian-path counts        1,...,45 over 16 attained values.            (21)
```

The telling case is (13): every one of its fifteen pairs is tied, so both
hard supports become the same transitive tournament, with score multiset
`(0,1,2,3,4,5)`, zero directed triangles, six singleton SCCs, and one
Hamiltonian path.  Pairwise tournament completion is therefore maximally
uninformative exactly where the proof needs the six-fold row product and the
cyclic square convolution.

The assumption challenge is substantive.  Runner vertices and pair edges
preserve only ratios two at a time.  Owner rows of (9) preserve the exact
common-sheet predicate; their product characters and power sums expose the
obstruction.  Passing from the matrix to a tournament destroys both.  The
faithful information ladder at scale seven is

```text
literal 42 owner-sheet obligations
  -> six multiplicative Latin rows in F_7^*
  -> row products G_o
  -> the quadratic-coset square convolution.                  (22)
```

## 6. Cross-scale meaning and scope

The scale-four through scale-six story was organized by 64 recurring
signed-doubling cycles and increasingly restrictive owner-obligation nerves.
At scale seven every one of the 924 supports is owner-locally feasible: each
single all-different row admits `6!` unit words.  The obstruction becomes
global and algebraic instead.  A degree-six product moment discards 922
supports, and a degree-two power sum kills the two subgroup cosets on which
that product moment is blind.  The change of carrier is the useful structural
lesson, not merely the empty finite bank.

This theorem closes only the primitive proper AP-centred common-scale-seven
Hamming-six face.  It does not close `c>=8`, the smooth-ramified H5 metric
bank, non-AP-centred/deep-sheet packets, radius-seven continuum-to-grid
issues, or global `n=12` sporadic emptiness.

## Verification

```bash
python3 \
  04-computation/lrc13_scale_seven_hamming_six_product_obstruction_codex_S63.py |
  cmp - \
  05-knowledge/results/lrc13_scale_seven_hamming_six_product_obstruction_codex_S63.out

python3 -O \
  04-computation/lrc13_scale_seven_hamming_six_product_obstruction_codex_S63.py |
  cmp - \
  05-knowledge/results/lrc13_scale_seven_hamming_six_product_obstruction_codex_S63.out
```

The verifier independently performs all of the following:

1. reconstructs every local mask from its CRT definition and checks (6)--(8)
   at all twelve owners, ratios, and units;
2. enumerates all `7^6` state words and verifies the five admissible strata;
3. classifies all 924 support product signatures and freezes their payload;
4. checks the square-sum contradiction on all `3^6` square words; and
5. replays the 93,312 exceptional unit contexts on both the symbolic rows and
   literal sheet masks.

The production Lean module
`TournamentH7/LRCScaleSevenSquareSum.lean` independently proves the terminal
algebraic implication used in Section 4.  It first proves that the three
opposite-pair equations force `S=0` in `ZMod 7`, then combines one resulting
opposite-pair zero sum with a kernel-`decide` proof that two nonzero squares
modulo seven cannot sum to zero.  This module does not claim the preceding
literal-mask-to-character reduction, which remains covered by the frozen
finite-exact certificate.

Frozen SHA-256 values:

```text
4a478d69444f31e52ae5168adb87d4646da091bc98d62ee37cd17e2e8e5b006e  Python source
e547e9df128af622b04eea0f34a4eb5d71500e87ff41134cbb32738d10249917  output
```
