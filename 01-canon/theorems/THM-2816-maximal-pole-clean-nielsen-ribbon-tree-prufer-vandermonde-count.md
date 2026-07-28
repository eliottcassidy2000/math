---
id: THM-2816
title: "Maximal-pole clean Nielsen atlas via ribbon trees and Prüfer--Vandermonde counting"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  For every
  e>=1, N>=2e, and ordered positive (e+1)-part pole
  partition p, the clean genus-zero Nielsen class with zero inertia
  2^e 1^(N-2e) and totally ramified third inertia has exactly
  (e-1)! binom(N-e-1,e-1) ordered classes, independently of p.  A
  noncrossing-chord/ribbon-tree bijection and Prüfer--Vandermonde
  cancellation prove the formula.  Consequently every maximal-pole
  balanced response chamber h=e+1 is finite with this exact count.  This
  does not give explicit accessory coordinates, Keller-chart entry,
  JC(2), or DC(2).
source: root/maximal-pole-clean-nielsen-2026-07-28
depends_on:
  - THM-2796-balanced-response-stieltjes-pade-normal-form-and-one-double-zero-classification
related:
  - THM-2808-three-pole-e2-maxwell-polynomial-and-finite-accessory-classification
script: 04-computation/jc_maximal_pole_clean_nielsen_ribbon_tree_thm2816.py
output: 05-knowledge/results/jc_maximal_pole_clean_nielsen_ribbon_tree_thm2816.out
script_sha256: 0fcb7f8b9f6a197e83fa4c1818c338a4e643b387c124f4b188b8bb1c0c605f6a
output_sha256: f1a3eca7773e09a4c6affef92656c73f21191e024a1038cd58f3d83fb6f812b8
hash_basis: LF-normalized bytes
---

# THM-2816 -- maximal-pole clean Nielsen atlas via ribbon trees and Prüfer--Vandermonde counting

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

THM-2808's four-gap count is the first nontrivial member of an all-degree
tree law.  The pole lengths affect the local corner compositions, but
Prüfer degree multiplicities cancel the ribbon-order multiplicities so
completely that the final ordered count depends only on `e` and `N`.

## 1. Exact Nielsen statement

Fix integers

```text
e>=1,                       N>=2e,
p=(p_1,...,p_(e+1)),        p_i>=1,       sum_i p_i=N. (1)
```

Fix the oriented full cycle

```text
rho=(0 1 ... N-1).
```

Let `tau` have cycle type

```text
2^e 1^(N-2e),                                      (2)
```

and label the cycles of `tau rho` by `1,...,e+1`.  Count simultaneous
conjugacy classes under the centralizer `<rho>`, requiring the labelled
cycle lengths of `tau rho` to be `(p_1,...,p_(e+1))`.

Then the exact number is

```text
H_e(N;p)
 =(e-1)! binom(N-e-1,e-1)
 =(N-e-1)!/(N-2e)!.                                  (3)
```

In particular it is independent of the ordered pole partition `p`.

## 2. Maximum cycle count is exactly noncrossing

Draw the `e` transpositions of `tau` as disjoint-endpoint chords on the
oriented `N`-gon.  The associated oriented chord diagram has one vertex,
`e` edges, and `f=#cycles(tau rho)` boundary faces.  Its Euler equation is

```text
1-e+f=2-2g,

f=e+1-2g.                                             (4)
```

Equivalently, in fixed-point-safe permutation form,

```text
c(tau)+c(rho)+c(tau rho)-N=2-2g,
```

because `c(tau)=N-e` and `c(rho)=1`.

Thus `tau rho` has the required `e+1` pole cycles exactly when `g=0`,
equivalently when the chord matching is noncrossing.  This also shows why
crossing data cannot be repaired by relabelling the pole cycles: it loses
two cycles per handle.

## 3. Dual ribbon-tree bijection

The planar dual of a noncrossing `e`-chord matching is a tree `T` on the
`e+1` labelled complementary regions.  Region `i` is the labelled pole
cycle of length `p_i`.  If its tree degree is `d_i`, then:

1. its incident dual edges have `(d_i-1)!` cyclic orders; and
2. its `p_i` boundary steps split into `d_i` positive corner arcs in
   `binom(p_i-1,d_i-1)` ways.

Conversely, a labelled tree, a cyclic order at every vertex, and these
positive corner compositions reconstruct the oriented polygonal chord
diagram by the contour traversal of the ribbon tree.  The only forgotten
datum is the initial polygon vertex, exactly the rotation by `<rho>`.
Hence this is a bijection, not only a count-preserving map.

For a fixed labelled abstract tree,

```text
# realizations
 =product_i (d_i-1)! binom(p_i-1,d_i-1).              (5)
```

The positivity in `(5)` is load-bearing.  It uses distinct chord endpoints;
allowing incident transpositions would create zero-length corners and a
different passport.

## 4. Prüfer cancellation and Vandermonde

Put

```text
k_i=d_i-1>=0,                  sum_i k_i=e-1.          (6)
```

The Prüfer degree formula gives

```text
# labelled trees with degrees (d_i)
 =(e-1)!/product_i k_i!.                               (7)
```

Multiplying `(7)` by `(5)` cancels every `k_i!`.  Summing over `(6)` gives

```text
H_e(N;p)
 =(e-1)! sum_(sum k_i=e-1) product_i binom(p_i-1,k_i)
 =(e-1)! binom(sum_i(p_i-1),e-1)
 =(e-1)! binom(N-e-1,e-1),                            (8)
```

where the middle equality is multivariate Vandermonde.  This proves `(3)`.

The first three rows are

```text
e=1:  1,
e=2:  N-3,
e=3:  (N-4)(N-5).                                    (9)
```

The `e=2` row is THM-2808's three-face/four-gap atlas.  The first genuinely
ternary row `e=3` is `2,6,12,20,30,...` for `N=6,7,8,9,10`.

## 5. Consequence for maximal-pole balanced responses

In THM-2796's balanced notation, take the maximal possible pole count

```text
h=e+1.                                                (10)
```

Then the third partition is `(N)`.  Every Nielsen class above is realized
by a genus-zero rational cover with branch values normalized to
`0,infinity,1`.  After ordering the poles, put the full-cycle preimage at
infinity, send `beta_1,beta_2` to `0,1`, and normalize

```text
F=B/D,             F(infinity)=1,             B-D=A, (11)
```

where

```text
B=S E^2,                  deg E=e,
D=product_i(x-beta_i)^(p_i),
T=product_i(x-beta_i).                                (12)
```

The derivative has zeros of total order `N-(e+1)` at the poles and one
simple zero at each of the `e` double zeros.  These already exhaust
`deg D'=N-1`, so

```text
D'=N D E/T.                                           (13)
```

The constant `A` is nonzero, and with

```text
C=-NA,
G=C E/(2DT),
V=4 S D T^2/C^2,                                     (14)
```

one has

```text
F'=2G,                    F=VG^2,
2VG'+V'G=2.                                         (15)
```

Formulas `(14)--(15)` use `kappa=1`.  For a prescribed `kappa!=0`,
replace `G` by `G/kappa` and `V` by `kappa^2 V`.

Conversely, every maximal-pole balanced response has precisely this
three-value passport by THM-2796 and hence appears in the Nielsen atlas.
Therefore `(3)` is also the exact number of ordered normalized
maximal-pole response maps.

Every such response is nonsplit.  If `N>2e`, then `deg S=N-2e>0`.  If
`N=2e`, the `e+1` positive pole parts cannot all be even, so THM-2796's
odd-pole squareclass criterion applies.

## 6. What the binary/ternary motif really preserves

The rigorous common object is the oriented ribbon tree:

- order-two inertia contributes the chord edges;
- the complementary pole cycles are its vertices;
- a trivalent vertex carries a cyclic order, but this does not by itself
  create a global `C_3` action; and
- the partition data live in positive corner compositions and are erased
  by a bare abstract-tree quotient.

Thus the binary and ternary grammars co-occur on one object, but identifying
them with the free factors of `PSL_2(Z)=C_2*C_3` still requires a faithful
global action and a preserved response predicate.  THM-2808's quartic
`D_4/V_4=C_2` and THM-2800's separate `A_4/V_4=C_3` remain the sharp
co-occurrence boundary.

## 7. Exact controls

The companion independently:

1. enumerates every noncrossing chord matching and labelled pole-cycle
   orbit in declared ranges through `e=5`;
2. verifies `#cycles(tau rho)=e+1`, the Catalan raw count, every ordered
   passport, and formula `(3)`;
3. independently enumerates Prüfer words, ribbon cyclic orders, and corner
   compositions;
4. checks the multivariate Vandermonde count through `e=8`; and
5. verifies normal, optimized, and stored transcript identity.

It uses exact integer/permutation arithmetic and contains no Python
`assert` node.  Run

```text
python 04-computation/jc_maximal_pole_clean_nielsen_ribbon_tree_thm2816.py
python -O 04-computation/jc_maximal_pole_clean_nielsen_ribbon_tree_thm2816.py
```

The finite controls support but do not replace the all-degree bijection.

## 8. Scope and boundaries

This theorem gives a complete Nielsen and normalized-response count, but
for `e>=3` it does not solve the accessory coordinates by radicals or
write their multivariate Maxwell ideal.  More importantly, it does not
show that a response enters a polynomial Keller chart or satisfies the
Faber flux equations.  It proves neither `JC(2)` nor `DC(2)`.

For repeated pole parts, unmarked maps are orbits under the actual
multiplicity-stabilizer action; division by the stabilizer order can fail.
Crossing chords, incident transpositions, unlabeled poles, and
positive-characteristic wild ramification are outside `(3)`.
