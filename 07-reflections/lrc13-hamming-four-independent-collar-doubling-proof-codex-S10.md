---
theorem: THM-815
artifact_type: independent proof companion
title: Independent collar/doubling proof of the scale-one Hamming-four closure
status: PROVED (collar-cycle classification, reciprocal finite reduction, and recursive doubling) + FINITE-EXACT (768,735 component-containment rows)
source: codex-2026-07-15-S10 Hamming-four continuation
depends_on:
  - LRCUpTo13  # only the settled eight-speed bound M>=1/9
  - THM-806    # universal collar and the symbolic two-/three-cycle exclusions
related: [THM-770, THM-795, THM-800, THM-804, THM-810, THM-816, HYP-6820]
verification:
  - 04-computation/lrc13_scale_one_hamming_four_collar_closure_codex_S10.cpp
  - 05-knowledge/results/lrc13_scale_one_hamming_four_collar_closure_codex_S10.out
---

# THM-815 independent proof companion — collar/doubling closure

The canonical theorem statement and the discrepancy/component proof are in
`THM-815-scale-one-hamming-four-safe-component-closure.md`.  This companion
retains a separately derived collar/doubling proof and its independent C++
certificate; it deliberately has no second theorem ID.

Put `delta=1/13` and `[12]={1,...,12}`.

## Theorem

Let `R` be a four-element subset of `[12]`, and for every `r in R` let
`h_r>=1`.  Then the proper scale-one four-replacement packet

```text
B=([12] minus R) union {u_r:r in R},
u_r=r+13h_r,                                             (1)
```

satisfies

```text
M(B)>1/13.                                               (2)
```

Consequently the common-scale alternative in THM-810 is empty.  Every
hypothetically tight four-replacement packet around a unit AP dilation would
have to lie on THM-810's all-order-three coset interface.  THM-816, proved
independently by a dynamic residual-comb recursion, closes that structured
`s=3` interface.  Thus THM-815 and THM-816 together close the entire
residue-preserving Hamming-four star at arbitrary height and AP scale.

## 1. The owner collar and the exact arrow convention

Write

```text
P=[12] minus R.
```

For an owner label `r`, choose `a_r r=1 (mod 13)` with `1<=a_r<=12` and put
`tau_r=a_r/13`.  The proof of THM-806's universal collar uses only the fact
that retained core speeds lie in `[12]`, not that there are nine of them.
It therefore gives here as well

```text
(tau_r-1/156,tau_r)
  subset {t:min_(p in P)||pt||>1/13}.                    (3)
```

If `u_r>24`, the left endpoint

```text
x_r=tau_r-2/(13u_r)                                     (4)
```

of the owner's aligned danger tooth lies strictly inside (3).  Immediately
to the left of `x_r`, the own replacement is safe.  Under hypothetical
tightness some other replacement must cover that left germ.

The convention in this theorem is:

```text
q -> r       means provider q covers owner r at x_r,
z=q r^(-1)  in F_13^*,
lambda=u_q/u_r.                                         (5)
```

The exact half-open handoff predicate is

```text
-1 < z-2lambda-13m <= 1                 for some m in Z, (6)
```

or, after putting `k=z-13m`,

```text
(k-1)/2 <= lambda < (k+1)/2.                            (7)
```

In particular, if `lambda<1`, then necessarily

```text
z=2,                         1/2<=lambda<1.              (8)
```

THM-806 proves directly from these same bands that a handoff digraph has no
directed two-cycle and no directed three-cycle.  Those proofs involve only
the vertices on the alleged cycle and therefore apply unchanged inside the
four-label graph.

## 2. An all-large row has one forced four-cycle

Assume first that every replacement exceeds `24`.  By (3)--(6), every owner
has positive cross-handoff indegree.  A finite loopless digraph with positive
indegree contains a directed cycle.  The two- and three-cycle exclusions
leave a directed four-cycle.

For its four directed edges let `(lambda_i,z_i,k_i)` have the meaning in
(5)--(7).  Going once around the cycle gives

```text
product lambda_i=1,                    product z_i=1 mod 13. (9)
```

Every `lambda_i>=1/2`: this is (8) below one and is automatic above one.
Equation (9) then also gives `lambda_i<=8`.  Consecutive cycle labels are
distinct, so `z_i!=1`.  The complete possible integer-centre bank is

```text
k_i in {2,3,...,12,15,16,17}.                           (10)
```

Multiplying the lower and upper inequalities in (7) yields

```text
product(k_i-1)<=16<product(k_i+1),
product k_i=1 mod 13.                                   (11)
```

### Finite band lemma

Among the `C(17,4)=2380` nondecreasing four-multisets drawn from (10), the
unique multiset satisfying (11) is

```text
{2,2,2,5}.                                              (12)
```

This is a tiny exact integer lemma; the replay enumerates all 2,380 rows.
It is also readily hand-checkable after observing from the first inequality
that `product(k_i-1)<=16`.

With the provider-to-owner convention (5), normalize the residue-five edge
to start at `a`.  The cycle and missing-label packet are

```text
a -> 8a -> 4a -> 2a -> a,
z-word       (5,2,2,2),
R=a{1,2,4,8}.                                           (13)
```

The cyclic rotation `(2,2,2,5)` is of course the same word.  The alternative
display `a{1,2,4,7}` describes the same twelve-set multiplicative orbit,
because

```text
2{1,2,4,7}={1,2,4,8} mod 13.                            (14)
```

Thus there is no convention-dependent second packet.

## 3. The exceptional four-cycle is uniformly finite

Write the four speeds around (13) as `(A,B,C,D)`, with the residue-five edge
`A -> B`.  Their exact bands are

```text
2B<=A<3B,
C<=2B<3C,
D<=2C<3D,
A<=2D<3A.                                               (15)
```

Every pair of the four speeds has ratio at most four.  Adjacent pairs follow
immediately from (15); for the two diagonals use, for example,
`A/C=(A/D)(D/C)<=4` and
`D/B=(D/C)(C/B)<=4`.  Hence, with `x=min(A,B,C,D)`,

```text
max(A,B,C,D)<=4x.                                       (16)
```

There is also an absolute bound on `x`.  Let `m=max(P)<=12`.  The settled
eight-runner Lonely Runner bound supplies a point at which every member of
`P` has clearance at least `1/9`.  The Lipschitz interval on which all eight
are strict-safe at `delta` has length

```text
L=2(1/9-1/13)/m=8/(117m).                               (17)
```

For one replacement speed `u`, THM-806's periodic-danger estimate says

```text
meas(I intersect D_u(delta))<=2L/13+2/(13u).             (18)
```

If the four danger combs covered the interval in (17), summing (18) would
give

```text
sum_(r in R) 1/u_r >= (5/2)L=20/(117m).                 (19)
```

The residue-five edge in (15) has one speed at least twice another.  Since
every speed is at least `x`,

```text
sum_(r in R) 1/u_r <= 3/x+1/(2x)=7/(2x).                (20)
```

Combining (19)--(20),

```text
x<=floor(819m/40)<=245.                                 (21)
```

Equations (13), (16), and (21) are a genuinely small uniform finite box for
the only all-large collar survivor.

## 4. Exact rejection of every all-large survivor

For each of the twelve label sets in (13), enumerate the positive lifts
satisfying (15), every speed greater than `24`, (21), and (16).  Sorting the
four replacements as

```text
x<v<w<z,                                                 (22)
```

put

```text
Q=P union {x,v,w},
E_Q={t:min_(q in Q)||qt||>1/13}.                         (23)
```

The strict-safe bands of one speed `q` are

```text
((13j+1)/(13q),(13(j+1)-1)/(13q)),        0<=j<q.        (24)
```

Their exact intersections give every component `(l,h)` of `E_Q`.  With
`c=(l+h)/2` and `eta=(h-l)/2`, the last danger comb contains this component
exactly when

```text
||zc||+z eta<=1/13.                                     (25)
```

The replay clears the endpoint denominators in (25); its hot loop uses only
integers.  It obtains

```text
all-large rows                         626,962
distinct eleven-speed cores             64,404
rows satisfying every (25)                    0.         (26)
```

The closest first failed containment has labels `(3,4,6,8)`, replacements
`(188,302,550,593)`, and component

```text
(14/117,857/7150).
```

Its denominator-cleared relative surplus is `1/550`, equivalently the
unscaled margin in (25) is `1/7150`.

Thus a hypothetical tight row cannot be all-large:

```text
min_(r in R) u_r<=24.                                   (27)
```

## 5. The collar recursively forces a doubling box

Order the four replacements as in (22).  From properness and (27),

```text
14<=x<=24.                                               (28)
```

The subunit rule (8) now recursively bounds the other three replacements.

First, if `v>2x`, then all of `v,w,z` exceed `24`, and `x` is too small to
handoff to any of their owner exits.  The three large vertices would have
positive indegree inside their induced digraph, forcing a forbidden directed
two- or three-cycle.  Hence

```text
v<=2x.                                                   (29)
```

Next, if `w>2v`, then `w,z>24`, and neither `x` nor `v` can reach either
owner exit.  The only possible handoffs are `z -> w` and `w -> z`, a
forbidden two-cycle.  Hence

```text
w<=2v.                                                   (30)
```

Finally, if `z>2w`, every other replacement has provider/owner ratio below
`1/2`, so the large owner `z` has no incoming handoff at all.  Hence

```text
z<=2w.                                                   (31)
```

The complete residual is therefore the tiny symbolic box

```text
14<=x<=24,       x<v<=2x,       v<w<=2v,       w<z<=2w. (32)
```

## 6. Exact rejection of the anchored box

For each missing-label quadruple, the least replacement `x<=24` is unique
and equals its label plus `13`.  Enumerate the remaining three residue lifts
under (32), form the eleven-speed core in (23), and apply the same exact
component-containment test (25).  The result is

```text
anchored rows                           141,773
distinct eleven-speed cores              38,196
rows satisfying every (25)                    0.         (33)
```

The closest first failed containment has missing labels `(2,4,7,9)`, ordered
replacements `(22,41,82,111)`, and component

```text
(7/65,58/533).
```

Its relative surplus is `1/41`, or `1/533` in (25).  Together, (26) and
(33) certify

```text
768,735 exact rows,                  zero tight rows.    (34)
```

As an independent low-height check, exact piecewise-linear maximin evaluation
of all `C(12,4)=495` height-one rows has minimum

```text
1/11,
```

uniquely for missing labels `(1,3,5,7)`.

The full packet is a complete nonzero residue transversal modulo `13`, so it
already has margin exactly `1/13` at every nonzero thirteenth.  Consequently
failure of even one containment (25) leaves a nonempty strict witness
interval and proves `M(B)>1/13`.  Equations (26), (33), and (34) prove the
theorem. ∎

## 7. Descent from THM-810's common-scale alternative

For a four-replacement packet around `c[12]`, THM-810 defines

```text
D_r=c/gcd(c,w_r).
```

In its all-order-one alternative, `D_r=1` gives `c|w_r` for all four
replacements.  Write `w_r=cu_r`.  Because `c` is a unit modulo `13`, the
residue condition gives `u_r=r (mod 13)`.  Any unchanged coordinates reduce
to the already closed Hamming-one, -two, or -three stars; if all four are
proper, the normalized packet is exactly (1).  Multiplication by `c` is
onto on the circle, so

```text
M(cB)=M(B)>1/13.
```

Therefore THM-810's common-scale radius-four branch is empty.  Its other
alternative is the order-three coset packet identified there as an `s=3`
deep interface; THM-816 now proves that alternative empty as well.  The two
theorems are logically complementary: THM-816 is not used in the proof of
the present scale-one closure.

## 8. Tournament Analysis and assumption challenge

The theorem-bearing local carrier is not a tournament on runners.  It is the
bipartite incidence between owner-collar exit obligations and provider danger
combs, decorated by `(z,lambda)`, the half-open boundary flag, and the owner
endpoint.  At finite depth the carrier changes to strict-safe components
versus last-speed danger teeth, decorated by component width as in (25).

For telemetry, use the method-limit row

```text
R={1,2,4,8},          (u_1,u_2,u_4,u_8)=(79,54,30,34).  (35)
```

It has exactly the live handoff cycle

```text
1 -> 8 -> 4 -> 2 -> 1
```

and silent diagonals `{1,4}`, `{2,8}`.  The pairwise observable is the
antisymmetric left-handoff difference.  Complete silent pairs first by
increasing label and then by decreasing label as the switch.  The shared
tie Hamiltonian path is `(1,8,4,2)`.  Both gauges have

```text
score histogram                     {1:2,2:2}
directed triangles                           2
SCC sizes                                  [4]
Hamiltonian paths                            5
edge flips under the switch                  2.          (36)
```

Yet (35) is loose with exact `M=3/19`.  Thus the live cycle is a finite-box
reduction, not the LRC predicate.

- Runner vertices with the handoff observable preserve a pair decision but
  destroy simultaneous owner-exit coverage.
- Residue vertices preserve the forced packet (13) but destroy `lambda`,
  hence both (21) and the recursive doubling box.
- Gaps or fixed circle sections destroy the moving tooth scale.
- Bare wall-crossing events retain order but lose owner identity and width.
- Fourier modes retain average danger mass but not the connected component
  whose containment is tested in (25).
- Bare proof-obligation vertices forget whether their witness intervals are
  compatible on the circle.

The challenged assumption is therefore explicit: tournament vertices need
not be runners or arcs.  Owner exits and strict-safe components are the two
correct vertex sorts at the two recursive stages; the tournament shadow is
useful only after their typed incidence and metric stalk have been retained.

## Reproduction and certificate digests

```bash
c++ -O3 -std=c++17 \
  04-computation/lrc13_scale_one_hamming_four_collar_closure_codex_S10.cpp \
  -o /tmp/lrc13_h4
/tmp/lrc13_h4
```

For each row, the replay hashes sixteen big-endian unsigned 64-bit fields:
four sorted labels, four ordered speeds, the first failed component index,
the cleared surplus, denominator and left side, and the two unreduced
component endpoints.  The branch digests are

```text
all-large  27c45d31f19370b8b3c30e79f378b5b3ed9b1f9538062ac2f80e7dd056a6a64e
anchored   07594ab0e69196583fdf667b4d54c8a048a1b4d2b2a87924d26a7da4d8bc7542
```

Independent `-O3` and `-O0` builds produce byte-identical stored output with
SHA-256

```text
f098acc358f534f4edf75e1affcaa03ff0bf9cda83f058d5fc86cfc984d2dca0.
```
