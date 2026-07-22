---
id: THM-2098
title: "The mixed-torus arity wall is an exact collision budget with a vertical gap"
status: >
  PROVED from THM-2097's pair-rigidity lemma and Hunter's forest inequality.
  For a rank-two guard/terminal character restriction with n transverse
  terminals covering the radius-1/7 guard
  complement by radius-1/14 bands, the zero-intersection graph has maximum
  degree two, its positive complement is connected, total pair-collision mass
  is at least 5(n-7)/49, and every maximum spanning tree has weight at most
  the same quantity. At n=7 the zero budget is impossible, reproving the
  mixed-torus escape; n=8 is the first honest positive-collision wall, with
  exact budget 5/49. With vertical characters present, a cover cannot have
  between one and six vertical and simultaneously between one and six
  transverse bands. This isolates the live rank-8--10 guarded coefficient-
  plane residual (and gives the same abstract ledger at n=11), but does not
  discharge it or prove LRC(14). The collision budget applies only to the
  pure-transverse lane; no such budget is claimed for the low mixed lanes.
source: codex-2026-07-22-LRC-mixed-torus-arity-wall
depends_on:
  - THM-2080
  - THM-2081
  - THM-2097
related:
  - THM-2091
  - THM-2094
  - THM-2096
  - THM-2099
  - THM-2114
  - THM-2116
  - MISTAKE-239
  - MISTAKE-242
---

# THM-2098 -- the mixed-torus collision budget

Let

```text
c_0,c_1,...,c_n in Z^2,             n>=7,              (1)
```

assume the characters `c_0,...,c_n` span `Q^2`, and suppose some parameter
direction `d in Z^2` gives an odd positive guard `c_0.d` and pairwise distinct
positive terminal speeds `c_i.d`. Write

```text
C={X in T^2:||c_0.X||>1/7},
B_i={X in C:||c_i.X||<1/14}.                           (2)
```

Boundaries are null throughout. The strict convention makes every nonempty
intersection used below open and hence of positive Haar measure.

## 1. The transverse zero graph has degree at most two

Assume first that no `c_i` is rationally proportional to `c_0`. Then every
integer torus map

```text
X |->(c_i.X,c_0.X)
```

has nonzero determinant and is surjective. Therefore

```text
measure(B_i)=(1/7)(5/7)=5/49.                          (3)
```

Define the zero graph `Z` on `{1,...,n}` by

```text
ij in E(Z) iff B_i intersect B_j is empty.             (4)
```

For a zero edge, simultaneous strict danger of `c_i,c_j` implies
`||c_0.X||<=1/7`. THM-2097's pair-rigidity lemma applies and gives

```text
c_0=epsilon c_i+eta c_j,
epsilon,eta in {+1,-1}.                                (5)
```

Fix `i`. Up to sign, equation (5) leaves only the two possible terminal
characters

```text
c_0+c_i,                 c_0-c_i.                      (6)
```

No two terminal danger bands coincide. Indeed `c_j=c_k` contradicts
distinct specialized speeds, while `c_j=-c_k` contradicts their simultaneous
positivity on `d`. Hence `i` has at most two zero neighbors, and

```text
Delta(Z)<=2.                                           (7)
```

Let `P` be the complement graph, whose edges are the positive pair
intersections. Then

```text
delta(P)>=n-3,
|E(P)|>=n(n-3)/2.                                      (8)
```

For `n>=7`, the minimum degree in (8) is greater than `(n-1)/2`, so `P` is
connected. Thus the transverse band family always has a spanning tree all of
whose edge intersections have positive measure.

The faithful graph is not merely the support of pair intersections. Retain
the exact weights

```text
w_ij=measure(B_i intersect B_j).                        (9)
```

## 2. Covering creates exactly the arity-excess budget

Suppose the terminal bands cover the guard complement:

```text
C subset union_(i=1)^n B_i                             (10)
```

up to a null set.

The direction in (10) is load-bearing.  Here `C=E_(c_0)^c`, so (10) is
exactly the corrected implication in MISTAKE-239: containment of the joint
safe set in the guard forces terminal dangers to cover the **complement** of
the guard.  Nothing here asserts that the danger bands cover the guard
itself.  If the strict mixed safe cell is empty, equality boundaries form a
finite null union, which is why (10) holds almost everywhere with the strict
bands in (2).

Define

```text
N(X)=sum_i 1_(B_i)(X).
```

Equations (3) and (10) give

```text
integral_C N=5n/49,
integral_C (N-1)=5(n-7)/49.                            (11)
```

Since `N>=1` and `binom(N,2)>=N-1` pointwise,

```text
sum_(i<j) w_ij
 =integral_C binom(N,2)
 >=5(n-7)/49.                                         (12)
```

Equality holds in (12) exactly when `N in {1,2}` almost everywhere on `C`.
Thus the first excess band above seven cannot disappear: it becomes literal
pair-collision mass.

Let `tau` be the maximum total weight of a spanning tree on the `n` labels.
Hunter's pointwise forest inequality and (3) give

```text
5/7=measure(C)
 <=sum_i measure(B_i)-tau
 =5n/49-tau.
```

Consequently every cover satisfies the matching upper invoice

```text
tau<=5(n-7)/49.                                        (13)
```

Equations (12)--(13) are the exact arity wall:

```text
total pair mass >= E_n,       maximum tree mass <= E_n,
E_n=5(n-7)/49.                                         (14)
```

At `n=7`, one has `E_7=0`, while Section 1 supplies a positive-weight
spanning tree. This contradiction is an alternate weighted-graph proof of
THM-2097's transverse seven-band escape.

At `n=8`, the common budget is

```text
E_8=5/49.                                              (15)
```

It is positive, so pair rigidity by itself no longer contradicts a cover.
This is the first honest higher-rank obstruction, not a failure of the
seven-band proof technique by loose constants.

## 3. Vertical/transverse middle counts are impossible

Now allow `r` terminal columns to be rationally proportional to `c_0`, and
put `t=n-r` for the number of transverse columns. Assume

```text
1<=r<=6,                 1<=t<=6.                      (16)
```

After a unimodular coordinate change, write the guard and vertical columns
as

```text
c_0=(0,b),             c_i=(0,a_i),                    (17)
```

with `b` odd and the `a_i` positive and distinct. THM-2080 gives

```text
measure({||by||>1/7} intersect {||a_i y||<1/14})
 <=5/42,                                                 (18)
```

with equality only at `a_i=6b`. Therefore the `r` vertical danger sets do
not cover the guard-safe `y` set. For `r<=5` this follows from
`5r/42<5/7`; at `r=6`, equality would force all six distinct `a_i` to equal
`6b`.

Choose a strict-safe `y`. Each of the `t<=6` remaining transverse danger
sets has `x`-measure exactly `1/7`, so their union does not cover the
`x`-circle. The resulting `(x,y)` is safe for the guard and every terminal
character. Hence (16) is incompatible with a mixed-threshold cover.

For `n=8,9,10,11`, any rank-two mixed-threshold cover must therefore have
vertical count in the following exact outer ledger:

```text
n=8:       r in {0,1,7},
n=9:       r in {0,1,2,7,8},
n=10:      r in {0,1,2,3,7,8,9},
n=11:      r in {0,1,2,3,4,7,8,9,10}.                 (19)
```

The value `r=n` is omitted because it violates the rank-two character-span
assumption. In coefficient-plane applications, a rank-one character image is
routed separately through the existing freeze bounds. The pure-transverse
case `r=0` is governed by (14). A low mixed row
`1<=r<=n-7` has at least seven transverse bands, but (14) does **not** follow:
the transverse subfamily need not cover `C`, and the vertical events do not
have mass `5/49`. These rows are count-isolated only.

On the high side `r>=7`, one has `t<=n-7<=4` in (19). In fact the vertical
subfamily must cover the one-dimensional guard complement `C_b` up to a null
set. Otherwise its leftover `y`-set has positive measure. For every such `y`,
the remaining `t<=6` transverse bands occupy at most `t/7<1` of the
`x`-circle, so Fubini gives a positive-measure mixed safe set, contradicting
the almost-everywhere cover. Thus the high side is a genuine one-dimensional
commensurability cover modulo endpoints. The forbidden middle is completely
removed.

In the live dyadic tower, THM-2078/2080 leave guarded terminal sizes
`7,8,9,10`; thus the new LRC lanes are the `n=8,9,10` rows. The `n=11` row is
a valid abstract extension, not a guarded depth-zero LRC residual.

## 4. Frontier effect

THM-2097 makes every rank-seven depth-four coefficient template finite. The
next guarded terminal ranks eight through ten are not an amorphous failure of
its proof. They split into the following lawful lanes:

```text
pure transverse: collision/tree budget E_n=5(n-7)/49;
low mixed:       at least seven transverse bands, but no inherited E_n budget;
high vertical:   at least seven vertical bands covering C_b almost everywhere. (20)
```

THM-2091/2096's multiplicity and tree moments are tailored to the first
branch. THM-2095's prime-power deck is tailored to the third, although its
current six-comb capacity must be strengthened before it handles seven or
more vertical terminals. A conditional-on-vertical-safe-set moment is still
needed for the low mixed rows. This theorem supplies the lawful handoff; it
does not claim any higher-rank lane is empty.

## 5. Assumption challenge and Tournament Analysis

The challenged assumption is that the relevant higher-arity datum is merely
the number of bands. Equation (14) shows what the eighth band becomes: a
located pair-collision budget simultaneously constrained by a graphic-
matroid maximum. Forgetting location keeps total mass and loses the tree
upper bound; keeping only the zero graph loses both rational weights.

Candidate tournament vertices were characters, zero edges, positive edges,
multiplicity levels, and spanning-tree obligations. Orienting positive edges
by weight gives a search scheduler, but its score sequence, SCCs, cycles, and
Hamiltonian paths do not determine `tau`. The zero graph itself is a union of
paths and cycles because of (7); an arbitrary orientation adds no theorem.
The faithful carrier is the weighted positive-intersection graph together
with the multiplicity excess `E_n` and the vertical-count sidecar. QED.
