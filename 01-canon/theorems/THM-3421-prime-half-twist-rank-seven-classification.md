---
id: THM-3421
title: "Prime half-twist rank-seven classification"
status: >
  RESERVED / UNPROVED PROOF CANDIDATE UNDER INDEPENDENT AUDIT.  The candidate
  classifies literal half-twist covers by at most seven blocks on a prime
  number of sheets.  Its proposed proof is elementary above p=511 and uses a
  standard-library exact normalized-mask census below that threshold.  It
  must not be used as a proved dependency until the audit and status promotion
  land.
source: root-prime-half-complete-2026-08-15
depends_on:
  - THM-3416-zero-mode-cochain-global-rank-six-support
  - THM-3420-prime-rank-seven-zero-and-half-twist-splitter-closures
related:
  - THM-3405-common-centre-gcd-gauge-and-boolean-half-twist
  - THM-3423-odd-interval-ratio-complement-and-dyadic-clique-law
script: 04-computation/lrc_prime_half_twist_rank7_classification_thm3421.py
output: 05-knowledge/results/lrc_prime_half_twist_rank7_classification_thm3421.out
script_sha256: c4193e53bb3199667681f1d9ac91fa1495d7013248988d41683922d989ee6a82
output_sha256: bf9de38f72c3e85aa6a27ecc434afbb81b9f29a878baa3ec0fd6b0e6734b0972
semantic_sha256: 9e27f37172631c8bf1133d58b769ffcbddb76e9cdbb5765865b96225db219c8c
hash_basis: LF-normalized bytes
---

# THM-3421 -- prime half-twist rank-seven classification

**RESERVED / UNPROVED PROOF CANDIDATE UNDER INDEPENDENT AUDIT.**

## 1. Candidate statement and scope

For a prime `p` and a transverse residue `r` modulo `2p`, retain the literal
half-twist mask from THM-3420,

```text
H_(p,r)={ell in Z/pZ: ||r(2ell+1)/(2p)||<1/14}.        (1)
```

The candidate theorem is

```text
some at-most-seven masks H_(p,r) cover Z/pZ
iff p is one of {11,13,23,29}.                         (2)
```

The exact minimum half-layer ranks at the four primes are respectively

```text
p:                 11, 13, 23, 29
r_(1/2)(p):          6,  7,  6,  7.                   (3)
```

Literal positive atoms are

```text
p=11: (1,2,3,5,7,9),
p=13: (1,2,3,5,7,9,11),
p=23: (1,4,5,7,9,11),
p=29: (1,5,7,8,12,13,22).                             (4)
```

The first three rows are exact sheet partitions.  The `p=29` row has
multiplicity one on 28 sheets and multiplicity three at the reflection-fixed
sheet.  All four rows pass the primitive gcd gate.  THM-3416 already proves
the rank-six claims at `11,23`, while THM-3420 proves the rank-seven claims at
`13,29`.

This is a theorem about literal Boolean half-twist masks and hence about the
corresponding prime primitive zero-mode-cochain layer.  It does **not**
classify composite rank seven, arbitrary physical common times, or LRC(14).

## 2. Inheritance and connection ledger

THM-3416 proves that a prime half-twist layer has rank at most six exactly at
`p=11,23`.  Thus, away from those two positive primes, a putative cover in
`(2)` can be assumed to consist of exactly seven essential blocks.  THM-3420
already closes the all-even layer through its fixed-zero prime theorem and
closes the `p=13 mod 14` critical layer directly.

| field | exact connection |
|---|---|
| source | seven prime-sheet half-twist danger blocks |
| target | multiplicative dilates of a short odd interval in `F_p^*` |
| map | encode a sheet by its odd word, remove the reflection-fixed word, and reduce modulo `p` |
| preserved | strict endpoints, numerator parity, block size, reflection pairs, intersection, and total overlap |
| destroyed | the original owner label and the fixed sheet after it is invoiced separately |
| required sidecar | the ratio complement `F_p^* minus O/O` |
| cheapest decisive tests | positive `p=29`, hostile `p=37`, and the first large rows `p=541,547,557,587` |

The key move is to classify the **missing ratios**, not the covered points.
This turns disjoint block packings into a small Cayley clique problem.

## 3. Capacity, parity, and the reflection-pair budget

Write

```text
p=14k+s,       s in {1,3,5,9,11,13},                  (5)
```

for `p>7`, and let `e` be the number of even-numerator blocks among seven.
The half-twist reflection is

```text
sigma(ell)=-1-ell.                                    (6)
```

It has one fixed sheet.  Every even-numerator block contains that sheet and
every odd-numerator block avoids it, so `e>=1`.  An even block has size
`2k+1`.  An odd block has size `2k` for `s=1,3,5` and size `2k+2` for
`s=9,11,13`.

Let `Omega` be total block mass minus `p`.  The fixed sheet already costs
`e-1` units of overlap, while every remaining overlap comes in reflection
pairs.  Direct substitution gives

```text
s in {1,3,5}:     Omega=e-s;
s in {9,11,13}:   Omega=14-s-e.                       (7)
```

Thus `s=3,5` are impossible because `Omega<e-1`.  For the four surviving
classes, let `B` be the number of off-fixed reflection-pair overlap units.
Then

```text
s=1:    B=0;
s=9:    1<=e<=3,   B=3-e;
s=11:   1<=e<=2,   B=2-e;
s=13:   e=1,        B=0.                              (8)
```

If `e=7`, the odd-word change of sheet identifies the seven even half blocks
with seven fixed-zero blocks.  THM-3420 therefore forces `p=29`.  Hence the
remaining argument may assume at least one odd block.

## 4. The exact odd-interval ratio complement

Represent sheets by odd words `x=2ell+1 mod 2p`, remove the fixed word `p`,
and reduce modulo `p`.  An even block has off-fixed part a dilate of

```text
E={+-1,...,+-k}.                                      (9)
```

An odd block is a dilate of

```text
O_L={x in Z: x odd and |x|<=L} subset F_p^*,          (10)
```

where

```text
s:       1   9   11  13
L:     2k-1 2k+1 2k+1 2k+1
c=p-7L:  8   2    4   6.                              (11)
```

Define the 24-element rational set

```text
D={+-a/b: gcd(a,b)=1, a,b have opposite parity,
             and a+b<=7}.                             (12)
```

We claim that for `p>=448+8c`, reduction modulo `p` gives the exact identity

```text
F_p^* minus (O_L/O_L) = D.                            (13)
```

### 4.1 Every missing ratio is short

Fix `lambda in F_p^*`.  Among the eight circular points

```text
0,lambda,2lambda,...,7lambda
```

two consecutive points have gap at most `floor(p/8)`.  Hence

```text
d lambda=+-a mod p,       1<=d<=7, 1<=a<=p/8.         (14)
```

Divide `a,d` by their gcd and rename the coprime pair `A,b`.  If `A,b` are
both odd, then `A,b<=L` and `(14)` already puts `lambda` in `O_L/O_L`.
Otherwise they have opposite parity.

Suppose `A+b>=8`.  The line

```text
b x-A y=p                                               (15)
```

has odd/odd integral solutions: because `A,b` have opposite parity, one
coordinate is automatically odd and shifting once changes the parity of the
other.  The odd/odd solutions are spaced by `(2A,2b)`.  Choose the solution
nearest the real point

```text
(x_0,y_0)=(p/(A+b),-p/(A+b)).                          (16)
```

Then

```text
max(|x|,|y|)<=p/(A+b)+max(A,b)<=p/8+8.                (17)
```

For completeness, the last inequality is not a rounding heuristic.  If
`A<=b`, it is immediate.  If `A>=b`, the function
`p/(A+b)+A` is convex on

```text
max(b,8-b)<=A<=p/8;
```

at both endpoints it is at most `p/8+8`.  Since

```text
p/8+8<=L=(p-c)/7  iff  p>=448+8c,                     (18)
```

equation `(15)` gives `lambda=x/y in O_L/O_L`.  Therefore a missing ratio
must have `A+b<=7`, so it lies in `D`.

### 4.2 Every short opposite-parity ratio is missing

Conversely, take `+-A/b in D` and suppose it equals `x/y` modulo `p` for
odd `x,y` with absolute values at most `L`.  Then one of

```text
b x-A y,       b x+A y
```

is an odd, hence nonzero, multiple of `p`, but

```text
|b x +/- A y|<=(A+b)L<=7L=p-c<p,                      (19)
```

a contradiction.  This proves `(13)`.

For the relevant classes, `c<=8`, so the uniform analytic threshold is

```text
p>=512.                                                (20)
```

## 5. The dyadic clique obstruction

Two odd-block dilates are disjoint exactly when the ratio of their scalars
lies in `D`.  Normalize one scalar to one.  Pairwise-disjoint odd blocks
therefore form a clique in the Cayley graph on `{1} union D` with adjacency

```text
u~v  iff  u/v in D.                                   (21)
```

This graph has clique number exactly three.  Indeed, every element of `D` has
two-adic valuation in `{+-1,+-2}`.  The colouring

```text
colour(u)=nu_2(u) mod 3                                (22)
```

is proper, while `{1,2,4}` is a triangle.

This rational graph is exactly the modular graph for every `p>217`.  The
exact companion checks all `25*25*24` comparisons and finds

```text
max |numerator(u/v-d)|=217
```

over unequal rational pairs `u/v,d`.  Thus a false modular equality would
force `p` to divide a nonzero integer of absolute value at most 217.  No
asymptotic or probabilistic graph statement is being used.

## 6. Spending overlap produces too many disjoint odd blocks

At any reflection orbit, if `t` odd blocks meet there, deleting `t-1` of
them removes every odd/odd collision on that orbit.  Doing this over all
orbits deletes at most the total off-fixed pair budget `B`; overlaps involving
even blocks only consume more of that same budget.  Therefore the
`7-e` odd blocks contain a pairwise-disjoint subfamily of size at least

```text
7-e-B.                                                (23)
```

For `s=9`, `(8)` makes this lower bound `4` for every `e=1,2,3`.  For
`s=11` it is `5`, and for `s=13` it is `6`.  Each contradicts the clique
bound three.

For `s=1`, the budget is zero but `(23)` alone would weaken when `e` grows.
Here any two even blocks already collide off the fixed sheet.  Indeed, for
`k>=13` one has `p<(k+1)^2`; the circular-gap argument on

```text
0,lambda,...,k lambda
```

gives `d lambda=+-a` with `1<=a,d<=k`, so `E/E=F_p^*`.  Because `B=0`, this
forces `e=1`.  The remaining six odd blocks are pairwise disjoint, again
contradicting clique number three.

Equations `(20)--(23)` exclude every prime `p>=512` outside the already
positive all-even `p=29` case.

## 7. Exact finite boundary and positive atoms

The companion performs an independent normalized-mask census for every one
of the `96` primes below `512`.  For each scalar-feasible value of `e`, common
multiplicative scaling fixes one even block, and a bit-mask DFS chooses the
remaining even and odd dilates while spending the exact overlap budget.  The
explicit universe is

```text
3<=p<512, p prime;
1<=e<=7;
197 scalar-feasible (p,e) profiles;
seven distinct normalized block dilates.              (24)
```

The positive support is exactly

```text
{11,13,23,29}.                                        (25)
```

There are twelve positive seven-block parity profiles: one at `11`, one at
`13`, three at `23`, and all seven values of `e` at `29`.  The `e=7` row at
`29` is the Paley fixed-zero splitter transported to the half layer.  Normal
and optimized runs are byte-identical.  The four literal witnesses in `(4)`
are replayed separately, including their full multiplicity profiles.

Together with THM-3416's all-prime rank-at-most-six classification, `(24)`
also covers a hypothetical shorter finite-boundary cover: away from `11,23`
it would have to be an essential seven-block cover and hence appears in the
normalized census.  This proves the finite half of `(2)` and completes the
candidate proof.

## 8. Equality boundary, losses, and non-consequences

- The analytic threshold is explicit: the worst odd-interval offset is
  `c=8`, giving `p>=512`; every smaller prime is in `(24)`.
- The clique bound is sharp at three, but the cover budget forces at least
  four disjoint odd blocks in every surviving large-prime class.
- The `s=1` lane is the only one that needs the even-block ratio sidecar;
  the other three classes close from odd-block overlap debt alone.
- The proof retains literal sheet union and strict endpoints.  It does not
  survive quotienting to density alone: all four residue classes pass scalar
  capacity.
- Composite orders can mix quotient orders and pullback degrees, so the prime
  ratio complement does not classify composite rank seven.
- The result does not lower the direct LRC(14) frontier.
