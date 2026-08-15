---
id: THM-3421
title: "Prime half-twist rank-seven classification"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.  Literal half-twist covers
  by at most seven blocks on a prime number of sheets exist exactly for
  p in {11,13,23,29}, with exact minimum ranks 6,7,6,7 respectively.  The
  proof is elementary above p=511 and uses a standard-library exact
  normalized-mask census below that threshold.  Composite sheet numbers,
  arbitrary physical common times, and LRC(14) remain outside the scope.
source: root-prime-half-complete-2026-08-15
audit: independent proof reconstruction; literal-mask bank conjugacy and direct optimization census; normal/optimized/stored-output replay; hash and documentation audit clean
depends_on:
  - THM-3416-zero-mode-cochain-global-rank-six-support
  - THM-3420-prime-rank-seven-zero-and-half-twist-splitter-closures
  - THM-3423-odd-interval-ratio-complement-and-dyadic-clique-law
related:
  - THM-3405-common-centre-gcd-gauge-and-boolean-half-twist
script: 04-computation/lrc_prime_half_twist_rank7_classification_thm3421.py
output: 05-knowledge/results/lrc_prime_half_twist_rank7_classification_thm3421.out
script_sha256: 693db87db8664f76fdc7bc6d520c8d3d36c51376aa6c9a67beef1aad6d84aae9
output_sha256: 77ad7cf9566f5aa02f9bfa7d9abbbce21b4a37db7137616b2771f1b8169775dc
semantic_sha256: f0c7ceb122e90fd3792d72e80a205a6c4bd823cfb5b04064764b95ddc9d9471d
hash_basis: LF-normalized bytes
---

# THM-3421 -- prime half-twist rank-seven classification

**PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.**

## 1. Statement and scope

For a prime `p` and a transverse residue `r` modulo `2p` (`p` does not divide
`r`), retain the literal half-twist mask from THM-3420,

```text
H_(p,r)={ell in Z/pZ: ||r(2ell+1)/(2p)||<1/14}.        (1)
```

The theorem is

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
For the even prime `p=2`, every literal half mask is empty; the companion
checks this boundary separately.

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
pairs.  If `B` is the sum, over nonfixed reflection orbits, of multiplicity
minus one, then

```text
Omega=(e-1)+2B.
```

Direct substitution also gives

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

## 4. THM-3423 gives the exact ratio complement and clique bound

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

THM-3423, specialized at `h=7`, gives the exact identity

```text
F_p^* minus (O_L/O_L) = D                             (13)
```

and the exact packing bound

```text
maximum number of pairwise-disjoint dilates of O_L =3. (14)
```

Indeed, the theorem's ratio threshold is

```text
p>=8(56+c),
```

whose worst case in `(11)` is `p>=512`, and its modular clique threshold is
`p>2*6^3=432`.  Its proof colours the rational Cayley graph by `nu_2 mod 3`;
the sharp clique is `{1,2,4}`.  The present companion independently replays
the `h=7` graph (`25` vertices, `84` edges, `60` triangles), its exact modular
stability bound `217`, and the first large prime in each residue class.

## 5. Spending overlap produces too many disjoint odd blocks

At any reflection orbit, if `t` odd blocks meet there, deleting `t-1` of
them removes every odd/odd collision on that orbit.  Doing this over all
orbits deletes at most the total off-fixed pair budget `B`; overlaps involving
even blocks only consume more of that same budget.  Therefore the
`7-e` odd blocks contain a pairwise-disjoint subfamily of size at least

```text
7-e-B.                                                (15)
```

For `s=9`, `(8)` makes this lower bound `4` for every `e=1,2,3`.  For
`s=11` it is `5`, and for `s=13` it is `6`.  Each contradicts the clique
bound three.

For `s=1`, the budget is zero but `(15)` alone would weaken when `e` grows.
Here any two even blocks already collide off the fixed sheet.  Indeed, for
`k>=13` one has `p<(k+1)^2`; the circular-gap argument on

```text
0,lambda,...,k lambda
```

gives `d lambda=+-a` with `1<=a,d<=k`, so `E/E=F_p^*`.  Because `B=0`, this
forces `e=1`.  The remaining six odd blocks are pairwise disjoint, again
contradicting clique number three.

Equations `(13)--(15)` exclude every prime `p>=512` outside the already
positive all-even `p=29` case.

## 6. Exact finite boundary and positive atoms

The companion performs an independent normalized-mask census for every one
of the `96` primes below `512`.  For each scalar-feasible value of `e`, common
multiplicative scaling fixes one even block, and a bit-mask DFS chooses the
remaining even and odd dilates while spending the exact overlap budget.  The
explicit universe is

```text
3<=p<512, p prime;
1<=e<=7;
197 scalar-feasible (p,e) profiles;
seven coefficient representatives, distinct modulo sign within each
parity bank; accidental coincident masks retained.    (16)
```

The coefficient-representative restriction loses no negative case.  If two
same-parity blocks repeat a coefficient modulo sign, their masks coincide;
deleting one gives a cover by at most six blocks, and THM-3416 already
classifies those shorter prime covers as `p=11,23`.  The census does not
identify accidental coincidences between inequivalent coefficient
representatives: it retains them.  In particular, its positive `p=11,e=2`
profile uses two inequivalent even coefficients whose masks are the same
reflection-fixed singleton.

The positive support is exactly

```text
{11,13,23,29}.                                        (17)
```

There are twelve positive seven-block parity profiles: one at `11`, one at
`13`, three at `23`, and all seven values of `e` at `29`.  The `e=7` row at
`29` is the Paley fixed-zero splitter transported to the half layer.  Normal
and optimized runs are byte-identical.  The four literal witnesses in `(4)`
are replayed separately, including their full multiplicity profiles.

Together with THM-3416's all-prime rank-at-most-six classification, `(16)`
also covers a hypothetical shorter finite-boundary cover: away from `11,23`
it would have to be an essential seven-block cover and hence appears in the
normalized census.  This proves the finite half of `(2)` and completes the
proof.

## 7. Equality boundary, losses, and non-consequences

- The analytic threshold is explicit: the worst odd-interval offset is
  `c=8`, giving `p>=512`; every smaller odd prime is in `(16)`, and `p=2`
  is checked separately.
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
