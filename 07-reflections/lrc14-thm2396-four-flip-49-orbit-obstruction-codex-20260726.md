# LRC14 THM-2396 four-flip `49`-orbit obstruction

**Status:** FINITE-EXACT INDEPENDENT REFEREE OF PROVED THM-2396.

This is an independent, strictly relaxed finite referee for proved and
twice independently hostile-audited THM-2396,
`THM-2396-common-core-forty-nine-orbit-word-incompatibility`. It is a
reproduction/sharpening sidecar, not a separate theorem dependency.

THM-2393 and THM-2394 reduce the no-clean branch to

```text
(C_1,C_2,c_1,c_2)=h(1,13,13,169)
```

and make each of the six top-safe septimal bins an exact six-address
transversal. THM-2391 makes the top-danger bin equally rigid: its
two-address guard is the disjoint union of the `h` and `13h` core
addresses, and every lower ordinary word lands at one of those two
addresses.

The new finite result is a sharp incompatibility between those seven
bins.

## 1. A deliberately relaxed universe

Work on `Z/49Z` and normalize the top-danger bin to

```text
Q=7Z/49Z.
```

An ordinary unit danger is any translate of a seven-term arithmetic
progression with unit step. A guard is any translate of a fourteen-term
unit-step progression. There are exactly

```text
1029 ordinary masks,
1029 guard masks.
```

Multiplication of the orbit coordinate by the unit `h` preserves `Q`
and normalizes the three common-core steps to

```text
A=D_h:       1,
B=D_(13h):   13^(-1)=34 mod 49,
C=D_(169h):  169^(-1)=29 mod 49.
```

The search allows the translates of `A`, `B`, and `C` to vary
**independently**. It also allows the guard and all four lower ordinary
words independent steps and translates. This strictly enlarges the
physical common-phase universe.

Only two canonical constraints are retained.

1. On `Q`, the guard is `A disjoint_union B`, and the four lower
   ordinary words lie at those two addresses.
2. On each of the six safe bins, guard plus the four lower words is an
   exact transversal of six of the seven addresses. Its hole is `B`,
   unless it is the exceptional common address `A=C`.

## 2. The four-flip lemma

Fix a guard. Among ordinary masks which use a guard address on `Q` and
avoid the guard on every safe bin, the candidate-count distribution is

```text
22 candidates: 294 guards,
25 candidates: 294 guards,
26 candidates: 147 guards,
33 candidates: 294 guards.
```

Enumerating four-word exact covers gives only `1029` distinct
guard/hole words in total:

```text
0 hole words: 588 guards,
1 hole word:  147 guards,
3 hole words: 294 guards.
```

The `1029` guards form only four affine orbits under the transformations
which preserve `Q`: two `294`-element orbits realize no hole word, one
`294`-element orbit realizes three, and one `147`-element orbit realizes
one. Thus the large-looking enumeration has a four-type symmetry
certificate.

Now range over independent translates of normalized `A` and `B` and
retain the top-bin partition. There are exactly

```text
100842
```

top-compatible `(A,B,guard,hole-word)` cases to test. In every case
where each non-`B` hole is at the corresponding `A` address, at least
four of the six safe bins have a non-`B` hole:

```text
minimum number of A-for-B flips =4,

survivors with at most three flips =0.               (1)
```

Thus the finite obstruction has one bin of slack beyond what the
physical geometry needs.

## 3. Why the common core cannot pay four flips

The exact independent-translate `A/C` intersection census is

```text
0 equal addresses:  490 pairs,
1 equal address:   1421 pairs,
2 equal addresses: 490 pairs.
```

In particular, `A` and `C` agree in at most two of all seven bins.
THM-2394 permits a non-`B` hole only at an `A=C` bin. The physical
common core can therefore pay at most two flips, while (1) requires at
least four. Hence the relaxed `49`-orbit universe is empty.

This proof does not need the common phase linking the three core
translates. It does not need THM-2395's successor-shell estimate or its
carry automaton. The obstruction is static but uses all seven bins
together.

## 4. Sharp boundary

Four flips are genuinely possible once the `A=C` budget is relaxed.
One exact positive control is

```text
A = {3,4,5,6,7,8,9},

B = {0,11,15,19,30,34,45},

guard
  ={0,4,7,10,17,20,23,27,30,33,36,40,43,46},

safe-bin holes
  =(1,1,6,1,0,0),

A-for-B flip bins
  ={1,2,5,6}.
```

The four lower masks are

```text
{1,7,13,19,25,38,44},
{0,15,16,31,32,47,48},
{3,7,18,22,26,37,41},
{2,7,12,24,29,34,39}.
```

They satisfy the top partition and produce the six safe-bin
transversals pointwise. This proves that the number four in (1) is a
real equality boundary, not search slack.

## 5. Scope

Canonical THM-2396 localizes the same static contradiction orbitwise.
Every base in the high-safe set of mass `66/91` has at least one clean
root among its `49` roots. Exact disintegration therefore gives

```text
common-core clean mass >=66/(91*49)=66/4459.
```

Together with THM-2392 this yields

```text
q_*-labelled charged cell >=33/115934,
complete blocker cell      >=11/57967,
owner-resolved tensor      >=33/753571.
```

THM-2393's complementary packets have floor `1/26754`, so the uniform
last-lane clean floor is now `1/26754`.

This does not by itself turn a clean-hole coefficient into the
canonical terminal target current:
THM-2392/2391 still retain a same-parent charged
`C_7 x C_13` tensor, while endpoint/owner transport remains a separate
interface. No ledger decrement is claimed in this reflection.

## Exact reproduction

Run

```bash
python3 04-computation/lrc14_thm2396_relaxed_49_orbit_independent_referee.py
python3 -O 04-computation/lrc14_thm2396_relaxed_49_orbit_independent_referee.py
```

and compare both transcripts after LF normalization with

```text
05-knowledge/results/lrc14_thm2396_relaxed_49_orbit_independent_referee.out
```

Both modes currently match the stored `16`-line transcript. LF-normalized
SHA256:

```text
script 1a397efff2315750296b7a494858bf070b9b4d045bbe87c79e525ce8069f925a
output 7a4d909d2fea0418babeaeeebdf054e3487914ff1d8f0e308564c010b2bd524f
```

An independent hostile audit accepted the physical embedding, strict
relaxation, mask completeness, exact-cover recursion, all finite
ledgers, sharp hostile, and quantitative localization. No defect
remains.
