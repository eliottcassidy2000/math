---
id: THM-814
title: At n=8 the black B2/B3 orbit codec has sixteen fixed-layer positional collisions, while every source-q stratum remains outward-dominant
status: PROVED (all-size black endpoint-q denominator polynomial) + FINITE-EXACT (atlas-free n=8 codec refuter and source-q flow)
source: codex-2026-07-15-S13
depends_on: [THM-643, THM-785, THM-790, THM-809, THM-811]
related: [THM-801, HYP-6880]
verification:
  - 04-computation/black_curvature_n8_boundary_codex_S13.cpp
  - 05-knowledge/results/black_curvature_n8_boundary_codex_S13.out
---

# THM-814 — the `n=8` black-orbit boundary

Two conjectural continuations of THM-811 separate at `n=8`.

1. Projected endpoint-node pair plus the reflection orbit of `(B2,B3)` is no
   longer a complete black-edge-orbit codec.  It has exactly sixteen double
   collisions among 523,264 black reflection orbits.
2. Source-normalized mixed-to-pure-black flow remains strictly larger than
   the reverse flow in every nonempty source-curvature stratum `q=0,...,6`.

Both decisions are atlas-free.  The first uses tournament isomorphism only
inside exact `(B2,B3)` candidate cells.  The second uses THM-643's direct
self-converse characterization of black mixed endpoints.

## 1. Exact orbit key without an `n=8` atlas

Let `t` be the apex-zero endpoint of a black complement line and `sigma` its
staircase reflection.  Rank `B2(t)` by the weak-composition ranks of the four
binary states `00,01,10,11` on layers `tau=3,...,7`, followed by the number of
ones in the three fixed `tau=8` positions.  The mixed radices are

```text
(4,4,10,10,20,4),             0 <= rank(B2) < 128000.   (1)
```

Rank `B3(t)` by the one-counts in strata

```text
(A,B,C,AB,AC,BC,ABC)
```

of sizes `(1,1,1,4,4,4,6)`, using radices

```text
(2,2,2,5,5,5,7),              0 <= rank(B3) < 7000.     (2)
```

Put

```text
lit(t)=7000 rank(B2(t))+rank(B3(t)),
K(e)=sort(lit(t),lit(sigma t)).                          (3)
```

This is an exact key below 60 bits, not a hash.  Streaming all black lines and
retaining one representative per free reflection orbit gives

```text
black lines                         1,046,528
black reflection orbits               523,264
K cells                               331,500
K collision cells                     144,500
K collision excess                    191,764
within-K candidate pairs              275,368
maximum K multiplicity                     10.           (4)
```

For each candidate pair, build only its four eight-vertex endpoint
tournaments.  Degree-bucket backtracking decides ordinary tournament
isomorphism, allowing the two line endpoints to swap.  Because (3) already
contains the whole staircase-reflection orbit, this exactly compares unordered
converse-merged endpoint-node pairs.  The refinement has

```text
cells                                 523,248
collision cells                            16
collision excess                           16.           (5)
```

No class atlas or canonical node identifier enters this proof.

## 2. The sixteen collisions are one positional kernel

The complete collision list is

```text
1725 9789       1981 10045       5821 13885       6077 14141
296613 304677   296869 304933    300709 308773    300965 309029
593565 601629   593821 601885    597661 605725    597917 605981
741086 749150   741342 749406    745182 753246    745438 753502. (6)
```

Every pair satisfies

```text
line_0 xor line_1 = 0x02080,                             (7)
```

which swaps exactly tile bits 7 and 13, namely `(7,2)` and `(6,3)`.  Both are
staircase-fixed `tau=8` positions in the common `ABC` stratum.  Consequently
raw `B2` sees only the unchanged fixed-layer one-count and `B3` sees only the
unchanged `ABC` count.  Their join has deliberately forgotten position.

The first positional moment of the fixed layer separates all sixteen pairs.
So does THM-809's lower-face `Lambda` address.  This proves that the failure is
the predicted count-without-position kernel, not a canonicalization accident.

All 32 endpoint tournaments in (6) are non-self-converse, so these are
pure-black edges.  None is a projected loop.  Every pair has

```text
(q0,q1,Delta C3)=(0,2,0),                               (8)
```

and the Smith defects are `-4,-2,0,2`, four collision pairs at each value.
Thus neither boundary curvature nor Smith current can repair the codec.  The
missing bit is positional.

## 3. Atlas-free source-category classifier

THM-643 implies, at a black tiling endpoint,

```text
mixed node       <=> endpoint tournament is self-converse,
pure-black node  <=> endpoint tournament is not self-converse. (9)
```

Self-converseness is tested directly by degree-complement-compatible
anti-isomorphism backtracking.  On all `2^21` fixed-path tilings the exact
census is

```text
score-complement-symmetric candidates       744,678
self-converse tilings                         58,712
blue tilings                                   4,096
black mixed endpoints                         54,616.    (10)
```

This independently recovers the 53,072 black cross-category lines: 12,584
are `C3` ties and 40,488 are non-tied.

## 4. The all-size black endpoint denominator

Put

```text
M=binom(n-1,2),       r=n-4,
R=(M+floor((n-1)/2))/2.
```

Let `E_K(x)` count black line endpoints by their own boundary curvature `q`.
Then

```text
E_K(x)=2^(M-2n+5)(4+(1+x)^2)(3+x)^r
       -2^(R-n+2) G_n(x),                               (11)
```

where

```text
G_n(x)=(3+x^2)^((n-2)/2)                    if n is even,
G_n(x)=(1+x)(3+x^2)^((n-3)/2)               if n is odd. (12)
```

Proof: specialize THM-801/811's black endpoint-pair polynomial first as
`K(x,1)` and then as `K(1,x)`, and add.  The first term in (11) is the total
endpoint specialization; (12) is the two blue endpoint specializations.

At `n=8`, (11) gives the exact coefficient vector

```text
(412992,718848,578880,282624,84416,14336,960),          (13)
```

which the direct self-converse/non-self-converse split reproduces.  Reflection
acts freely on black lines and preserves `q`, `C3`, and self-converse status,
so every denominator and every directed count below is even.

## 5. Source-`q` flow survives `n=8`

Orient each non-tied mixed--pure-black line toward increasing `C3`.  The exact
source-disintegrated table is

| `q` | outward / mixed mass | reverse / pure-black mass |
|---:|---:|---:|
| 0 | `5296/8282 = 2648/4141` | `8620/404710 = 862/40471` |
| 1 | `5166/16240 = 369/1160` | `9086/702608 = 4543/351304` |
| 2 | `1326/19160 = 663/9580` | `1356/559720 = 339/139930` |
| 3 | `1820/7844 = 455/1961` | `2496/274780 = 624/68695` |
| 4 | `2158/2734 = 1079/1367` | `2422/81682 = 1211/40841` |
| 5 | `180/180 = 1` | `368/14156 = 92/3539` |
| 6 | `172/176 = 43/44` | `22/784 = 11/392` |

The left rate is strictly larger in all seven rows.  Yet raw counts point in
the opposite direction:

```text
mixed -> pure-black       16,118
pure-black -> mixed       24,370.                        (14)
```

Thus the black flow reversal is quantitatively a phase-volume effect through
the first untested size.  Curvature `q` controls that disintegration even
though (8) proves it is not an edge identity codec.

## 6. Tournament Analysis and preservation boundary

Tournament Analysis takes the seven `q` strata as vertices.  Its pairwise
observable is which stratum has the larger outward-minus-reverse signal; the
two switches are raw count bias and source-normalized rate bias, with
increasing `q` as the tie Hamiltonian path.  Both gauges are transitive but
eight pair orientations flip between them.

The challenged assumption is that sufficient aggregate face counts identify
a line orbit.  They do through `n=7` and fail at `n=8` exactly when the fixed
layer acquires multiple positions invisible to both counts.  The safe object
must retain a fixed-layer word or positional moment, in addition to the
reflection orbit and endpoint nodes.

These carriers preserve tournament edge identity or black-flow information;
they do not preserve owner labels, metric wall position, centered-CF phase,
or the LRC loneliness predicate.  THM-813's `Q` sidecar and THM-808's
owner/root stalk remain differently typed obligations.
