---
id: THM-4214
title: "Two-newcomer Pascal-complete eleven-body Haar charts"
status: >
  PROVED ORGANIZATIONAL SYNTHESIS RELATIVE TO THM-4150/4156/4191/4203/4207/
  4211 +
  VERIFIED-EXACT + INDEPENDENTLY AUDITED ZERO-LAYER CONTROL. The safe faces
  through rank eleven contain the join of the displayed eighteen-label
  simplex with the fixed-50 newcomer star. For every outsider r, all
  C(20,11)=167,960 eleven-subsets of C union {50,r} have Haar-safe mass at
  least 4/63, partitioned by newcomer count as
  31,824+87,516+48,620. More generally the finite m-leaf star has
  f_j=C(18,j)+(m+1)C(18,j-1)+mC(18,j-2) certified rank-j faces for
  0<=j<=11. For Q={50,51}, {49,50}, or {50,r} with r>=5682, all
  C(32,11)=129,024,480 eleven-subsets of P union Q are also safe. A labelled
  Pascal-block/combinadic rank bijects each displayed complete rank-eleven
  chart with an initial interval of the natural numbers. This theorem proves
  no arbitrary-pair entry or LRC(14).
source: codex-lrc-pascal-complete-20260826
depends_on:
  - THM-4150-safe-set-haar-measure-universal-odd-tail-lrc14-transfer
  - THM-4156-divisor-complete-anchor-pool-haar-odd-tail-transfer
  - THM-4191-complete-full-pool-newcomer-haar-transfer
  - THM-4203-fixed-pool-seventeen-body-depth-eight-haar-completion
  - THM-4207-two-newcomer-sharp-depth-transition-base-surplus-composition-and-variable-pool-chart-number
  - THM-4211-fixed-fifty-cofinal-two-newcomer-haar-tail-and-eighteen-label-chart
script: 04-computation/lrc14_two_newcomer_pascal_complete_chart_thm4214.cpp
output: 05-knowledge/results/lrc14_two_newcomer_pascal_complete_chart_thm4214.out
independent_audit_script: 04-computation/lrc14_two_newcomer_pascal_complete_chart_independent_audit_thm4214.cpp
independent_audit_output: 05-knowledge/results/lrc14_two_newcomer_pascal_complete_chart_independent_audit_thm4214.out
script_sha256: 2f2a4107adaf89ddef15c7b6b76fca8b8093980cc64ae9f2d466f49fbadfa9c6
output_sha256: 6ac42b701fb06adcafdc3044e08ef8f56eb4577bc3e0200c383051cc74268e37
independent_audit_script_sha256: 3597158a0ce7acd991e728fb68441838e2a0dc8ceda73162b637e81672d9e9cd
independent_audit_output_sha256: edc7f977894bdfbabc1c052e5a3eaaa217d09c9427c5957d1cbe8a2c79b809e3
hash_basis: raw LF bytes
primary_audit: >
  PASS. Failure-mask atom integration plus a subset zeta transform directly
  checks all 31,824 zero-newcomer eleven-bodies, the exact Pascal partition,
  the mod-six control split, and the bare-Bonferroni hostile. This is a
  redundant verification of the THM-4156 layer, not a new dependency.
independent_audit: >
  ACCEPT. A separate per-body labelled safe-comb event sweep uses no Boolean
  transform or shared mass table. It agrees on every zero-layer mass, the
  strict minimizer, the semantic ledger, and the rank-nine Bonferroni hostile.
---

# THM-4214 -- two-newcomer Pascal-complete eleven-body Haar charts

**PROVED ORGANIZATIONAL SYNTHESIS RELATIVE TO THM-4150/4156/4191/4203/4207/
4211 +
VERIFIED-EXACT + INDEPENDENTLY AUDITED ZERO-LAYER CONTROL; LRC(14) REMAINS
OPEN.**

## 1. Statement and inheritance pass

Retain the thirty-label pool

```text
P={8,10,15,16,20,30,40,42,60,63,80,84,85,88,95,
   120,126,132,143,145,168,170,176,190,193,240,252,
   264,286,290},                                        (1)
```

its eighteen-label subchart

```text
C={8,16,40,42,80,84,85,88,95,
   120,126,143,145,168,193,240,252,286},                (2)
```

and, for a finite positive set `S`,

```text
G_S={y in R/Z:min_(s in S)||sy||>=1/14},
G_empty=R/Z,                    alpha=4/63.              (3)
```

Let `R` be any finite set of positive labels disjoint from
`P union {50}`, and put `m=|R|`. On newcomer vertices `{50} union R`, let
`Sigma_R` be the simplicial star whose only two-vertex faces are

```text
{50,r},                         r in R.                 (4)
```

Thus `Sigma_R` has one empty face, `m+1` vertices, and `m` edges. Define

```text
K_R={L union T:L subset C, T a face of Sigma_R,
                 |L|+|T|<=11}.                         (5)
```

> **Pascal-star theorem.** Every `F in K_R` satisfies
>
> ```text
> mu(G_F)>=alpha.                                       (6)
> ```
>
> With the convention `binom(n,k)=0` for `k<0` or `k>n`, the exact number
> of certified cardinality-`j` faces is
>
> ```text
> f_j=binom(18,j)+(m+1)binom(18,j-1)+m binom(18,j-2),
>                                             0<=j<=11. (7)
> ```

For one leaf `R={r}`, where necessarily `r notin P union {50}`, the star is
the complete simplex on `{50,r}`. Hence `(5)` is the entire rank-at-most-
eleven downset of the twenty-label set

```text
A_r=C union {50,r}.                                    (8)
```

In particular, every `H in binom(A_r,11)` satisfies `(6)`. Its exact
newcomer-count partition is

```text
zero newcomers:  binom(18,11)   = 31,824,
one newcomer:  2 binom(18,10)   = 87,516,
two newcomers:   binom(18,9)    = 48,620,
                                     -------
                                     167,960=binom(20,11). (9)
```

There is also a full-pool completion on the proved pair families. If

```text
Q={50,51},  Q={49,50},  or  Q={50,r} with r>=5682,      (9a)
```

then every `H in binom(P union Q,11)` satisfies `(6)`. Its four labelled
outsider-incidence blocks have sizes

```text
binom(30,11), binom(30,10), binom(30,10), binom(30,9),
```

and hence contain exactly

```text
54,627,300+30,045,015+30,045,015+14,307,150
  =129,024,480=binom(32,11)                              (9b)
```

eleven-bodies for each pair in `(9a)`.

For each rank-eleven `H` in `(5)`, or each
`H in binom(P union Q,11)` with `Q` from `(9a)`, every positive integer `c`
and every two distinct positive odd integers `a,b`, THM-4150 supplies some
`x in R/Z` with

```text
min_(v in 2cH union {a,b})||vx||>=1/14.                (10)
```

The closest proved mechanisms are THM-4156/4203's small/full zero-newcomer
layers, THM-4191's complete one-newcomer layer, and THM-4207/4211's fixed-pair
two-newcomer charts. The canonical hostile is the exact rank-nine body in
Section 5, where the bare two-danger-set Bonferroni bound is negative. The
corrected near miss is to call the inherited rows separate charts and forget
that their newcomer incidences glue. The least-used sidecar is precisely the
newcomer count together with the central-versus-leaf incidence retained by
`Sigma_R`.

The live concept board was

```text
common-pool simplex | one-newcomer cone | fixed-50 star
full-pool Pascal blocks | combinadic rank | mod-six hostile.             (11)
```

The selected method cards are “search the statement before the method” and
“controlled forgetting requires a sidecar”: the theorem is obtained by
recovering the current body theorems, then retaining which newcomer vertices
each one actually permits.

## 2. Proof by newcomer count

Take `F=L union T` as in `(5)`.

If `|T|=0`, then `F subset C subset P`. THM-4156 proves the stronger common
pool statement

```text
mu(G_P)>alpha.                                         (12)
```

Safe-set inclusion reverses label inclusion, so `G_P subset G_F`, proving
`(6)` in the zero-newcomer layer.

If `|T|=1`, write `T={q}`. Since `|L|<=10`, extend `L` inside `P` to some
`B in binom(P,10)`. Every newcomer vertex of `Sigma_R` lies outside `P`,
including `50`, so THM-4191 applies and gives

```text
mu(G_(B union {q}))>=alpha.                            (13)
```

Because `F subset B union {q}`, `(13)` and monotonicity prove `(6)`.

If `|T|=2`, the star definition forces

```text
T={50,r}                         for some r in R.       (14)
```

Now `|L|<=9`. Extend `L` inside `C` to some `K in binom(C,9)`.
THM-4211 applies because `r notin P union {50}` and gives

```text
mu(G_(K union {50,r}))>=alpha.                         (15)
```

Again `F subset K union {50,r}`, so `(15)` proves `(6)`. These are all faces
of the newcomer star, completing the proof.

This also explains downward closure rather than merely asserting it. Each
lower face is extended *within the same permitted newcomer face*: to ten
pool labels in the one-newcomer row and to nine chart labels in the central
two-newcomer row. No forbidden leaf-leaf edge is introduced.

For the full-pool assertion `(9a)`, partition an eleven-set `H` by
`|H intersect Q|`. With zero outsiders, THM-4203 applies because `11<=17`.
With one labelled outsider, THM-4191 applies to the unique ten-subset of `P`.
With both outsiders, THM-4207 applies to `{50,51}`, while THM-4211 applies
to `{49,50}` and to `{50,r}` for `r>=5682`. These four cases are disjoint and
exhaustive. This proof does not infer a new pair by interpolation in its
labels.

## 3. Exact face enumeration and the Pascal completion

The cardinality enumerator of the full simplex on `C` is `(1+z)^18`. The
newcomer-star enumerator is

```text
1+(m+1)z+mz^2.                                         (16)
```

Faces of a simplicial join are disjoint unions of one face from each factor,
so the coefficient of `z^j` in the product of these two polynomials is
exactly `(7)`. The rank bound `j<=11` is imposed only after taking that
coefficient.

When `m=1`, `(16)` is `(1+z)^2`, and hence

```text
(1+z)^18(1+2z+z^2)=(1+z)^20.                           (17)
```

The coefficient of `z^11` is the three-row decomposition `(9)`. This is the
meaning of **Pascal-complete** here: for each fixed leaf `r`, the zero-, one-,
and two-newcomer theorems fill every cell in the eleventh Pascal row of the
same twenty-label chart, and their monotone extensions fill its complete
rank-at-most-eleven downset.

For `m>1`, `(16)` is not `(1+z)^(m+1)`. The missing terms are the exact
leaf-leaf and higher-newcomer debt, not a counting accident.

## 4. Independent exact zero-layer control

The proof already owns the zero-newcomer layer through the stronger
THM-4156 statement `(12)`. Two new exact paths nevertheless audit this layer
directly, so a transcription error in `(2)` or the Pascal gluing cannot hide
behind the inherited citation.

The primary path builds the complete wall arrangement of `C`, accumulates
failure-mask atoms, and applies a subset zeta transform. It obtains

```text
D=18,241,159,416,480,        cells=4,371,
binom(18,11)=31,824 bodies,   failures=equalities=0.    (18)
```

The strict minimum is

```text
H_min={16,40,42,85,88,95,145,193,240,252,286},
D mu(G_(H_min))=3,306,338,342,934,
63D mu(G_(H_min))-4D=135,334,677,938,922>0.             (19)
```

The independent path enumerates every eleven-body separately and intersects
its labelled safe-comb intervals by an event sweep. It has no Boolean-lattice
transform or shared mass table and agrees on every mass, `(18)--(19)`, and
the semantic ledger

```text
FNV1A64_LE=f594af0b446a5bf2.                           (20)
```

This is redundant verification, not a claim that the zero layer was open
before this theorem.

## 5. Why the center label is load-bearing

For any nine-body `K` and arbitrary new labels `q,r`, the union bound gives
only

```text
mu(G_(K union {q,r}))>=mu(G_K)-2/7.                    (21)
```

To reach `alpha`, this scalar argument would require

```text
mu(G_K)>=22/63.                                        (22)
```

Both exact paths find the hostile

```text
K_0={8,80,85,88,95,143,145,168,193},
D mu(G_(K_0))=4,286,977,525,566,
63D mu(G_(K_0))-22D=-131,225,923,051,902<0.            (23)
```

Thus the first attempted implication

```text
large zero/one-newcomer downsets
  => arbitrary two-newcomer completion by bare Bonferroni               (24)
```

already fails on the displayed chart. Equation `(23)` does **not** say that
`K_0 union {q,r}` is unsafe for any pair. It says only that the generic
measure lower bound has discarded the overlap coordinate. THM-4207's exact
base-surplus lemma diagnoses the same loss abstractly. THM-4211 restores the
needed literal geometry only on central edges `{50,r}`; this theorem has no
lawful reason to add an edge `{r,s}` between two leaves.

## 6. Divisor-clock audit

MISTAKE-511 requires every advertised physical family to be tested against
the inherited mod-six clock. The chart `C` has seven multiples of `6` and
eleven nonmultiples, while `50` is a nonmultiple. Therefore the twenty-label
chart `A_r` has exactly twelve nonmultiples if `6|r`, and thirteen otherwise.
Among its `167,960` rank-eleven bodies, the exact number containing no
multiple of `6` is consequently

```text
binom(12,11)=12                       if 6|r,
binom(13,11)=78                       if 6 does not divide r.             (25)
```

For these scale-one bodies, the physical row `2H union {a,b}` is already
closed at `x=1/12` for every distinct odd pair: the body clearance is at least
`1/6` and the tail clearance at least `1/12`. The remaining `167,948` or
`167,882` bodies contain a multiple of `6` and therefore defeat that one
clock, but no claim is made that they evade every earlier certificate. The
Haar conclusions `(6)` and `(10)` remain valid in all cases and for every
allowed common scale.

## 7. Full-pool Pascal blocks and ordinal rank

Fix an order `Q=(q_0,q_1)` on any pair in `(9a)`, order `P` increasingly,
and let `rho_(n,k)` be a fixed zero-based combinadic rank on the `k`-subsets
of an ordered `n`-set. Define

```text
rank_Q(H)=
  rho_(30,11)(H),                                      H intersect Q=empty;
  binom(30,11)+rho_(30,10)(H intersect P),             H intersect Q={q_0};
  binom(30,11)+binom(30,10)+rho_(30,10)(H intersect P),
                                                       H intersect Q={q_1};
  binom(30,11)+2binom(30,10)+rho_(30,9)(H intersect P),
                                                       H intersect Q=Q. (26)
```

The block offsets are the partial sums in `(9b)`, so `(26)` is a bijection

```text
binom(P union Q,11) <-> {0,1,...,129,024,479}.          (27)
```

Replacing `P,30` by `C,18` gives the analogous bijection from the complete
one-leaf chart in `(9)` to `{0,1,...,167,959}`. The block recovers the exact
outsider-incidence word and the combinadic coordinate recovers the pool
subset. Thus the rank forgets the proof certificate but not the discrete
body. Coincidences after common scaling or between different pair families
are not asserted to be disjoint. Adding one to either zero-based rank gives
the literal one-based ordinal `1,...,N`.

Before selecting a second outsider, THM-4203 and THM-4191 similarly fill

```text
binom(30,11)+binom(30,10)=84,672,315=binom(31,11)       (28)
```

eleven-bodies in every one-newcomer pool `P union {q}` with `q notin P`.
The connection contract for the full-pool completion is

```text
source:       H in binom(P union Q,11)
target:       the zero-, labelled-one-, or two-outsider body theorem
map:          H -> (H intersect Q,H intersect P)
preserved:    every label, outsider identity, cardinality, Haar threshold
destroyed:    repair edge, wall address, and chosen safe phase
sidecar:      the two-bit outsider-incidence word
decisive test: all four labelled Pascal blocks are covered.              (29)
```

## 8. Exact artifacts and replay

Primary:

```text
04-computation/lrc14_two_newcomer_pascal_complete_chart_thm4214.cpp
sha256 2f2a4107adaf89ddef15c7b6b76fca8b8093980cc64ae9f2d466f49fbadfa9c6

05-knowledge/results/lrc14_two_newcomer_pascal_complete_chart_thm4214.out
sha256 6ac42b701fb06adcafdc3044e08ef8f56eb4577bc3e0200c383051cc74268e37
```

Independent:

```text
04-computation/lrc14_two_newcomer_pascal_complete_chart_independent_audit_thm4214.cpp
sha256 3597158a0ce7acd991e728fb68441838e2a0dc8ceda73162b637e81672d9e9cd

05-knowledge/results/lrc14_two_newcomer_pascal_complete_chart_independent_audit_thm4214.out
sha256 edc7f977894bdfbabc1c052e5a3eaaa217d09c9427c5957d1cbe8a2c79b809e3
```

Replay from the repository root:

```bash
mkdir -p /tmp/thm4214-replay

g++ -std=c++20 -O2 -Wall -Wextra -Wpedantic -Werror \
  04-computation/lrc14_two_newcomer_pascal_complete_chart_thm4214.cpp \
  -o /tmp/thm4214-replay/primary
/tmp/thm4214-replay/primary | diff -u \
  05-knowledge/results/lrc14_two_newcomer_pascal_complete_chart_thm4214.out -

g++ -std=c++20 -O2 -Wall -Wextra -Wpedantic -Werror \
  04-computation/lrc14_two_newcomer_pascal_complete_chart_independent_audit_thm4214.cpp \
  -o /tmp/thm4214-replay/independent
/tmp/thm4214-replay/independent | diff -u \
  05-knowledge/results/lrc14_two_newcomer_pascal_complete_chart_independent_audit_thm4214.out -
```

Both programs also byte-match at `-O0` and `-O3`. Clang undefined-behavior
sanitizer runs emit the same streams with no diagnostics.

## 9. Strict scope

This theorem proves two compatible organizations. On `C`, THM-4156 owns the
zero-newcomer simplex, THM-4191 the one-newcomer cone, and THM-4211 the
fixed-50 star edges. On `P` and the pairs `(9a)`, THM-4203 owns the zero block,
THM-4191 the labelled one blocks, and THM-4207/4211 the two block. The new
content is that these inputs glue without losing newcomer incidence and admit
the exact ordinal ranks `(26)--(28)`.

It does **not** prove a leaf-leaf edge `{r,s}` in the arbitrary star, any
three-newcomer face, replacement of the center `50`, arbitrary two-newcomer
safety on `C` or `P`, entry of an arbitrary LRC(14) instance into the charts,
or full LRC(14).
