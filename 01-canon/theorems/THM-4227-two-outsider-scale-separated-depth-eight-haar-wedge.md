---
id: THM-4227
title: "Two-outsider scale-separated depth-eight Haar wedge"
status: >
  PROVED RELATIVE TO THM-4150/4170 + FINITE-EXACT + INDEPENDENTLY AUDITED.
  Every nine-body in the displayed thirty-label pool is Haar-safe after
  adjoining any ordered pair of outsiders in an explicit scale-separated
  wedge. The proof combines a complete depth-eight no-nine-cover census with
  two sequential component-discrepancy estimates. This gives infinite
  two-outsider LRC(14) families but not arbitrary finite pair entry or full
  LRC(14).
source: codex-lrc14-two-outsider-wedge-20260826
depends_on:
  - THM-4150-safe-set-haar-measure-universal-odd-tail-lrc14-transfer
  - THM-4170-triple-deletion-matching-eventual-haar-odd-tail-transfer
related:
  - THM-4188-all-newcomer-zero-original-anchor-hierarchy-and-resonance-filtration
  - THM-4203-fixed-pool-seventeen-body-depth-eight-haar-completion
  - THM-4207-two-newcomer-sharp-depth-transition-base-surplus-composition-and-variable-pool-chart-number
  - THM-4211-fixed-fifty-cofinal-two-newcomer-haar-tail-and-eighteen-label-chart
  - THM-4214-two-newcomer-pascal-complete-eleven-body-haar-charts
  - THM-4231-arbitrary-pair-cofinal-depth-six-haar-repair-and-exact-outsider-lift
script: 04-computation/lrc14_two_outsider_scale_separated_prefix_certificate_thm4227.cpp
output: 05-knowledge/results/lrc14_two_outsider_scale_separated_prefix_certificate_thm4227.out
full_layer_scout_script: 04-computation/lrc14_two_newcomer_deletion_staircase_probe_20260826.cpp
full_layer_scout_output: 05-knowledge/results/lrc14_two_newcomer_deletion_staircase_probe_20260826.out
independent_audit_script: 04-computation/lrc14_two_outsider_scale_separated_wedge_independent_audit_thm4227.cpp
independent_audit_output: 05-knowledge/results/lrc14_two_outsider_scale_separated_wedge_independent_audit_thm4227.out
script_sha256: 2ff07bd8ba229af5364426d2e3ab1313584b78a179e4a336261d64c669a1c4a4
output_sha256: 20262f49f297892bccc16a6601ff64fde669483c9c93a8fd189bfed1d47b0e21
full_layer_scout_script_sha256: 53ffe1642457241153a858c3da3c20cf31e8cc79ab6009f9321895c2f20b06c2
full_layer_scout_output_sha256: 850f67053a7bfcbe73a0e117d5a2613836e32555e60cc8b49df66eb18ba36f5c
independent_audit_script_sha256: 7c2c2378ca538cdae25973865b348c8495f687d6ac78ace83b25fd2f63f878eb
independent_audit_output_sha256: 5e67948795614b35b8f117f0a88e2cf2862c1e6ee9ffc37f86dd9c52ebc21555
primary_include_sha256: 25a6978484c7ab122fdc6c8e1593cfa2ad3468f7184a156045ea7e6cb2efc45d
independent_include_sha256: 58817d5f5e1a8cc07384f3ea82a1feb221af37ab0907afde890ab4fbdd949137
hash_basis: raw LF bytes
audit: >
  PASS/ACCEPT. The primary pool-cell path and an independent explicit
  fixed-wall path both enumerate all 5,852,925 depth-eight deletions and all
  14,307,150 nine-bodies. They agree on 5,267,460 strict double-limit edges,
  zero equalities, minimum numerator 944,928 at the same labelled edge, and
  zero nine-covers. A component/slack-ordered 859,307-edge prefix already has
  no nine-cover and carries the sharp certificate envelope used here. The
  independent path has frozen semantic ledgers and a different geometry. Both
  quality-enabled sources replay at O3; the inherited full-layer/549 controls
  also byte-match at O0/O3 and under Clang UBSan.
---

# THM-4227 -- two-outsider scale-separated depth-eight Haar wedge

**PROVED RELATIVE TO THM-4150/4170 + FINITE-EXACT + INDEPENDENTLY AUDITED;
LRC(14) REMAINS OPEN.**

## 1. Statement

Retain the thirty-label pool

```text
P={8,10,15,16,20,30,40,42,60,63,80,84,85,88,95,
   120,126,132,143,145,168,170,176,190,193,240,252,
   264,286,290}.                                        (1)
```

For a finite positive label set `S`, put

```text
G_S={y in R/Z:min_(s in S)||sy||>=1/14},
alpha=4/63.                                              (2)
```

Define the exact integers

```text
Q=3,391,
N=321,902,813,232,
K=10,633,545,731.                                       (3)
```

Say that an ordered positive pair `(q,r)` lies in the directed wedge `W` when

```text
q>=Q,
Kr>=N(q+130).                                           (4)
```

> **Two-outsider wedge theorem.** Let `q,r` be distinct positive integers
> outside `P`. If `(q,r)` or `(r,q)` lies in `W`, then for every
> `B in binom(P,9)`,
>
> ```text
> mu(G_(B union {q,r}))>=4/63.                          (5)
> ```

Consequently, for every positive integer `c` and every two distinct positive
odd integers `a,b`, THM-4150 supplies some `x in R/Z` with

```text
min_(v in 2c(B union {q,r}) union {a,b})||vx||>=1/14.  (6)
```

At the first permitted `q`, an exact unsplit error allocation sharpens the
linear-wedge corner to

```text
q=3,391,
r>=106,033.                                             (7)
```

Equation `(4)` gives the clean cofinal linear tail, whose own corner is
`106,590`; the same assertions hold after swapping `q,r`. Equations `(5)--(7)` give
infinite thirteen-speed LRC(14) families with two genuinely variable body
outsiders. THM-4231 separately covers every distinct `q,r>=17548`; this wedge
can enter with one coordinate as low as `3391` when the other is sufficiently
large. Neither theorem gives arbitrary finite pair entry or full LRC(14).

## 2. Inheritance pass and faithful coordinate

The closest proved mechanism is THM-4170's exact discrepancy estimate for
intersecting a finite union of intervals with one runner comb. The canonical
hostile is THM-4207's comparable pair `(50,51)`: separate marginal repair
decks already fail to compose at depth four, so two labels cannot be treated
as independent merely because both are large. The corrected near miss is the
recurring implication “repair-deck cover implies unsafe body,” refuted by
positive direct margins in THM-4203/4207/4211. The least-used sidecar is the
component count after the **first** outsider is adjoined.

The live concept board is

```text
base depth-eight repair | double-limit slack | directed scale order
component birth under first comb | literal comparable-pair deck.          (8)
```

The directed order is load-bearing. The first comb costs `O(c_R/q)` in mass;
it can create `q` new interval components, so the second costs
`O((c_R+q)/r)`. Scale separation makes both losses small. Comparable pairs need
their native joint-wall geometry instead.

## 3. The exact double-limit repair hypergraph

The common wall denominator of `P` is

```text
D=18,241,159,416,480.                                  (9)
```

For `R in binom(P,8)`, write

```text
mu(G_(P\R))=m_R/D,
Delta_R=81m_R-7D,                                      (10)
E_infty={R:Delta_R>0}.                                 (11)
```

The normalization is not cosmetic:

```text
(36/49)mu(G_(P\R))>4/63  iff  Delta_R>0.               (12)
```

Two independent exact implementations prove

```text
|binom(P,8)|=5,852,925,
|E_infty|=5,267,460,
#{R:Delta_R=0}=0,
min_(R in E_infty) Delta_R=944,928,                    (13)
tau(E_infty)>9.                                        (14)
```

The minimum in `(13)` occurs at

```text
R_min={15,16,20,40,170,190,193,240}.                   (15)
```

For `(14)`, both programs inspect every one of the

```text
binom(30,9)=14,307,150                                 (16)
```

labelled nine-bodies and find a disjoint edge in `E_infty`; there are zero
covers. The independent fixed-wall path uses a different edge order, performs
`413,986,509` incidence checks, and freezes the complete deletion ledger
`ed44e240af3d3a59`. Thus `(13)--(14)` are **FINITE-EXACT**, not a sampled
asymptotic claim.

Equations `(10)--(13)` give the uniform strict slack

```text
epsilon_R=(36/49)mu(G_(P\R))-4/63
         =4Delta_R/(441D)
         >=epsilon=4/8,513,189,685.                    (17)
```

The proof uses a component-aware subdeck. For each `R in E_infty`, let `c_R`
be the number of positive-length circular interval components of `G_(P\R)`,
ignoring isolated boundary points, and define its half-slack activation

```text
Q_R=ceil(162 c_R D/(7 Delta_R)).                        (17a)
```

Encode `R` by its thirty-bit pool mask and order `E_infty` by

```text
(Q_R ascending, Delta_R descending, mask ascending).   (17b)
```

Let `E_*` be the first `859,307` edges. Two independent implementations prove

```text
|E_*|=859,307,                 tau(E_*)>9,
max_(R in E_*) Q_R=3,391,
prefix ledger=232b38cfc255acfb.                        (17c)
```

The prefix component/slack envelope has the dominant edge

```text
R_0={60,84,170,240,252,264,286,290},
Delta_0=16,269,324,968,430,
c_0=130,
Q_(R_0)=3,374.                                         (17d)
```

At `Q=3,391`, and hence for every `q>=Q`, it satisfies

```text
Delta_R>=Delta_0,
Delta_0(q+c_R)<=Delta_R(q+c_0) for every R in E_*,

27D/Delta_0=321,902,813,232/10,633,545,731.            (17e)
```

The exact corner ceiling from `(17e)` is `106,590`. If the two discrepancy
losses are optimized jointly rather than split in half, the corner improves
to `106,033`; the clean cofinal wedge `(4)` uses the half-split envelope.

The ordering boundary is also exact. The body

```text
B_*={88,95,170,193,240,252,264,286,290}                (17f)
```

covers every edge with `Q_R<=3,390`; its first missed ordered edge has
`Q_R=3,391`. Thus `3,391` is sharp for this explicitly ordered sufficient
filtration, not a literal outsider threshold or an unsafe-body witness.

## 4. Sequential component-discrepancy proof

Fix `B in binom(P,9)`. By `(17c)`, choose `R in E_*` disjoint from `B` and
put

```text
U=G_(P\R),   M=mu(U).                                   (18)
```

THM-4170's discrepancy lemma gives

```text
mu(U intersect G_q)>=(6/7)M-6c_R/(49q).                (20)
```

The comb `G_q` has `q` interval components. Intersecting unions of `c_R` and
`q` circular intervals has at most `c_R+q` positive-length components: the
boundary of the intersection is contained in the union of the two boundary
sets. Applying the same lemma a second time, after discarding null isolated
points, yields

```text
mu(U intersect G_q intersect G_r)
 >=(36/49)M-36c_R/(343q)-6(c_R+q)/(49r).               (21)
```

For this edge, `(17)` gives

```text
epsilon_R=4Delta_R/(441D).                              (22)
```

By `(17c)` and `q>=Q`, the first error in `(21)` consumes at most half of
this slack:

```text
36c_R/(343q)<=epsilon_R/2.                              (23)
```

The second inequality in `(17e)`, verified at `q=Q`, remains true with `q`
replacing `Q`, because its difference grows by
`(Delta_R-Delta_0)(q-Q)>=0`. Hence `(4)` and `(17e)` give

```text
r >=(321,902,813,232/10,633,545,731)(q+130)
   =(27D/Delta_0)(q+c_0)
  >=(27D/Delta_R)(q+c_R).                              (23a)
```

Rearranging the last inequality shows that the second error in `(21)` is at
most `epsilon_R/2`. Combining `(21)--(23a)` therefore gives

```text
mu(G_((P union {q,r})\R))
 =mu(U intersect G_q intersect G_r)>=alpha.             (24)
```

At the single boundary `q=Q`, keeping the unused first-stage slack gives the
exact edgewise condition

```text
r >= ceil(189Dq(q+c_R)/
          (2(7Delta_R q-81c_R D))).                    (24a)
```

The complete `E_*` census has maximum right side `106,033`, attained at
`R_0`; this proves the sharper corner `(7)`.

Because `R intersect B=empty`,

```text
B union {q,r} subset (P union {q,r})\R.                (25)
```

Safe-set inclusion reverses label inclusion, so `(24)--(25)` prove `(5)`.
The multiplication map `y |-> cy` preserves Haar measure, so
`mu(G_(cH))=mu(G_H)` for every positive integer `c`; THM-4150 then proves
`(6)`. The eleven distinct even body labels cannot collide with the two
distinct odd tails, so `(6)` has exactly thirteen nonzero relative speeds.
Equality in `(4)` is admissible throughout because the target threshold is
the non-strict inequality `mu>=alpha`. **QED.**

## 5. Connection contract and failure boundary

The new bridge is

```text
source:       robust prefix repair R in E_*
target:       a lawful joint repair after outsiders q,r
map:          R -> (M_R,c_R) -> sequential discrepancy bound (21)
preserved:    deletion labels, exact base mass, component envelope, threshold
destroyed:    component addresses, gcd/phase correlation of q and r
sidecar:      the directed outsider order; literal joint walls when comparable
decisive test: tau(E_*)>9 together with inequalities (17c)--(17e).       (26)
```

The boundary also explains the theorem's scope. If `q/r` stays near one,
the second error term in `(21)` does not tend to zero: the first outsider has
created order-`q` boundary complexity. This is a failure of the sequential
envelope, not evidence of an unsafe pair. THM-4207 and THM-4211 prove several
comparable-pair charts by retaining literal phase geometry.

Neither `Q` nor `(7)` is asserted to be a minimal literal entry threshold.
They are sharp only within this explicitly ordered prefix certificate: `Q` is
the boundary of its half-slack activation filtration, while `(7)` is the
optimized joint-loss corner on the resulting prefix. Exact joint geometry or
a different subdeck can improve them. The theorem says nothing about pairs
outside its wedge, a pair with one label inside `P`, physical entry of a
general LRC(14) row into this chart, or LRC(14). THM-4231 independently covers
the cofinal northeast quadrant, but finite pair entry remains open.

## 6. Exact audit and replay

The primary program inherits THM-4188's pool-cell implementation, computes
every edge's exact component/slack activation, and uses Gosper masks. The
independent program instead inherits the explicit fixed-wall referee, rebuilds
the low-arity mass and adjacency atoms, recursively enumerates both subset
universes, and reproduces the quality order and ledger.

Replay from the repository root:

```bash
g++ -std=c++20 -O3 -DNDEBUG -Wall -Wextra -Wpedantic -Werror \
  04-computation/lrc14_two_outsider_scale_separated_prefix_certificate_thm4227.cpp \
  -o /tmp/thm4227-primary
/tmp/thm4227-primary

g++ -std=c++20 -O3 -DNDEBUG -Wall -Wextra -Wpedantic -Werror \
  04-computation/lrc14_two_outsider_scale_separated_wedge_independent_audit_thm4227.cpp \
  -o /tmp/thm4227-independent
/tmp/thm4227-independent
```

Compare stdout with the two outputs named in the frontmatter. Both final
quality-enabled sources were replayed at `-O3`. Their inherited full-layer and
549-prefix controls also byte-match at GCC `-O0/-O3` and in a Clang
undefined-behavior sanitizer build. Include hashes are recorded because each
source deliberately reuses a separately audited fixed-pool geometry. All
hashes use raw LF bytes.
