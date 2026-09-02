---
id: THM-4211
title: "Fixed-fifty cofinal two-newcomer Haar tail and eighteen-label chart"
status: >
  PROVED RELATIVE TO THM-4150/4170 + VERIFIED-EXACT + INDEPENDENTLY
  AUDITED. For the displayed thirty-label pool P, every nine-body B subset P
  becomes Haar-safe after adjoining 50 and every second newcomer r>=5682.
  The number 5682 is the exact transition only for the sufficient
  component-discrepancy activation deck, not a minimal literal newcomer
  threshold. The displayed eighteen-label chart is stronger in a different
  direction: all C(18,9)=48,620 chart bodies are safe with 50 and every
  positive outsider r. An exact (49,50) depth-eight checkpoint closes all
  C(30,9) bodies for that fixed pair. THM-4150 transfers these body statements
  to odd-tail LRC(14) families; arbitrary pair entry and LRC(14) remain OPEN.
source: codex-pair-entry-jc-mixed-20260826
depends_on:
  - THM-4150-safe-set-haar-measure-universal-odd-tail-lrc14-transfer
  - THM-4170-triple-deletion-matching-eventual-haar-odd-tail-transfer
related:
  - THM-4154-mod-six-fixed-clock-and-haar-pool-inheritance-correction
  - THM-4191-complete-full-pool-newcomer-haar-transfer
  - THM-4203-fixed-pool-seventeen-body-depth-eight-haar-completion
  - THM-4207-two-newcomer-sharp-depth-transition-base-surplus-composition-and-variable-pool-chart-number
  - THM-4214-two-newcomer-pascal-complete-eleven-body-haar-charts
  - THM-4326-lrc14-rank-two-wall-graph-complete-typed-universe-closure
  - THM-4329-lrc14-complete-thirty-label-fixed-outsider-and-thirty-two-label-pascal-chart
cofinal_primary_script: 04-computation/lrc14_fixed50_pair_and_cofinal_chart_primary_thm4211.cpp
cofinal_primary_output: 05-knowledge/results/lrc14_fixed50_depth9_cofinal_chart_primary_thm4211.out
depth_eight_output: 05-knowledge/results/lrc14_fixed50_depth8_reserve_obstruction_thm4211.out
pair_prefix_output: 05-knowledge/results/lrc14_fixed50_pair49_depth8_independent_audit_thm4211.out
pair_joint_scout_script: 04-computation/lrc14_two_newcomer_deletion_staircase_probe_20260826.cpp
pair_joint_scout_output: 05-knowledge/results/lrc14_two_newcomer_deletion_staircase_probe_20260826.out
cofinal_referee_script: 04-computation/lrc14_fixed50_depth9_cofinal_referee_audit_thm4211.cpp
cofinal_referee_output: 05-knowledge/results/lrc14_fixed50_depth9_cofinal_referee_audit_thm4211.out
literal_boundary_script: 04-computation/lrc14_fixed50_depth9_boundary_literal_audit_thm4211.cpp
literal_boundary_output: 05-knowledge/results/lrc14_fixed50_depth9_boundary_literal_audit_thm4211.out
chart_script: 04-computation/lrc14_fixed50_chart18_two_newcomer_thm4211.cpp
chart_output: 05-knowledge/results/lrc14_fixed50_chart18_two_newcomer_thm4211.out
chart_independent_script: 04-computation/lrc14_fixed50_chart18_two_newcomer_independent_audit_thm4211.cpp
chart_independent_output: 05-knowledge/results/lrc14_fixed50_chart18_two_newcomer_independent_audit_thm4211.out
cofinal_primary_script_sha256: 8ba7e706e6597c73e1b6c2f506b351591a9dfc1de79880071dfb91d610011f4b
cofinal_primary_output_sha256: 2cf53357584151f074ab58eb7642b852e468c3178773d8d85438afde36524f38
depth_eight_output_sha256: a1161ee3f6b84b85d3b122ffa85905a9053748fd4672443c485fb8c5079c407c
pair_prefix_output_sha256: 7cd424794fe4c587f0662fbbc3086b5eb6fbea9bfdd4486900754a72d012305f
pair_joint_scout_script_sha256: 53ffe1642457241153a858c3da3c20cf31e8cc79ab6009f9321895c2f20b06c2
pair_joint_scout_output_sha256: 850f67053a7bfcbe73a0e117d5a2613836e32555e60cc8b49df66eb18ba36f5c
cofinal_referee_script_sha256: 0fa6d5e07f12f079094fa8bf5ddc3abfaa09e476251142a8fd8dfd1bd151ef09
cofinal_referee_output_sha256: c5dd7af091021ff01c4d99e3100ecd13847c2e593bf35cbb79a5157903034f76
literal_boundary_script_sha256: a8c39f426de51d7deb89c84a84b6e5d22ef8ee25989eecec020a651de6e4e483
literal_boundary_output_sha256: 362fe7218578e5686cbd9b273fd663493d7b4507135dc52112f827da7298af53
chart_script_sha256: d51bdfc33b526b85d64ae79003ff9e0f6a967ed55526c29c497131833a648fe1
chart_output_sha256: 404cd6ba6d6d8241fbfd5dcf4ca83da440671dc4a797843fe39c8ccc11598b97
chart_independent_script_sha256: 821583715d1373030c0c842678f2f88553c0ec5e34ce14db1350d6e92d6f9af3
chart_independent_output_sha256: 117f59b743792489bc261412cca91d64d2ed4c32e4237d0cdb8e79caad0bcad6
hash_basis: raw LF bytes
primary_audit: >
  PASS. Fixed-prefix atom integration constructs every strict depth-nine
  reserve, computes its exact component-discrepancy activation, and exhausts
  all C(30,9) bodies at the 5681/5682 transition. A literal boundary program
  separately proves that the displayed activating edge is already lawful at
  r=5681, enforcing the scope firewall. A prior literal joint-wall scout and
  the new fixed-prefix path independently verify the (49,50) control; the
  fixed-prefix path also verifies the depth-eight reserve obstruction.
independent_audit: >
  ACCEPT. A fresh fixed-prefix literal-wall/colex reverse-incidence implementation
  reproduces the common denominator, compressed atoms, all 4,673,779 strict
  depth-nine reserves, both filtered edge counts, the 5681 cover, the 5682
  no-cover verdict, and the activating witness. The eighteen-label chart has
  a separate fixed-prefix complement-scatter audit which checks every body at
  every finite outsider, importing neither the primary zeta transform nor its
  body-filtering schedule.
---

# THM-4211 -- fixed-fifty cofinal two-newcomer Haar tail and eighteen-label chart

**PROVED RELATIVE TO THM-4150/4170 + VERIFIED-EXACT + INDEPENDENTLY
AUDITED; LRC(14) REMAINS OPEN.**

## 1. Statement and inheritance pass

Retain the fixed pool

```text
P={8,10,15,16,20,30,40,42,60,63,80,84,85,88,95,
   120,126,132,143,145,168,170,176,190,193,240,252,
   264,286,290}.                                        (1)
```

For a finite positive set `S`, put

```text
G_S={y in R/Z:min_(s in S)||sy||>=1/14},
alpha=4/63.                                              (2)
```

> **Cofinal fixed-fifty theorem.** For every integer `r>=5682` and every
> `B in binom(P,9)`,
>
> ```text
> mu(G_(B union {50,r}))>=alpha.                         (3)
> ```
>
> **Universal eighteen-label chart.** Put
>
> ```text
> C={8,16,40,42,80,84,85,88,95,
>    120,126,143,145,168,193,240,252,286}.               (4)
> ```
>
> For every positive `r notin P union {50}` and every
> `K in binom(C,9)`,
>
> ```text
> mu(G_(K union {50,r}))>=alpha.                         (5)
> ```

Since every element of `P` is below `5682`, the outsider qualification in
`(3)` is automatic. In both statements, THM-4150 gives the genuine odd-tail
consequence: for every positive integer `c` and every two distinct positive
odd integers `a,b`, some `x in R/Z` satisfies

```text
min_(v in 2c(B union {50,r}) union {a,b})||vx||>=1/14,  (6)
```

and likewise with `K` in place of `B`. These are infinite thirteen-speed
LRC(14) families, not arbitrary entry and not full LRC(14).

There is also a finite exact checkpoint:

```text
mu(G_(B union {49,50}))>=alpha  for every B in binom(P,9). (7)
```

The closest proved mechanisms are THM-4191's complete one-newcomer deck and
THM-4207's exact joint deck for the fixed pair `(50,51)`. The present theorem
changes the operation: freeze the first newcomer, retain each repair's base
mass and component count, and activate it uniformly in the second newcomer.
THM-4207 and `(3)` therefore have different quantifiers and different lawful
edge predicates.

The canonical hostile is a depth-eight reserve cover whose body is itself
strictly safe in the limiting calculation. The corrected near miss is to
read a sufficient-deck cover as body danger. The least-used relevant sidecar
is the number of connected components of the unprojected safe set. The live
concept board was

```text
exact joint deck | strict reserve deck | activation filtration
failure singletons | restricted pool chart | literal boundary. (8)
```

## 2. Component-discrepancy activation deck

For `R in binom(P,9)`, define the fixed-first-newcomer base

```text
U_R=G_((P union {50})\R),
mu(U_R)=m_R/D_50,
D_50=91,205,797,082,400,                             (9)
```

and let `c_R` be the number of connected components of `U_R` on the circle.
THM-4170 proves, for every positive integer `r`,

```text
mu(U_R intersect G_{r})
 >=(6/7)(m_R/D_50)-6c_R/(49r).                         (10)
```

Call `R` a strict reserve when

```text
s_R=27m_R-2D_50>0.                                    (11)
```

This is exactly the condition that the first term on the right of `(10)` is
strictly larger than `alpha`. Define its sufficient activation time by

```text
kappa(R)=ceil(27c_R D_50/(7s_R)),                      (12)
```

and the cutoff-filtered deck

```text
A_9(Q)={R in binom(P,9):s_R>0 and kappa(R)<=Q}.         (13)
```

Equations `(10)--(12)` give the load-bearing implication

```text
r>=Q and R in A_9(Q)
 => mu(G_((P union {50,r})\R))>=alpha.                 (14)
```

Thus `A_9(Q)` is a sufficient subdeck of the literal joint deck for every
second newcomer at least `Q`. It need not contain every literally lawful
repair. Keeping that distinction is essential in Section 4.

The exact connection contract is

```text
source:       base safe set U_R at fixed newcomer 50
target:       a lawful repair for every r>=Q
map:          (m_R,c_R) -> kappa(R)
preserved:    labelled deletion R, base mass, component count, threshold
destroyed:    component addresses and r-specific phase alignment
sidecar:      literal joint-wall mass at the boundary
decisive test: tau(A_9(Q))>9.                           (15)
```

## 3. Exact transition at 5682

The primary exact calculation constructs all strict reserves and obtains

```text
|{R in binom(P,9):s_R>0}|=4,673,779,
threshold equalities=0,
min kappa=730,
max kappa=18,323,165,763.                              (16)
```

At the adjacent cutoffs it finds

| `Q` | `|A_9(Q)|` | exact nine-body verdict |
|---:|---:|:---|
| `5681` | `2,010,675` | has a cover |
| `5682` | `2,011,013` | no cover among all `14,307,150` bodies |

One cover at `5681` is

```text
B#={16,85,88,95,143,168,193,240,290}.                 (17)
```

At `5682`, the exhaustive scan proves

```text
tau(A_9(5682))>9.                                     (18)
```

In particular, every `B in binom(P,9)` misses a deletion
`R in A_9(5682)`. For every `r>=5682`, `(14)` makes that deletion lawful and

```text
B union {50,r} subset (P union {50,r})\R,
G_((P union {50,r})\R) subset G_(B union {50,r}).      (19)
```

Safe-set monotonicity proves `(3)`.

The edge which first destroys the displayed cover is

```text
R#={40,63,80,120,126,145,176,190,252},
s_(R#)=10,402,250,150,088,
c_(R#)=168,
kappa(R#)=5682.                                       (20)
```

The primary ordering uses `1,410,296,496` incidence checks at `5682`. The
independent referee orders edges differently and uses `1,361,570,163`; it
nevertheless reproduces every count and witness in `(16)--(20)`. This
difference is a positive control that the verdict is not a frozen traversal
trace.

## 4. Why 5682 is not a literal entry threshold

The activation time `(12)` is obtained from a worst-case component
discrepancy bound. It is sufficient, not necessary. A separate literal
joint-wall integration of the edge `R#` gives

```text
r=5681:
D=699,244,444,298,400,
m=46,923,230,464,562,
63m-4D=159,185,742,073,806>0;                          (21)

r=5682:
D=28,790,629,945,677,600,
m=1,935,981,420,245,616,
63m-4D=6,804,309,692,763,408>0.                        (22)
```

So `R#` is already literally lawful at `5681`, even though it is absent from
`A_9(5681)`. The cover `(17)` therefore proves only

```text
tau(A_9(5681))<=9.                                    (23)
```

It does not prove a cover of the literal joint deck, an unsafe body, or
failure of `(3)` at `5681`. The strongest exact boundary statement is:

```text
5682 is minimal for the activation-filtered certificate (13);
the minimal literal cofinal threshold remains OPEN.                    (24)
```

## 5. The depth-eight obstruction is a certificate obstruction

For `R in binom(P,8)`, again put

```text
U_R=G_((P union {50})\R),        mu(U_R)=m_R/D_50,
L_8={R:27m_R-2D_50>0}.                                 (25a)
```

Even the complete strict-limit reserve deck has

```text
|L_8|=576,900, threshold equalities=0,                 (25)
```

and the nine-cover

```text
B*={16,85,88,95,143,145,168,193,240}.                 (26)
```

The best depth-eight deletion candidate disjoint from `B*` has

```text
27m_R-2D_50=-2,181,834,281,286<0,                     (27)
```

so it is not a member of `L_8`. Define instead

```text
mu(G_(B* union {50}))=M_(B*)/D_50.                    (27a)
```

Direct integration of this body gives the positive limiting margin

```text
54M_(B*)-4D_50=573,504,307,025,796>0.                 (28)
```

Thus depth eight fails robustly as a universal reserve-deck proof while its
hostile body points in the safe direction. The first failed implication is

```text
body covers the reserve deck  !=>  body is unsafe.     (29)
```

The strongest survivor is the sharp activation-deck transition at depth
nine. This obstruction is distinct from THM-4207's exact fixed-pair equality
`tau(E_8(50,51))=8`.

## 6. An eighteen-label chart valid for every outsider

For `K in binom(C,9)`, put

```text
V_K=G_(K union {50}),
mu(V_K)=n_K/D_50,
d_K=#components(V_K).                                 (30)
```

The same discrepancy lemma gives

```text
mu(V_K intersect G_r)
 >=(6/7)(n_K/D_50)-6d_K/(49r).                         (31)
```

All

```text
binom(18,9)=48,620                                    (32)
```

bases have strict limiting reserve. The least exact reserve is

```text
54n_K-4D_50=517,215,373,867,152>0                     (33)
```

at

```text
K={8,80,85,88,95,143,145,168,193}.                   (34)
```

The individual discrepancy cutoffs range from `124` through `448`. The
primary path performs the exact literal joint-wall scan only below each
body's own cutoff: `11,592,512` comparisons over the `416` lawful outsiders
below `448`, with zero failures and zero threshold equalities. The closest
normalized finite row is

```text
r=96,
K={8,80,88,95,120,145,168,193,286},
63mu(G_(K union {50,96}))-4
 =407,867,606,407/85,159,474,400>0.                   (35)
```

The independent implementation takes the hostile route: it checks every one
of the `48,620` bodies at every one of those `416` outsiders, for
`20,225,920` literal comparisons, without the primary path's upward zeta
transform or cutoff pruning. It again finds zero failures and zero
equalities, and directly checks the boundary `r=448`. Equations `(31)--(35)`
prove `(5)` for every positive outsider.

Define the fixed-fifty universal chart number by

```text
chi_50=max{|C'|:C' subset P and
  mu(G_(K union {50,r}))>=alpha
  for every K in binom(C',9) and every positive
  r notin P union {50}}.                               (36)
```

The proved conclusion is exactly

```text
chi_50>=18.                                           (37)
```

No maximality claim is made. The chart contains seven multiples of `6`, so
only `binom(11,9)=55` of its nine-bodies avoid all multiples of `6`.
Consequently `48,565` chart bodies lie outside THM-4154's simple fixed-clock
hypothesis before the second newcomer is even inspected. This is a genuine
Haar/deletion contribution, not merely a repackaging of that mod-six clock.

## 7. The fixed pair `(49,50)`

For literal joint newcomer decks, write

```text
E_d(q,r)={R in binom(P,d):
  mu(G_((P\R) union {q,r}))>=alpha}.                   (38)
```

A fixed-prefix exact integration at `(q,r)=(50,49)` obtains

```text
|E_8(50,49)|=1,536,023,
threshold equalities=0,                               (39)
```

and exhausts all `14,307,150` nine-bodies with no cover. Its closest body is

```text
{80,88,95,143,168,193,240,252,290},                   (40)
```

which still misses

```text
{8,10,16,85,132,145,170,264}.                         (41)
```

Thus `tau(E_8(50,49))>9`, and the same monotonicity argument as `(19)` proves
`(7)`. This finite checkpoint does not interpolate between `49` and the
cofinal range; literal newcomer dependence is not monotone in the label.

## 8. Exact failure-singleton sunflower lemma

The cofinal and restricted-chart proofs keep full deletion edges. A dual
entry view keeps failure atoms instead. Let `B` be a nonempty finite set of
positive labels with cardinality `m`, let `q,r` be distinct positive labels
outside `B`, and inside `G_{q,r}` define

```text
w(F)=mu({y:the set of failed labels of B is exactly F}),
F subset B.                                            (42)
```

For `p in B`, set

```text
H_p=(B\{p}) union {q,r}.                               (43)
```

Then exactly the empty and singleton-`p` atoms survive:

```text
mu(G_(H_p))=w(empty)+w({p}).                           (44)
```

For distinct `p,p'`, the two sets intersect in the common core

```text
G_(H_p) intersect G_(H_p')=G_(B union {q,r}),
mu(G_(B union {q,r}))=w(empty),                        (45)
```

and their residual petals have the disjoint masses `w({p})`. Hence

```text
sum_(p in B)mu(G_(H_p))
 =m w(empty)+sum_(p in B)w({p}).                       (46)
```

In particular, a successful one-label exchange exists if and only if

```text
w(empty)+max_(p in B)w({p})>=alpha,                   (47)
```

and the average condition

```text
m w(empty)+sum_p w({p})>=m alpha                      (48)
```

is sufficient. This is an analytic identity, not an empirical heuristic.
It also diagnoses the missing coordinate in a mass-only averaging argument:
THM-4191 lower-bounds selected exchanged masses, but does not by itself
control how much of that mass is the common empty atom and how much is a
labelled singleton petal.

## 9. Audit and replay

The primary cofinal program forms exact fixed-prefix atoms, zeta-sums every
depth-nine deletion, computes `(12)`, and performs the complete cover tests.
The independent referee rebuilds the fixed-prefix literal base walls and uses
colex reverse incidence. It imports neither the primary source nor its edge
order. The separate boundary audit alone builds the literal joint wall at
`r=5681/5682`.
The two chart programs are independent for the same reason: upward subset
zeta versus direct atom-to-body complement scatter.

Replay from the repository root:

```bash
mkdir -p /tmp/thm4211-replay

g++ -std=c++20 -O3 -DNDEBUG -Wall -Wextra -Wpedantic -Werror \
  04-computation/lrc14_fixed50_pair_and_cofinal_chart_primary_thm4211.cpp \
  -o /tmp/thm4211-replay/primary
/tmp/thm4211-replay/primary limit-chart-e9
/tmp/thm4211-replay/primary limit-chart
/tmp/thm4211-replay/primary 49 --all-covers

g++ -std=c++20 -O3 -DNDEBUG -Wall -Wextra -Wpedantic -Werror \
  04-computation/lrc14_fixed50_depth9_cofinal_referee_audit_thm4211.cpp \
  -o /tmp/thm4211-replay/referee
/tmp/thm4211-replay/referee

g++ -std=c++20 -O3 -DNDEBUG -Wall -Wextra -Wpedantic -Werror \
  04-computation/lrc14_fixed50_depth9_boundary_literal_audit_thm4211.cpp \
  -o /tmp/thm4211-replay/boundary
/tmp/thm4211-replay/boundary

g++ -std=c++20 -O3 -DNDEBUG -Wall -Wextra -Wpedantic -Werror \
  04-computation/lrc14_fixed50_chart18_two_newcomer_thm4211.cpp \
  -o /tmp/thm4211-replay/chart-primary
/tmp/thm4211-replay/chart-primary

g++ -std=c++20 -O3 -DNDEBUG -Wall -Wextra -Wpedantic -Werror \
  04-computation/lrc14_fixed50_chart18_two_newcomer_independent_audit_thm4211.cpp \
  -o /tmp/thm4211-replay/chart-independent
/tmp/thm4211-replay/chart-independent
```

Compare each stdout stream with its matching frozen output. The cofinal
primary source emits three separate frozen streams for the depth-nine chart,
the depth-eight hostile, and `(49,50)`.

## 10. Strict scope

This theorem proves:

1. the cofinal fixed-`50` assertion `(3)` for every full-pool nine-body;
2. the all-outsider chart assertion `(5)` on the displayed eighteen labels;
3. the fixed literal-pair checkpoint `(7)`;
4. the exact failure-singleton identities `(44)--(48)`; and
5. the THM-4150 odd-tail transfers of `(3)`, `(5)`, and `(7)`.

It does **not** prove that `5682` is a minimal literal threshold, safety for
every second newcomer on the full pool, an arbitrary two-newcomer theorem,
entry of an arbitrary LRC(14) instance into these families, arbitrary
thirteen-speed safety, or LRC(14).

**Forward consolidation (2026-09-01).** THM-4326 later proves the arbitrary
two-outsider theorem over this full fixed pool, and THM-4329 completes every
eleven-face of the resulting thirty-two-label Pascal chart. Those results
supersede the fixed-chart question left open here, but still do not supply
arbitrary-row entry or LRC(14).
