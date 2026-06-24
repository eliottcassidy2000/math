---
id: HYP-2979
title: LRC14 holistic route atlas and source-kernel convergence
status: PROOF-SYNTHESIS / historical route atlas and theorem target, not a proof
source: codex-2026-06-24-S159
script: 04-computation/lrc14_holistic_route_atlas_codex_s159.py
result: 05-knowledge/results/lrc14_holistic_route_atlas_codex_s159.out
related:
  - HYP-2978
  - HYP-2976
  - HYP-2977
  - HYP-2975
  - HYP-2974
  - HYP-2973
  - HYP-2972
  - HYP-2971
  - HYP-2970
  - HYP-2969
  - HYP-2968
  - HYP-2967
  - HYP-2966
  - HYP-2965
  - HYP-2964
  - HYP-2963
  - HYP-2962
  - HYP-2961
  - HYP-2956
  - HYP-2955
  - HYP-2954
  - HYP-2953
  - HYP-2952
  - HYP-2951
  - HYP-2950
  - HYP-2949
  - HYP-2948
  - HYP-2947
  - HYP-2942
  - HYP-2937
  - HYP-2924
  - HYP-2920
  - HYP-2908
  - THM-358
  - THM-360
  - THM-365
  - THM-366
  - THM-379
  - THM-381
  - THM-398
  - THM-501
  - THM-523
  - THM-566
  - THM-568
  - THM-569
  - THM-571
  - THM-572
  - OPEN-Q-108
---

# HYP-2979: LRC14 Holistic Route Atlas And Source-Kernel Convergence

Integration note after rebasing over codex-2026-06-24-S160, the spectral
shadow dual route, and the Ramanujan-divisor quotient guardrail lane: HYP-2976
is now the completed lineage synthesis, HYP-2977 is the spectral-shadow dual
route, HYP-2978 is the Ramanujan quotient guardrail, and HYP-2975 is the taut
bridge graph curvature lane.  This HYP-2979 file is the computed
route-atlas/proof-kernel complement to that lineage request.

This synthesis looks back across the repo's LRC work and records the current
understanding as a single proof object.  The claim is not that LRC14 is proved.
The claim is that the many attempts have converged on a precise missing theorem:

```text
There is no primitive strict LRC14 source kernel after qdiv, exact M/Farey,
Haar/Baire boundary, endpoint-owner, C27/unital/K33, fixed-margin packet,
lift, moment, twist, and Fourier-dual data are retained.
```

Equivalently, a strict counterexample must be more than a row, more than a
tournament class, more than a scalar gap, and more than a zero Haar-open
chart.  It must be a labelled source packet that defeats every known primal
and dual certificate while avoiding the AP/Goddyn-Wong boundary atom and the
K33/H=7 state-lift endpoint.

## Computed Atlas

The audit script:

```text
04-computation/lrc14_holistic_route_atlas_codex_s159.py
```

Stored output:

```text
05-knowledge/results/lrc14_holistic_route_atlas_codex_s159.out
```

It scans theorem, hypothesis, reflection, and forum markdown artifacts with an
LRC signal.  The main corpus readout is:

```text
markdown artifacts scanned with LRC signal: 2245
canon theorem artifacts: 430
hypothesis artifacts: 638
forum post/comment artifacts: 48
```

The most frequent route families in the scan are:

```text
endpoint_core              895
moon_core_apex             660
singular_series            607
raw_tournament_guardrail   371
fixed_margin_packet        316
wide_decorrelation         314
farey_exact_m              295
ap_gw_boundary             218
moment_dual                206
c27_unital_k33             181
qdiv_small_denominator     147
boundary_moment_gk8        135
observer_source            112
haar_baire                  98
lift_packet                 45
nork_pinch                  39
fourier_toeplitz            27
twist_ladder                18
```

The top co-occurrences say the same thing qualitatively: endpoint/core,
Moon/apex, singular-series, fixed-margin packets, AP/GW, Farey, C27/K33, and
moment data are not rival threads.  They are projections of one retained
packet.

## How The Understanding Changed

### 1. Endpoint and qdiv arithmetic became the first invariant

The earliest durable theorem layer is not a computation-heavy search.  It is
endpoint arithmetic:

```text
THM-358: initial segment has exactly the unit endpoint skeleton.
THM-360: a unit endpoint a/n is protected only by a speed divisible by n.
THM-365: a counterexample needs a labelled endpoint-protection cycle.
THM-366: every counterexample must hit every small denominator 2..n.
THM-379: owner-realized endpoint cores force pressure cycles.
THM-381: LRC loneliness is observer-source reachability in a marked tournament.
```

The first lesson was already the modern lesson in miniature.  A bare endpoint
cycle is too topological and too weak; abstract cyclic arc covers exist.  A
useful cycle must carry arithmetic labels: owner speeds, exact endpoint
coordinates, strict integer inequalities, and later winding/credit data.

THM-366 becomes the qdiv gate used everywhere later:

```text
qdiv(S) <= 14  =>  direct rational witness.
strict LRC14 counterexample  =>  qdiv(S) > 14.
```

This is why AP and GW are equality atoms, not strict counterexample candidates.

### 2. Measure and density were powerful but not the closed-gap theorem

THM-398 reduced one broad route to the C-prime style multiple-of-n branch and
developed endpoint-cover, cap, and circuit-to-gap functionals.  THM-501 and
the singular-series work then gave a density-side view:

```text
L(S) = limiting covering-deficit density, controlled by exact additive
relation lattices and sinc weights.
```

This was useful because it identified AP-like relation energy as the suppressor
of witness density and explained why generic/Sidon-like rows are easy.  It also
corrected the analytic fantasy: the LRC singular series is an archimedean
singular integral, not a nontrivial Euler product over primes.

But THM-523 forced a crucial correction:

```text
positive open safe measure proves M(S)>1/14,
but zero open safe measure does not mean M(S)<1/14.
```

The conjecture is about the closed threshold gap `M(S) >= 1/14`.  AP and GW
have zero strict open mass and are still valid equality rows.  Therefore the
proof must carry closed endpoint debt and open-cover failure, not only density.

### 3. THM-523 made covering sets the true strict branch

THM-523 gave the constructive q-witness reduction:

```text
if S omits a multiple of some q in {2,...,14},
then t=1/q is a lonely point.
```

Therefore a strict counterexample must be a covering set, containing a multiple
of every `q=2,...,14`.  The local AP/GW tight rows are not covering sets,
because they have no multiple of 14.  The true strict branch is:

```text
qdiv(S) > 14.
```

THM-566 then refuted the naive denominator-ladder closure:

```text
for every fixed B, the row {1,...,11,13,84*lcm(1,...,B)}
kills every rational witness with denominator <= B.
```

The current twist-ladder work does not contradict THM-566.  It changes the
role of a finite ladder: a successful twist is a primal certificate in a
bounded packet bank; a failed ladder is a blocker hypergraph and a recursion
resource, not proof failure.

### 4. AP/GW became a boundary stratum, not two examples

HYP-2920 split AP/GW resemblance into necessary layers:

```text
q=14 punctured divisor cover,
apex unit cover,
odd skeleton and shell profile,
cofinite residue support,
maximal antipodal zsum,
literal AP complement binders,
single even dipole,
minimal one-petal acceleration,
top-petal gate 12 -> 24.
```

The residue liar `12->26` refuted residue-only proofs: it has AP-like residues
but is loose by q-threshold.  The near row `12->36` refuted coarse apex proofs:
it looks close to GW but escapes at the Farey child `3/41`.

HYP-2924 made the tournament lesson exact.  Once a vertex set, observable,
switch, and tie Hamiltonian path are fixed, one can compute achievable
tournament isomorphism classes.  But the first summit atlas also showed the
guardrail:

```text
AP and residue-liars can share a raw apex-winding tournament class.
Tournament class alone is magnitude-blind.
```

HYP-2937, HYP-2942, and HYP-2947 then refined the low frontier:

```text
C=27 shell transfer:
  AP                       perfect lower transversal
  GW 12->24                H[12:g3] -> D[3:g3]
  near/K33 12->36          H[12:g3] -> D[9:g9]
  unit petals 10/13        unit-visible hole -> unit double

q=3 unital:
  raw global GW+K33 lift fails by repeated pair;
  branch-local charts work;
  calibrated AP/GW pair-completion separates the linked K33 splice.

measurable rank:
  rank 0 = AP/GW and unit-petal discharges;
  rank 1 = nonunit K33 state-lift obligations.
```

HYP-2948 through HYP-2952 added the Haar/Baire and derived-tournament boundary
front: AP/GW are closed-boundary atoms with zero strict Haar mass, while
near/K33 and unit petals already have positive open intervals.  The six-unit
apex-pressure tournament is a necessary front filter, not a classifier.

### 5. The missing scalar became a missing functor

HYP-2953 and HYP-2954 changed the language from "find the invariant" to
"preserve the source":

```text
SourceSpec(S)
  = phase source-spectrum
    x exact M/Farey branch
    x Haar/Baire boundary status
    x C27/unital/K33/gK8 packet labels.
```

The proof object became a pullback/functor:

```text
primitive reduced residual
-> exact M/Farey branch
-> Haar-Baire open-or-boundary front
-> C27/unital/K33 owner address
-> discharge, AP/GW boundary atom, covering strictness, or TournamentStateLift.
```

This explains why so many earlier projections were useful but incomplete:

```text
raw scalar M              forgets owner/source label;
raw Haar boundary         forgets C27/GW hidden transfer;
raw C27 shell             forgets exact M/Farey scale;
raw tournament class      forgets qdiv and magnitude;
one denominator chart     forgets continuous safe intervals;
endpoint first-current    cancels on covering rows;
singular density          forgets closed threshold endpoint debt.
```

### 6. Families and sporadics became fixed-margin labelled packets

HYP-2956, HYP-2961, HYP-2962, and HYP-2963 turned the route into a classifier.
The strict counterexample grammar became:

```text
D0 qdiv witness rows                 discharged
D1 scale-separated one-large rows    discharged
D2 positive Haar-open rows           discharged
D3 unit-petal/GW-strip rows          discharged
D4 K33-labelled rows                 state-lift obligation
D5 wide/gK8-positive rows            discharged or live wide theorem

Live families:
L1 apex/lift residual
L2 genuine-wide zero-moment
L3 bounded covering core
L4 K33 zero-open state lift
L5 unnamed source kernel
```

HYP-2962 made "family" precise:

```text
family = fixed-margin labelled packet class
sporadic = bounded singleton after family parameters are removed.
```

HYP-2963 then emitted theorem-facing packets retaining:

```text
exact M,
binding denominators,
q_threshold,
Farey branch and excess,
strict Haar/Baire safe mass,
boundary debt,
C27 transfer,
S145 route/rank,
K33/state-lift flag,
covering/source family.
```

The bounded HYP-2963 bank had:

```text
audited rows          21913
below threshold       0
tight rows            2: AP, GW
unknown packets       0
Q-WITNESS             14676
BOUNDARY-AP-GW        2
BOUNDARY-PETAL        4
K33-STATE-LIFT        3
COVERING-MOMENT       7228
```

This did not prove LRC14, but it made the falsifier exact:

```text
SOURCE-SPECTRUM-UNKNOWN = 0 in the bounded bank.
```

The global theorem is to make that zero structural.

### 7. THM-571 shrank the covering branch to the Moon core

THM-568 proved the local apex-shell divisibility lemma:

```text
tight optimum denominator D has 14 | D,
and D divides an active pair sum.
```

It did not prove primitive shell collapse.  THM-569 formalized the denominator
14 unit-grid split.  THM-571 then closed the apex-majority branch:

```text
if at least seven speeds are multiples of 14,
then the row is strictly lonely,
modulo the accepted below-frontier LRC input.
```

Thus HYP-2964's Moon core is smaller:

```text
qdiv(S) > 14,
top-balanced,
|S cap 14Z| <= 6,
zero strict-Haar open front,
not wide/gK8 discharged,
not K33 state-lift discharged,
not fixed-margin packet exhausted.
```

HYP-2965 through HYP-2969 explored this core in several ways:

```text
boundary-gap packets:
  1187 qdiv>14 covering rows audited, all positive-open, zero zero-open.
  first endpoint current cancels, so raw divergence is false.

NORK pinch atlas:
  141351 exact qdiv>=14 rows, no non-AP/GW F6/NORK packet.

few-apex lift packet:
  8190 qdiv>14 structured rows, no zero lift mass.

boundary-moment ledger:
  one all-covered denominator chart is not an obstruction;
  covering rows can be chart-covered and still Haar-open positive.
```

This is the local positive-front side of the current route.

### 8. The newest dual routes converge on the same negation

HYP-2970 through HYP-2974 are deliberately different proof angles, but they
now point at the same object.

Endpoint-credit winding:

```text
strict failure = open danger arcs cover the circle
               = positive-winding cycle in G_open(S).

dual certificate = potential Phi with Phi(b)+epsilon <= Phi(a).

AP/GW = closed zero-credit boundary cycles.
K33 = retained positive-credit state-lift exit.
```

Multiplicity/danger-count moment duals:

```text
X_S(t) or N_S(t) counts dangerous runners.
If a strict counterexample exists, count is never zero.
Admissible polynomials give exact dual certificates for safe_mu>0.
AP/GW are the only zero-safe rows in the audited banks.
Positive rows are separated by degree 7..9 count duals in the hard frontier.
```

Twist ladder:

```text
safe rational twist a/q is a primal certificate.
failed ladder is a blocker hypergraph.
q<=42 certifies the HYP-2963 bank.
q<=27 misses only lcm-tail rows, rescued by 17/41.
```

Fourier-Toeplitz PSD:

```text
strict cover implies F_S(t)=C_S(t)-1 >= 0 a.e.
all finite Toeplitz matrices of Fourier coefficients must be PSD.
a negative eigenvalue gives a trigonometric-square dual certificate.
```

These are not separate proof philosophies.  They are primal/dual shadows of
the same statement: outside AP/GW boundary equality and K33 state-lift debt,
an open cover should be impossible.

## Refuted Or Demoted Routes

The following ideas remain useful only with labels attached, or were refuted as
complete proofs.

```text
raw endpoint-cycle topology:
  all-protected abstract arc systems exist; arithmetic labels are necessary.

raw runner or residue tournament class:
  AP and residue-liars can collapse to the same class; magnitude/qdiv is lost.

fixed bounded denominator witness:
  THM-566 builds covering rows killing every denominator <= B.

raw scalar Haar mass:
  AP/GW have zero strict Haar mass but are equality, not counterexamples.

one denominator all-covered chart:
  HYP-2969 found covering rows all-covered in the chart but Haar-open positive.

endpoint first-current divergence:
  HYP-2965 found zero net first current in all audited covering rows.

shell/C27 label alone:
  shell aliases recur in loose rows unless exact M/Farey is attached.

singular-series Euler product:
  local prime-power limits reconstruct the same archimedean integral.

local apex aperture alone:
  rows like 12->84 are blocked at every denominator-14 apex and escape off-axis.
```

## Tournament Analysis

Assumption challenge: the vertices do not need to be runners.  The atlas
considered:

```text
runners,
gaps,
danger arcs,
endpoint transitions,
safe components,
fixed denominator sections,
boundary owner events,
wall crossings,
residues,
cover arcs,
Farey nodes,
Fourier modes,
moment modes,
C27 shell transfers,
unital blocks,
K33 incidence packets,
lift packets,
fixed-margin packet families,
proof obligations.
```

Chosen vertices for the S159 synthesis are proof carriers.  The pairwise
observable is the retention vector:

```text
strict predicate,
exact arithmetic,
owner/packet labels,
adaptability to unbounded families,
dual-certificate strength,
executable auditability,
anti-scalarization guard.
```

The stored tournament fingerprint is:

```text
score_hist={0:1,1:1,2:1,3:1,4:1,5:1,6:1,7:1,8:1,9:1,10:1,11:1,12:1,13:1}
directed_3cycles=0
SCC_sizes=(1,1,1,1,1,1,1,1,1,1,1,1,1,1)
Hamiltonian_path_count=1
```

Hamiltonian path:

```text
endpoint_winding_dual
> boundary_gap_lift_packet
> twist_ladder_blocker
> fixed_margin_labelled_packet
> source_spectrum_pullback
> qdiv_gate
> danger_count_moment_dual
> c27_unital_k33_state
> exact_M_Farey
> haar_baire_boundary
> fourier_toeplitz_psd
> singular_series_relation_lattice
> raw_tournament_class
> raw_scalar_mass
```

This transitivity is not a proof.  It is a retention-order statement: under
this conservative gauge, the modern endpoint/packet/dual objects preserve more
of the LRC counterexample predicate than raw tournaments or scalar summaries.

## Current Necessary Conditions For A Strict Counterexample

A strict LRC14 counterexample must satisfy all of these, simultaneously:

1. `qdiv(S)>14`; otherwise THM-523 gives a rational witness.
2. `|S cap 14Z|<=6`; otherwise THM-571 discharges the apex-majority branch.
3. It is top-balanced after one-large peeling; otherwise HYP-2906/THM-398 style interval dodge applies.
4. Its strict regular-open safe set has zero Haar measure.
5. It is not AP/GW, because AP/GW have `qdiv=14` equality rather than strict failure.
6. It has no positive boundary bridge in the endpoint arrangement.
7. Its endpoint transition graph has an open unit-winding cycle and no Farkas potential.
8. The open cycle cannot be a closed zero-credit AP/GW circuit.
9. Any positive-credit nonunit cycle must avoid constructing the K33/H=7 state-lift.
10. It must defeat admissible multiplicity and danger-count moment duals.
11. It must defeat all relevant adaptive rational twist witnesses, while its blocker hypergraph gives no descent or next-rung explanation.
12. It must satisfy the Fourier-Toeplitz PSD necessary conditions for `C_S(t)-1>=0`.
13. It must not have a positive few-apex lift packet.
14. It must not have a positive gK8/L_y boundary-moment image.
15. It must not lie in a fixed-margin labelled packet component connected to a known discharge.
16. It must not expose a new Johnson-harmonic sector that has already been named F7; if it does, that sector becomes the next theorem/falsifier.

This stack is intentionally overdetermined.  A row that satisfies all of it is
the first genuinely new object the repo is hunting:

```text
a primitive qdiv>14, few-apex/top-balanced, zero-open,
non-AP/GW, non-K33-lifted, moment-invisible, twist-blocked,
PSD-compatible, fixed-margin source-spectrum unknown kernel.
```

No such row is known.

## Candidate Breakthrough Theorem

The clean theorem target is:

```text
LRC14 Source-Kernel Exclusion Theorem.

Every primitive 13-speed row S emits a labelled source packet P(S).
If S is strict-bad, then P(S) lies in the qdiv>14 Moon core.
Inside that core, one of the following holds:

1. endpoint/lift/twist/moment/Fourier dual data gives a positive witness;
2. P(S) is a closed AP/GW boundary atom, impossible because qdiv>14;
3. P(S) carries a K33/H=7 state-lift debt, contradicted by THM-572 after lift construction;
4. P(S) exposes a new fixed-margin Johnson-harmonic sector F7.

The theorem is proved once case 4 is ruled out or classified and discharged.
```

The summit is therefore not "find one more clever invariant."  It is proving
that the retained packet map is complete enough that every quotient has a named
exit.

## Next Work

Highest-value next steps:

1. Implement HYP-2974 and compare Toeplitz PSD-blind rows against HYP-2971,
   HYP-2972, and HYP-2973.
2. Build an endpoint-credit graph emitter for the HYP-2963 bank and record
   positive-winding SCC fingerprints, if any.
3. Add danger-count degree, twist witness, and endpoint-potential fields to
   the fixed-margin packet classifier.
4. Search directly for the overdetermined source-kernel row described above.
5. Try to prove a meta-dual theorem: any row defeating one dual must expose the
   structure that makes another dual succeed, unless it is AP/GW or K33.

That last step is the likely breakthrough shape: not one certificate, but a
finite incompatibility theorem among the ways certificates can fail.
