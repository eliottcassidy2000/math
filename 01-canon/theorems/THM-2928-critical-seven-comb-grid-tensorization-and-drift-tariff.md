---
id: THM-2928
title: Critical seven-comb grid tensorization and drift tariff
status: >
  PROVED.  On a nonempty 1/L-grid carrier, aligned danger combs factor
  through the normalized multiplier circle.  Seven aligned combs leave at
  least 15/154 of the carrier.  After k aligned combs, every pointwise
  (7-k)-comb drift cover pays a strict reciprocal tariff controlled by the
  normalized safe mass and tooth-component count.  For six aligned combs,
  the remaining speed obeys the 39/61 cone and an exact finite
  phase-address intersection.  On a literal six-body carrier, the two
  six-comb safe floors collapse that intersection to finitely many clocks,
  all of which fail; hence the fully aligned and one-drift branches are
  empty.  A component-free phase-load tariff bounds one slope in every
  k>=2 mixed branch, and every fixed five-aligned/two-drift chart reduces
  to an exact finite pair-clock quotient with the sharper second-slope
  bound c_2/d_2<13 max A.  At the 2/25 deletion threshold, either that
  bound improves to c_2/d_2<(25/3) max A or all five aligned multipliers
  lie in an explicit body-uniform finite staircase.  An independent body-
  projection sieve kills 240,560 of 251,536 body/divisor rows.  Its fixed-
  phase address capacity reduces 3,066,274 denominator-pair occurrences to
  23,755.  An exact relaxed arithmetic-progression mask kills all 22,813
  survivors having a denominator at most 1,000; an arbitrary-class load
  relaxation kills 940 of the remaining 942; and a quotient-fiber
  transversal obstruction kills the last two diagonal pairs.  Hence the
  literal five-aligned/two-drift branch is empty uniformly.  At four
  aligned combs and three drifts, support and arbitrary-class loads reduce
  298,255,882 denominator-triple occurrences to 544,571.  All-divisor
  fibre overload and common-phase activity reduce this to 29,364; an exact
  eight-state threshold transport over the divisor tree leaves 19
  occurrences; and local unit-needle set cover on Z/99Z and Z/98Z kills
  all 19.  Hence the literal four-aligned/three-drift branch is empty
  uniformly; the earlier septimal Z/49Z diagonal proof remains an
  independent structural witness.  Uniformly in k aligned and p=7-k drift
  tails, support transfer gives |S_D|/D <=(1-u_p)/u_k; an independently
  audited exact census records the surviving body/divisor rows and a
  divisor-Mobius formula counts their unordered denominator multisets.
  For every upward status event, the exact real one-marginal transport
  optimum is min(1,tau), where tau is the weighted fractional-cover value
  of its minimal true clutter.  In the k=2,p=5 quotient every support-hard
  resolving modulus is septimal.  Peeling D=7q splits denominators into
  transverse one-per-fibre sections and vertical two-level spikes; retaining
  every deterministic floor and optional spike bit reduces the exact
  denominator ledger from 951,545,890,235 to 200,389,247,292 occurrences.
  Restoring located phase closes the sole-denominator-six c=4 family and
  leaves 200,141,092,521 occurrences.
  In the four-drift quotient, the same septimal split plus an exact
  one-spike superlevel law reduces 21,357,714,101 raw occurrences to
  2,548,901,482 necessary occurrences.  These are necessary relaxations,
  not sector closures.  The critical
  two-torus carrier integral is positive on every fixed seven-tail affine ray;
  uniformly, seven tails in one common canonical-ruler quotient block are
  safe once the quotient is at least 315,586.  Branches with four or more
  drifts remain open; this does not prove LRC(14).
source: root-2026-07-29 with independent hostile audits by seven-wall-tensor-audit, critical-residue-tree, mixed-residual-atlas, mixed-overlap-bound, and mixed-lorenz-sieve
depends_on:
  - THM-594-pair-overlap-law-mirsky-newman-floor
  - THM-1094-exact-two-comb-component-theorem
  - THM-1166-seven-wall-fano-gcd-discrepancy
  - THM-1221-seven-wall-strict-spectrum-hunter-floor
  - THM-1234-sharp-five-comb-compatibility-floor
  - THM-2162-signed-endpoint-cocycle-and-bv-component-split
  - THM-2182-endpoint-grid-product-and-tail-overlap-sidecar
  - THM-2184-two-scale-tail-continuation-profile
  - LRC(<=13)
related:
  - THM-1135-r6-harmonic-tail-finite-box
  - THM-2184-two-scale-tail-continuation-profile
  - THM-1132-sharp-measure-horn-constant-dissolves-r6-wall
  - THM-1176-seven-wall-slow-gap-harmonic-crowding
  - THM-2893-complement-cap-finite-core-flag-lemma
  - THM-2941-critical-seven-slot-scalar-wall-and-balanced-boundary
  - HYP-7678
  - HYP-7870
script: 04-computation/lrc14_critical_one_drift_clock_thm2928.py
output: 05-knowledge/results/lrc14_critical_one_drift_clock_thm2928.out
support_script: 04-computation/lrc14_two_drift_body_projection_support_thm2928.py
support_output: 05-knowledge/results/lrc14_two_drift_body_projection_support_thm2928.out
address_script: 04-computation/lrc14_two_drift_relaxed_address_mask_thm2928.py
address_output: 05-knowledge/results/lrc14_two_drift_relaxed_address_mask_thm2928.out
address_completion_script: 04-computation/lrc14_two_drift_arbitrary_class_fiber_thm2928.py
address_completion_output: 05-knowledge/results/lrc14_two_drift_arbitrary_class_fiber_thm2928.out
three_drift_script: 04-computation/lrc14_three_drift_body_projection_fiber_thm2928.py
three_drift_output: 05-knowledge/results/lrc14_three_drift_body_projection_fiber_thm2928.out
three_drift_diagonal_script: 04-computation/lrc14_three_drift_diagonal_crt_thm2928.py
three_drift_diagonal_output: 05-knowledge/results/lrc14_three_drift_diagonal_crt_thm2928.out
three_drift_activity_script: 04-computation/lrc14_three_drift_mixed_activity_fiber_thm2928.py
three_drift_activity_output: 05-knowledge/results/lrc14_three_drift_mixed_activity_fiber_thm2928.out
three_drift_lorenz_script: 04-computation/lrc14_three_drift_mixed_lorenz_activity_thm2928.py
three_drift_lorenz_output: 05-knowledge/results/lrc14_three_drift_mixed_lorenz_activity_thm2928.out
three_drift_terminal_script: 04-computation/lrc14_three_drift_threshold_transport_terminal_thm2928.py
three_drift_terminal_output: 05-knowledge/results/lrc14_three_drift_threshold_transport_terminal_thm2928.out
three_drift_terminal_audit_script: 04-computation/lrc14_three_drift_terminal_local_incidence_audit_thm2928.py
three_drift_terminal_audit_output: 05-knowledge/results/lrc14_three_drift_terminal_local_incidence_audit_thm2928.out
support_transfer_ladder_script: 04-computation/lrc14_aligned_drift_support_transfer_ladder_thm2928.py
support_transfer_ladder_output: 05-knowledge/results/lrc14_aligned_drift_support_transfer_ladder_thm2928.out
support_transfer_ladder_audit_script: 04-computation/lrc14_aligned_drift_support_transfer_ladder_independent_audit_thm2928.py
support_transfer_ladder_audit_output: 05-knowledge/results/lrc14_aligned_drift_support_transfer_ladder_independent_audit_thm2928.out
upward_status_cover_script: 04-computation/lrc14_upward_status_fractional_cover_audit_thm2928.py
upward_status_cover_output: 05-knowledge/results/lrc14_upward_status_fractional_cover_audit_thm2928.out
k2_septimal_script: 04-computation/lrc14_k2_septimal_floor_exception_gf_thm2928.py
k2_septimal_engine: 04-computation/lrc14_k2_septimal_floor_exception_engine_thm2928.cpp
k2_septimal_output: 05-knowledge/results/lrc14_k2_septimal_floor_exception_gf_thm2928.out
k2_d6_located_script: 04-computation/lrc14_k2_d6_located_phase_closure_thm2928.py
k2_d6_located_output: 05-knowledge/results/lrc14_k2_d6_located_phase_closure_thm2928.out
k3_septimal_base_script: 04-computation/lrc14_k3_four_drift_divisor_status_gf_thm2928.py
k3_septimal_q7_script: 04-computation/lrc14_k3_four_drift_q7_all_D_gf_thm2928.py
k3_septimal_script: 04-computation/lrc14_k3_four_drift_expected_spike_gf_thm2928.py
k3_septimal_output: 05-knowledge/results/lrc14_k3_four_drift_expected_spike_gf_thm2928.out
k3_septimal_audit_script: 04-computation/lrc14_k3_four_drift_expected_spike_literal_audit_thm2928.py
k3_septimal_audit_output: 05-knowledge/results/lrc14_k3_four_drift_expected_spike_literal_audit_thm2928.out
affine_profile_script: 04-computation/lrc14_j7_affine_profile_min_frequency_thm2928.py
affine_profile_output: 05-knowledge/results/lrc14_j7_affine_profile_min_frequency_thm2928.out
---

# THM-2928 -- critical seven-comb grid tensorization and drift tariff

> **INDEPENDENT SECOND PROOF.**  THM-2941 closes the same literal
> five-aligned/two-drift face by the lossless projected residual
> `P_(E,Z)=phi_L(C_E minus union_(z in Z)D_z)`.  That drift-first projection
> is dual to the aligned-safe/body-address projection below and gives a
> separate exact reconstruction.

## 1. Grid carrier and the inherited drift chart

Write

```text
T=R/Z,                 D_w={t in T: ||wt||<1/14}.       (1)
```

Fix a positive integer `L`, a nonempty set of cell addresses
`J subset Z/LZ`, and the carrier

```text
C_J=union_(j in J) (j/L,(j+1)/L),       h=mu(C_J)=|J|/L. (2)
```

Grid endpoints never affect the measure statements.  Once at least one
aligned comb has been removed, they do not affect the pointwise statements
either: every grid endpoint belongs to every `D_(La)`.

For an arbitrary positive speed, write

```text
w=La+b,                 0<=b<L.
```

On cell `j`, put `t=(j+u)/L`.  Then

```text
||wt||=||(a+b/L)u+bj/L||.                               (3)
```

This is the `Q=L` specialization of the exact midpoint-cell coordinate in
[THM-2184, two-scale tail continuation profile][thm2184].  It identifies the
mixed-residue datum precisely:

```text
slope detuning:       b/L;
cell phase word:      bj mod L.                         (4)
```

When `b=0`, both disappear.  The aligned locus is therefore a genuine fixed
locus of the cell transport, not a scalar approximation to a general tail.

[thm2184]: THM-2184-two-scale-tail-continuation-profile.md

## 2. Boolean endpoint-grid tensorization

Let `A` be any finite set of positive integers.  For every bounded function
`F` of the Boolean danger word,

```text
integral_(C_J) F((1_(D_(La))(t))_(a in A)) dt
 =h integral_T F((1_(D_a)(u))_(a in A)) du.             (5)
```

Indeed, on each selected cell,

```text
La(j+u)/L=aj+au=au mod 1,
```

and `dt=du/L`.  Every selected cell therefore carries the same complete
Boolean word, and summing its `1/L` copy over `|J|` cells proves `(5)`.

For a product of safe indicators, `(5)` is exactly the endpoint-grid product
mechanism of
[THM-2182, endpoint-grid product][thm2182].  Equation `(5)` records its
Boolean extension; it is not a claim that the base product mechanism is new.

[thm2182]: THM-2182-endpoint-grid-product-and-tail-overlap-sidecar.md

Put

```text
R_A=T minus union_(a in A)D_a,       u_A=mu(R_A).        (6)
```

Then

```text
mu(C_J minus union_(a in A)D_(La))=h u_A,                (7)
```

and, inside every selected grid cell, the residual is a copy of `R_A`
scaled by `1/L`.  This is the exact toothpick self-similarity used below.

## 3. Tooth count and its scale law

Suppose `A={a_1,...,a_k}` has distinct entries and `1<=k<=7`.  Let

```text
N_A=sum_i a_i
```

and let `nu_A` be the number of **all** connected components of `R_A`,
including possible singleton components.  In every case used below `u_A>0`:
for `k<=6` this follows from the union bound, and for `k=7` it follows from
Section 4.

Each `D_(a_i)` has `a_i` open teeth.  The `k` teeth containing zero all
overlap, so those `k` components contribute only one component to their
union.  Further mergers can only lower the count.  Hence

```text
1<=nu_A<=N_A-k+1.                                      (8)
```

The common origin tooth also means every component of `R_A` has a
nonwrapping representative inside `(0,1)`.

For every positive integer `d`,

```text
R_(dA)={t:dt mod 1 belongs to R_A},
u_(dA)=u_A,
nu_(dA)=d nu_A.                                        (9)
```

The first identity is immediate from `D_(da)={t:dt in D_a}`; the degree-`d`
circle covering gives the other two.  Thus the tooth functional

```text
H(A)=7k u_A/(6nu_A)                                    (10)
```

has the exact homogeneity

```text
H(dA)=H(A)/d.                                          (11)
```

For one multiplier, `u_{a}=6/7`, `nu_{a}=a`, and therefore

```text
H({a})=1/a.                                            (12)
```

This is the recursive functional form behind the ordinary slow-gap harmonic
invoice.

## 4. The critical seven-comb fixed locus is quantitatively safe

For every seven distinct positive integers `A`, the global strict-spectrum
Hunter theorem
[THM-1221][thm1221] gives

```text
u_A>=15/154.                                           (13)
```

Combining `(7)` and `(13)`,

```text
mu(C_J minus union_(a in A)D_(La))>=15h/154>0.          (14)
```

Thus seven distinct aligned combs cannot cover any nonempty grid carrier.
More locally, every selected cell contains an aligned-safe interval of
length at least

```text
15/(154L nu_A)>=15/(154L(N_A-6)).                      (15)
```

The independent divisor-minimal Fourier argument of
[THM-594][thm594] gives the weaker but purely spectral floor

```text
u_A>=2sin^2(pi/7)/(7pi^2).                             (16)
```

The stronger value `(13)`, not `(16)`, is used in `(14)`--`(15)`.

[thm1221]: THM-1221-seven-wall-strict-spectrum-hunter-floor.md
[thm594]: THM-594-pair-overlap-law-mirsky-newman-floor.md

There is an important equality boundary.  For every positive integer `a`,
`(5)` gives

```text
mu(C_J intersect D_(La))=h/7.                          (17)
```

Consequently the sum of seven aligned singleton coverages is exactly `h`.
The strict scalar first-apex gate at the `p=7` wall really does stop at
equality, as anticipated by
[THM-2893][thm2893].  What repairs the aligned branch is the joint Boolean
law and its positive safe defect `(13)`, not a stricter singleton estimate.

[thm2893]: THM-2893-complement-cap-finite-core-flag-lemma.md

## 5. The strict reciprocal tariff after aligned deletions

Let `1<=k<=6`, let `1<=r<=6`, and suppose positive speeds

```text
z_1,...,z_r
```

cover the aligned residual pointwise:

```text
C_J minus union_(a in A)D_(La)
 subset union_(q=1)^r D_(z_q).                         (18)
```

Then the general tariff is

```text
sum_(q=1)^r 1/z_q
 >7(7-r)u_A/(6Lnu_A).                                 (19)
```

At the critical seven-wall cardinality `r=7-k`, this becomes

```text
sum_(q=1)^(7-k) 1/z_q
 >H(A)/L
 =7k u_A/(6Lnu_A)                                     (20)
 >=k(7-k)/(6Lnu_A)
 >=k(7-k)/(6L(N_A-k+1)).                              (21)
```

### Proof

Choose a largest positive-length component `K` of `R_A`.  It is a compact
interval with

```text
|K|>=u_A/nu_A.                                        (22a)
```

In any selected grid cell, its copy `I` is a compact interval of length

```text
ell=|K|/L>=u_A/(Lnu_A).                               (22b)
```

The sharp periodic discrepancy estimate in
[THM-1094, equation (10)][thm1094] says, for every interval `I`,

```text
mu(I intersect D_z)<=ell/7+6/(49z).                    (23)
```

The union of the `r` danger combs is open and contains the compact interval
`I`.  Hence a slightly longer interval `I'` of length `ell'>ell` is still
covered.  Applying `(23)` to `I'` gives

```text
ell'<=r ell'/7+(6/49)sum_q 1/z_q.
```

Rearranging and using `(22b)` proves the strict inequality
`(19)`.  For `r=7-k`, the union bound
`u_A>=1-k/7=(7-k)/7` and `(8)` give `(20)`--`(21)`.  QED.

For an almost-everywhere cover, the compact enlargement is unavailable and
the same proof gives the corresponding non-strict inequality.

[thm1094]: THM-1094-exact-two-comb-component-theorem.md

At the critical cardinality `r=7-k`, at least one drift speed satisfies

```text
min_q z_q
 <6rLnu_A/(7ku_A)
 <=6Lnu_A/k
 <=6L(N_A-k+1)/k.                                     (24)
```

For `k=1`, equations `(12)` and `(19)` recover exactly the strict slow-gap
pressure

```text
(La) sum_(q=1)^6 1/z_q>1,                              (25)
```

the cellwise starting rung of
[THM-1176, slow-gap harmonic crowding][thm1176].  For larger `k`, `(19)` is
its recursive residual version; the state is the pair `(u_A,nu_A)`, not
mass alone.

[thm1176]: THM-1176-seven-wall-slow-gap-harmonic-crowding.md

### Current all-cardinality safe floors

The crude union floor in `(21)` is not the strongest proved input.  For
distinct multiplier sets, the current uniform safe floors are

```text
k             1       2        3         4          5          6
u_A >=       6/7    66/91    55/91    558/1183   478/1365    61/273. (25a)
```

The `k=2,3,4` entries follow from the pair/triple and spanning-tree
consequences of [THM-1166][thm1166].  For `k=5`, the
`R_5>=44/273` floor of THM-1234 and the uniform spanning-tree average give
`478/1365`; its six-comb quadratic consequence gives `61/273`.
Substituting `(25a)` into `(20)` strengthens the tariff coefficient
`7ku_A/6` to

```text
c_k = 1, 22/13, 55/26, 372/169, 239/117, 61/39,       (25b)

sum 1/z_q > c_k/(Lnu_A)
             >=c_k/(L(N_A-k+1)).                      (25c)
```

[thm1166]: THM-1166-seven-wall-fano-gcd-discrepancy.md

### The component-free phase-load tariff

The largest-tooth tariff `(19)` is strongest when `nu_A` is controlled.
There is a complementary estimate that forgets components but retains the
complete load on one grid cell.  For `j in J`, put

```text
E_q(j)={u in [0,1]: ||z_q(j+u)/L||<1/14}.             (25d)
```

Apply THM-1094's interval discrepancy to the physical cell
`[j/L,(j+1)/L]` and rescale by `L`.  It gives

```text
mu(E_q(j))<=1/7+6L/(49z_q).                           (25e)
```

The proper compact set `R_A` lies in the open union of the `E_q(j)`.
Consequently that union has measure strictly larger than `u_A`, and

```text
L sum_(q=1)^r 1/z_q
 >(49/6)(u_A-r/7).                                   (25f)
```

For an almost-everywhere cover the same inequality is non-strict.  At the
critical cardinality `r=7-k`, the safe floors `(25a)` give

```text
k                 1      2      3        4         5        6
b_k              0     7/78   7/26    119/338   308/585   77/117

L sum 1/z_q > b_k.                                    (25g)
```

For `k>=2`, some drift therefore satisfies

```text
min_q z_q/L < (7-k)/b_k.                              (25h)
```

In particular, five aligned combs and two drifts force

```text
min(z_1,z_2)/L<585/154.                               (25i)
```

This is a slope bound independent of `nu_A`; it does not control the
denominator of that rational slope uniformly as `L` varies.  For `k=1`,
the coefficient vanishes and the component-sensitive tariff `(19)` remains
the live estimate.

## 6. Six aligned combs and one drift: the `39/61` cone

Now take `k=6`, `r=1`, and write the aligned speeds as

```text
w_i=La_i.
```

The six-comb consequence of
[THM-1234][thm1234] gives

```text
u_A>=61/273.                                           (26)
```

Equations `(19)` and `(26)` imply

```text
z<Lnu_A/(7u_A)
 <=(39/61)Lnu_A
 <=(39/61)L(N_A-5)                                    (27)
 =(39/61)(sum_i w_i-5L).
```

In integral form,

```text
z<=ceil(39L(N_A-5)/61)-1.                              (28)
```

This is a strict additive cone for the one-drift branch.  It is only a
uniform corollary: using the exact `(u_A,nu_A)` in the first inequality of
`(27)` can be substantially stronger.

[thm1234]: THM-1234-sharp-five-comb-compatibility-floor.md

## 7. Exact one-drift phase and gcd sidecar

Continue with `k=6`.  List **all** connected components of `R_A`, including
any singleton components, as

```text
K_s=[m_s-rho_s,m_s+rho_s] subset (0,1),       rho_s>=0. (29)
```

Then the one drift comb covers the entire aligned residual pointwise if and
only if

```text
||z(j+m_s)/L||+z rho_s/L<1/14
             for every j in J and every s.             (30)
```

Indeed, the copy `(j+K_s)/L` is connected and compact.  It lies in `D_z`
exactly when it lies strictly inside one open `z`-tooth; comparing its
centre and half-length gives `(30)`.

Equation `(30)` has a finite address form.  On the circle `R/LZ`, define

```text
A_z=intersection_s
 {x: ||(x+z m_s)/L||<1/14-z rho_s/L}.                  (31)
```

An individual set in `(31)` is empty when its displayed radius is
nonpositive.  Otherwise it is an open arc.  The exact condition is

```text
zJ mod L subset A_z.                                   (32)
```

This intersection retains the relative positions of *all* residual teeth;
replacing it by one width is a lossy relaxation.

Let `rho_max` be the largest half-length of a positive component and let
`span_L` denote the length of the shortest circular arc containing a finite
residue set.  Since

```text
2rho_max>=u_A/nu_A,
```

condition `(32)` forces

```text
span_L(zJ mod L)
 <L/7-2zrho_max
 <=L/7-z u_A/nu_A.                                    (33)
```

In particular the right side must be positive.  If `g=gcd(z,L)`, the image
`zJ` lies in the subgroup `gZ/LZ` and multiplication by `z` has fibers of
size `g`.  Hence the exact and width-only cardinality screens are

```text
|J|<=g # (A_z intersect (gZ/LZ)),                      (34)

|J|<=g ceil((L/7-z u_A/nu_A)/g).                       (35)
```

Thus the metric cone `(27)` does not leave an unstructured finite box.  It
leaves a simultaneous rational tooth-address intersection with an explicit
gcd compression.

There is a useful density form of this compression.  Put

```text
g=gcd(z,L),       d=L/g,       c=z/g,       P=J mod d. (35a)
```

Then `gcd(c,d)=1`, every residue modulo `d` has exactly `g` lifts modulo
`L`, and `(30)` is equivalent to

```text
P subset B_(c,d)(A)
 ={r in Z/dZ:
   ||c(r+m_s)/d||+c rho_s/d<1/14 for every s}.         (35b)
```

The largest-component argument proving `(35)` now gives

```text
h<=|P|/d
 <=ceil(d/7-c u_A/nu_A)/d.                            (35c)
```

In particular, since `ceil(x)<x+1`,

```text
g>L(h-1/7).                                            (35d)
```

Thus every carrier with `h>1/7` forces a bounded quotient clock.  If
`h>=61/273`, then

```text
g>22L/273,                 d=L/g<273/22,              (35e)
d<=12.
```

For fixed `A`, positivity in `(35c)` also gives

```text
c/d=z/L<nu_A/(7u_A),                                  (35f)
```

so only finitely many coprime phase clocks `(c,d)` remain.

## 8. Literal six-body one-drift closure

For a fixed finite speed set `F`, put

```text
G_F=T minus union_(v in F)D_v,
L_F=lcm(14v:v in F).                                   (36)
```

Every boundary of `G_F` has the form

```text
(14q+/-1)/(14v),
```

so, modulo its finitely many endpoints, `G_F` is a union of `1/L_F` cells.
Whenever `mu(G_F)>0`, the preceding theorem applies with `L=L_F` (or with
any smaller denominator that exactly resolves the same carrier).

Now suppose `F` consists of six distinct body speeds, six of the seven tail
speeds are the distinct aligned speeds `La_i`, and the last tail speed is
`z`.  Take `J` to be the exact set of open `1/L`-cells contained in `G_F`.
The finitely many grid endpoints are all in every aligned danger comb and
therefore do not occur in the residual cover problem.

Apply the six-comb quadratic consequence of THM-1234 twice, first to `F`
and then to `A`.  It gives

```text
h=mu(G_F)>=61/273,             u_A>=61/273.            (36a)
```

Assume for contradiction that `D_z` covers the aligned residual.  Equations
`(35c)` and `(35e)` give `d<=12`, while

```text
h<=ceil(d/7-c u_A/nu_A)/d
  <=ceil(d/7)/d.                                      (36b)
```

Among `1<=d<=12`, comparison with `61/273` leaves only

```text
d in {1,2,3,4,8}.                                     (36c)
```

The carrier address set is invariant under circle reflection:

```text
j in J  iff  L-1-j in J.                              (36d)
```

For `d=2,4`, equation `(36b)` permits only one occupied residue, but
reflection sends it to `-1-r mod d` and has no fixed residue.  Hence these
two clocks are impossible.  For `d=1,3`, the unique occupied residue is,
respectively, `0,1`.  For `d=8`, density forces exactly two occupied
residues, necessarily one of the reflected pairs

```text
{r,7-r},                  0<=r<=3.                    (36e)
```

It remains to test the phase windows, not the original speed box.  For an
occupied residue `r`, `(30)` says

```text
R_A subset E_(c,d,r)
 ={u in [0,1]: ||c(u+r)/d||<1/14}.                    (36f)
```

After the change of variable `y=(u+r)/d`, THM-1094's exact interval
discrepancy gives

```text
u_A<=mu(E_(c,d,r))<=1/7+6d/(49c).                     (36g)
```

Together with `(36a)`, this forces

```text
c<=819d/539.                                           (36h)
```

Since `gcd(c,d)=1`, only the following exact endpoint table remains.  The
displayed `d=8` masses are the maxima over the four reflected pairs in
`(36e)`.

```text
d  occupied phase(s)  possible c          simultaneous window masses
1       {0}               1               1/7
3       {1}             1,2,4             0, 3/14, 3/28
8     {r,7-r}       1,3,5,7,9,11          1/7,1/21,1/35,1/49,2/63,2/77.
                                                               (36i)
```

The table follows by subdividing `[0,1]` at the rational endpoints

```text
u=(d/c)(q+/-1/14)-r.
```

Its largest entry is `3/14<61/273`, contradicting `(36a)` and `(36f)`.
Therefore:

> **Literal one-drift closure.**  Six aligned tail combs and one arbitrary
> tail comb cannot cover the safe carrier of six distinct body combs.

The fully aligned case has the stronger quantitative conclusion

```text
mu(G_F minus union_(i=1)^7 D_(La_i))
 >=(61/273)(15/154)=305/14014.                         (36j)
```

The exact rational table `(36i)` is reproduced by two independent interval
integrations in
`04-computation/lrc14_critical_one_drift_clock_thm2928.py`.  Ordinary and
`python3 -O` output must agree byte for byte with
`05-knowledge/results/lrc14_critical_one_drift_clock_thm2928.out`.
Their frozen SHA-256 hashes are

```text
script  16e7a630eb4f721836dfc95c0adf97bf9a1a252518b119450b16c71e08305dcd
output  6eb2ad0d031a94d7e85db37a30cb48a2e0abb03e4caa080f05f870a23d449746.
```

## 9. Five aligned combs and two drifts are empty

Keep the literal six-body carrier `G_F`, take its canonical resolving
denominator `L=L_F` from `(36)`, and let `J` be its full selected-cell set.
Now let `A` contain five distinct aligned
multipliers and suppose two drift combs cover the residual.  Relabel them so
that `z_1<=z_2`.  Equation `(25i)` gives

```text
z_1/L<585/154.                                        (37a)
```

Write each drift uniquely as

```text
z_q=(L/d_q)c_q,       d_q=L/gcd(z_q,L),       gcd(c_q,d_q)=1. (37b)
```

Thus `d_q` divides `L`, and `(37a)` already leaves finitely many first
clocks:

```text
d_1 divides L,             c_1/d_1<585/154.           (37c)
```

For a residue `r in S_(d_1):=J mod d_1`, define

```text
E_(c_1,d_1,r)
 ={u in [0,1]: ||c_1(u+r)/d_1||<1/14},

B_r=R_A minus E_(c_1,d_1,r).                          (37d)
```

At least one `B_r` has positive measure.  One proof uses the cited LRC
through thirteen: the partial family

```text
F union LA union {z_1}
```

has twelve nonzero speeds, so it has a time at which all distances are at
least `1/13`.  The strict `1/13-1/14` margin supplies a positive safe
interval.  That interval cannot meet a grid endpoint, where every aligned
speed vanishes, and therefore normalizes into one positive part of some
`B_r`.  Independently, the exact address/reflection audit in the companion
script proves that one drift cannot cover `R_A` even almost everywhere:
the only surviving five-aligned clocks are `(d,c)=(3,1),(8,1),(8,3)`,
whose simultaneous phase masses are respectively `0`, at most `1/7`, and
at most `1/21`, all below `478/1365`.

There is a uniform second-clock bound.  Put

```text
a_*=max A,                 N_A=sum_(a in A)a.          (37e)
```

Because five distinct positive multipliers have `a_*>=5`, equation `(37a)`
gives `z_1<5L<=La_*`.  Also `v<=L/14` for every `v in F`.  Thus `La_*` is
the largest speed in the twelve-speed partial family.  The `1/13` witness
is therefore safe at level `1/14` on a closed physical circle arc `I` of
length

```text
|I|=1/(91La_*).                                       (37f)
```

Indeed, the distance-to-integer function for a speed `v` is
`v`-Lipschitz, so the closed radius-`1/(182La_*)` arc around the witness
stays safe at level `1/14` for every member of the partial family.  Under a
pointwise countercover, every point of `I`, including both endpoints, must
therefore lie in `D_(z_2)`.  The connected arc `I` lies in one component of
that open set, whose teeth have length `1/(7z_2)`.  Closed containment in an
open tooth is strict, and hence

```text
1/(91La_*)<1/(7z_2),

c_2/d_2=z_2/L<13a_*.                                  (37g)
```

There is also a useful, but weaker, BV sidecar.  On one normalized cell,
the first drift is the integer `c_1` comb on an interval of length `1/d_1`.
THM-1094's tooth count gives

```text
number of meeting teeth
 <=c_1/d_1+8/7
 <585/154+8/7
 =5327/1078<5.                                        (37h)
```

Thus at most four first-drift teeth cut `R_A`.  Since `(8)` gives
`nu_A<=N_A-4`, every positive `B_r` has at most `N_A` positive-length BV
components.  The arc `I` cannot meet a body-grid endpoint, because every
aligned comb is dangerous there.  It therefore lies in one selected cell
and normalizes to a subinterval of some `B_r` of length `1/(91a_*)`.
Writing the mass and positive-length BV component count of that `B_r` as
`mu_r` and `K_r`, we have `mu_r>=1/(91a_*)` and `K_r<=N_A`.  Applying
THM-2162 to `(r_2+B_r)/d_2` gives

```text
mu_r/d_2<=mu_r/(7d_2)+6K_r/(49c_2),
```

and recovers the older almost-everywhere estimate

```text
c_2/d_2<=13N_Aa_*,
```

but the connected closed-arc argument `(37f)`--`(37g)` is stronger for the
pointwise problem.  Hence both clocks lie in the explicit finite box

```text
d_1,d_2 divide L,
gcd(c_q,d_q)=1,
c_1/d_1<585/154,
c_2/d_2<13 max A.                                     (37i)
```

Row-specific exact values of `K_r/mu_r` retain additional BV information,
although the universal pointwise face `(37g)` is usually stronger.

On the direct exactly-six-in-window LRC(14) boundary,
`F subset {1,...,14}`.  In that case

```text
L_F divides L_14:=14 lcm(1,...,14)=5,045,040.          (37ia)
```

Thus all possible clock denominators, and the first numerator through
`(37c)`, already lie in one body-uniform finite set.  The remaining
unbounded scale is `max A`; its projective high slope
`c_2/(d_2 max A)` lies in the bounded interval `(0,13)` and remains
coupled to the five-multiplier shape.  After the fully aligned and
one-drift branches have been discharged, a genuinely two-drift row also
has `d_1,d_2>1`; a clock with `d_q=1` is aligned and belongs to an earlier
branch.

### The `2/25` deletion dichotomy and reciprocal staircase

Let

```text
P=F union LA union {z_1}
```

be the twelve-speed family obtained by deleting `z_2`, and write

```text
M(P)=max_t min_(v in P)||vt||.
```

There are two branches.

If `M(P)>=2/25`, the same Lipschitz argument as `(37f)` gives a closed arc
safe for `P` at level `1/14` of length at least

```text
2(2/25-1/14)/(La_*)=3/(175La_*).
```

It must lie strictly inside one open `z_2`-tooth.  Therefore

```text
z_2/L=c_2/d_2<(25/3)a_*.                              (37ib)
```

Suppose instead that `M(P)<2/25`.  Order

```text
a_1<a_2<a_3<a_4<a_5
```

and put

```text
p=4/25,          eta=p(1-p)=84/625,
C_0=585/154,
C_s=max(C_0,a_s)                 (1<=s<=4).            (37ic)
```

For `0<=s<=4`, let

```text
P_s=F union {z_1} union {La_1,...,La_s}.
```

This has `7+s` speeds and maximum at most `LC_s`.  Cited LRC through
thirteen supplies a point with clearance at least `1/(8+s)`, and
Lipschitz fattening supplies a closed `2/25`-safe arc of length at least

```text
ell_s=
 2(1/(8+s)-2/25)/(LC_s).                              (37id)
```

The remaining `r=5-s` aligned combs must cover that arc.  The exact
interval discrepancy at radius `2/25` is

```text
mu(I intersect D_w^(2/25))
 <=p|I|+eta/w.                                        (37ie)
```

Indeed, the centered period-one indicator of an interval of length `p`
has a primitive of range exactly `p(1-p)=eta`; scaling by `w` proves
`(37ie)`.  Summing it over the remaining aligned speeds and using
`(37id)` gives the reciprocal invoice

```text
sum_(i=s+1)^5 1/a_i >= k_s/C_s,

k_s=
 2(1/(8+s)-2/25)(1-(5-s)p)/eta.                       (37if)
```

The five exact tariffs are

```text
s             0        1        2        3        4
k_s         15/112    1/6     13/84    17/154    1/24
(5-s)/k_s  112/3      24      252/13   308/17     24. (37ig)
```

Since the remaining reciprocals are each at most `1/a_(s+1)`,
`(37if)` recursively yields

```text
a_1 <=       141,
a_2 <=     3,384,
a_3 <=    65,597,
a_4 <= 1,188,463,
a_5 <=28,523,112.                                     (37ih)
```

The harmonic inequalities `(37if)`, rather than only their rectangular
consequence `(37ih)`, are the sharper finite search object.  Together with
`L|5,045,040`, `(37a)`, and `(37g)`, they make the entire
`M(P)<2/25` subbranch body-uniformly finite.  No emptiness census of that
large box is claimed.  The unbounded branch is confined instead by the
strict projective cone `(37ib)`.

### The body projection and relaxed address-mask sieve

The quotient has a second exact projection that forgets the aligned shape
but retains the complete body address support.  Put

```text
Y_D(A)=union_(r in S_D)(r+R_A)/D.                    (37ii)
```

The copies in `(37ii)` lie in distinct `1/D` cells.  Consequently

```text
mu(Y_D(A))=(|S_D|/D)u_A.                             (37iii)
```

If `(37k)` holds, then

```text
Y_D(A) subset D_(a_1) union D_(a_2).
```

The two quotient speeds are distinct.  The repaired universal pair floor
from THM-1166 (via LEM-042) therefore gives

```text
mu(D_(a_1) union D_(a_2))
 =2/7-rho(a_1,a_2)
 <=2/7-1/91
 =25/91.
```

Together with `u_A>=478/1365`, this proves the support obstruction

```text
|S_D|/D<=375/478.                                    (37iv)
```

The exact all-body census ranges over

```text
F in C({1,...,14},6),             D divides L_F.
```

There are `251,536` such body/divisor rows.  Condition `(37iv)` eliminates
`240,560`, leaving `10,976`.  Every surviving divisor satisfies `D>=42`,
and each body has at most six surviving divisors.  The unique row at the
minimum divisor is

```text
F=(1,2,3,4,6,12),       L_F=168,       D=42,
|S_D|/D=32/42=16/21.                                 (37v)
```

There is a numerator-free denominator screen on every remaining row.  Fix
one `u in R_A` and define

```text
M_i(u)={r in Z/DZ: ||c_i(r+u)/d_i||<1/14}.
```

Multiplication by `c_i` permutes `Z/d_i Z`.  An open arc of length `1/7`
contains at most `ceil(d_i/7)` points of that equally spaced grid, and each
such residue has `D/d_i` lifts.  Hence

```text
|M_i(u)|<=(D/d_i)ceil(d_i/7),

|S_D|/D
 <=ceil(d_1/7)/d_1+ceil(d_2/7)/d_2.                 (37vi)
```

Among the `3,066,274` divisor pairs with

```text
d_1,d_2>1,                  lcm(d_1,d_2)=D,
```

on the support-hard rows, `(37vi)` leaves exactly `23,755` pair
occurrences on `6,292` rows.

Reflection supplies a useful visible subcase.  The involution

```text
sigma(r)=D-1-r
```

preserves `S_D`.  When `d_1=2`, it swaps the two parity classes, so each
contains exactly half of `S_D`; the first clock can meet at most one class.
The other denominator `d_2` must therefore satisfy

```text
|S_D|<=2(D/d_2)ceil(d_2/7).                          (37vii)
```

This kills `6,754` of the `6,756` cardinality survivors with denominator
two.  The two remaining parity rows also fail the stronger address test
below.

For a general denominator `d`, the actual fixed-`u` mask is the lift of at
most `ceil(d/7)` consecutive residues after multiplication by a unit modulo
`d`.  Enlarge it to exactly that many residues and allow the two clocks to
choose their phases independently.  This is an upper relaxation of the
literal mask problem.  For every pair with

```text
min(d_1,d_2)<=1000,
```

the exact verifier maximizes the first enlarged arithmetic-progression load
on `S_D` and grants the second its full cardinality capacity.  That already
kills `22,811` of the `22,813` pairs.  On each of the two survivors, for
both possible first parity masks, project the residual to
`P subset Z/d_2 Z`.  Containment in a `k_2=ceil(d_2/7)` term cyclic
arithmetic progression would necessarily imply

```text
|P|<=k_2,                    |P-P|<=2k_2-1.           (37viii)
```

Both rows violate the difference-set inequality.  Thus no denominator pair
with a denominator at most `1,000` survives even this relaxed necessary
test.  The entire remaining denominator ledger consists of `942`
occurrences, both denominators greater than `1,000`, on exactly two rows:

```text
F                         D=L_F     |S_D|/D       pairs
(1,4,5,7,9,11)           194040     2308/8085      371
(1,5,7,8,9,11)           388080     3029/10780     571. (37ix)
```

There is also an actual-slope consequence near the top of the support
window.  If

```text
x=|S_D|/D>2535/3346,
delta=x(478/1365)-13/49>0,
```

then a cover would force

```text
rho(a_1,a_2)<=1/49-delta.
```

Write `a_i=g alpha_i` with coprime `alpha_1<alpha_2`.  LEM-043 gives

```text
rho(a_1,a_2)>=1/49-1/(7alpha_2),
alpha_2<=1/(7delta).                                 (37ixa)
```

Moreover `g<=a_1<(585/154)D`, so `(37ixa)` bounds the actual second slope
`a_2=g alpha_2`, not only its projective ratio.  The exact census places
`1,150` quotient rows in this body-uniform finite-slope sector; its largest
displayed integer cap is

```text
a_2<=26,927,449,547.
```

The separate identity

```text
gcd(g,D)
 =gcd(D/d_1,D/d_2)
 =D/lcm(d_1,d_2)
 =1
```

is a useful clock sidecar but is not itself the size bound.

For this symmetric denominator ledger, reorder the two clocks so that
`d_1<=d_2`; this relabelling is independent of the earlier ordering by
speed.  The remaining `942` pairs admit a still coarser, and therefore
cheaper, completion.  Put

```text
k_d=ceil(d/7),
lambda_d(s)=#{r in S_D:r=s mod d},
Lambda_d=sum of the k_d largest values of lambda_d.  (37ixb)
```

The first actual mask is contained in the lift of some `k_(d_1)` residue
classes, even if all arithmetic-progression structure is forgotten.
Therefore it meets at most `Lambda_(d_1)` points of `S_D`.  Granting the
second mask its full ambient capacity gives the necessary inequality

```text
|S_D|<=Lambda_(d_1)+(D/d_2)k_(d_2).                  (37ixc)
```

The exact load ledger violates `(37ixc)` for `940` of the `942` pairs.  Its
only survivors are the two diagonal pairs

```text
F                         D          (d_1,d_2)       two capacities
(1,4,5,7,9,11)           194040     (D,D)           27720+27720
(1,5,7,8,9,11)           388080     (D,D)           55440+55440. (37ixd)
```

They fail for a structural reason hidden by cardinality.  Write `D=7k`.
Every enlarged diagonal mask has the form

```text
B(s,h)={s+jh mod D:0<=j<k},             gcd(h,D)=1.  (37ixe)
```

Because `h` is also a unit modulo `k`, projection of `(37ixe)` to
`Z/kZ` is a bijection.  Thus one mask contains exactly one point in each
modulo-`k` fibre, and two masks contain at most two.  The first row in
`(37ixd)` has four support points

```text
29701, 57421, 112861, 168301 = 1981 mod 27720,
```

while the second has three

```text
59401, 114841, 225721 = 3961 mod 55440.
```

So neither support word can be covered by two diagonal masks.  As a hostile
control, the complete fibre-multiplicity histograms are

```text
row 1: N_0=3960, N_1=0, N_2=17472, N_3=4704, N_4=1584;
row 2: N_0=7920, N_1=0, N_2=33516, N_3=14004.       (37ixf)
```

This finishes the necessary-mask ledger: all `23,755` denominator-capacity
survivors are impossible.

The two support implementations (cyclic bitsets and merged integer arcs),
the denominator and parity ledgers, and every mask relaxation are replayed
with checks active under optimized Python.  Ordinary and optimized outputs
are byte-identical.  The frozen hashes are

```text
support script  778842c0e8e7172835ca6ae673fb6156f212d4296e672bce4e7cc2815195bf1a
support output  648327d3b9b5b9a50c7760f0afd89a7a33161f57fa98c1b9e181d6b5b791a25f
address script  870498c4f0a2d97a2d42bce593c44283c77a141fb08669b4a91133e39db5c276
address output  74f7c270034dc40b4de3d33b9abf67481435d1e97eb6e52f2448d0a152cb68d7
completion script c0a07747c300c7e82d3da27b4f498425ffb72dc950f797381d1f0d7e9096655c
completion output 3ff1c54a2818b9a0f061912758584200ffa4dd5b549ea91ba2c6c7e6a92f3638.
```

To reconnect the relaxed ledger to the exact pointwise cover, put

```text
D=lcm(d_1,d_2),          S_D=J mod D,
a_q=c_qD/d_q,            X_r=(r+R_A)/D.               (37j)
```

Then `D` divides `L`, the `a_q` are integers, and the original residual is
covered pointwise if and only if

```text
X_r subset D_(a_1) union D_(a_2)       for every r in S_D. (37k)
```

This is a finite rational interval-arrangement test.  Isolated points pay
zero in `(37h)` but must remain in the final pointwise check `(37k)`.  If
`(37k)` held, then for any fixed `u in R_A` its values at
`(r+u)/D`, `r in S_D`, would give

```text
S_D subset M_1(u) union M_2(u).
```

The carrier `R_A` has positive mass, while the completed necessary-mask
ledger proves that no such inclusion exists for any admissible denominator
pair.  This contradiction proves:

> **Literal two-drift closure.**  Five aligned tail combs and two arbitrary
> tail combs cannot cover the safe carrier of six distinct literal body
> combs.

For reuse in higher-drift branches, the exact pair sidecar is
carrier-local.  If
`a_1=g alpha`, `a_2=g beta`, with `gcd(alpha,beta)=1`, and `P_(alpha,beta)`
is a periodic primitive of

```text
1_(D_alpha)1_(D_beta)-rho(alpha,beta),
```

then for interval components `[L_s,R_s]` of `X_r`,

```text
mu(X_r intersect D_(a_1) intersect D_(a_2))
 =mu(X_r)rho(alpha,beta)
  +(1/g)sum_s[P_(alpha,beta)(gR_s)-P_(alpha,beta)(gL_s)]. (37l)
```

The endpoint term depends on the primitive pair and the located carrier.
It cannot be replaced by the global overlap `rho(a_1,a_2)`: THM-1166's
hostile interval `[1/7,6/7]` has zero local overlap for every consecutive
pair `(N,N+1)`, although its global overlap is positive.

The local-current identity remains valid, but the two-drift branch no longer
needs its growing numerator box: the coarser address quotient is already
empty uniformly over all literal bodies and aligned shapes.

### The all-`k` support-transfer ladder

The support projection used above does not belong only to the two- and
three-drift sectors.  Let `k` tails be aligned, let `p=7-k` be genuine
drifts, and write their reduced clocks as

```text
z_i/L=c_i/d_i,       gcd(c_i,d_i)=1,       d_i>1,
D=lcm(d_1,...,d_p),       a_i=c_iD/d_i.
```

Then `D|L`, and the quotient speeds `a_i=Dz_i/L` are distinct.  If `R_A`
is the common safe set of the `k` aligned multipliers, the body-address
projection is

```text
Y_D(A)=union_(r in S_D)(r+R_A)/D,
mu(Y_D(A))=(|S_D|/D)mu(R_A).                         (37l+)
```

The pieces in `(37l+)` occupy distinct `1/D`-cells up to null endpoints.
A cover forces `Y_D(A)` into the union of the `p` quotient danger combs.
Writing `u_j` for the proved uniform safe-mass floor for `j` distinct
combs gives

```text
mu(R_A)>=u_k,       mu(union_i D_(a_i))<=1-u_p,
|S_D|/D <=(1-u_p)/u_k.                               (37l++)
```

The inequality is non-strict.  For `k=0`, take `R_empty=T` and `u_0=1`;
the `p=0` all-aligned branch is separate.  With

```text
(u_0,...,u_7)
=(1,6/7,66/91,55/91,558/1183,478/1365,61/273,15/154),
```

the exact necessary cutoffs and literal body/divisor census are

```text
k                0        1        2        3       4       5      6   7
p                7        6        5        4       3       2      1   0
cutoff       139/154  106/117  887/990  125/143   26/31  375/478 39/61 0
rows          27,210   27,240   27,163   26,970  13,778   10,976  6,237 0.
                                                               (37l+++)
```

These filters are not nested: for example `106/117>139/154`.  They are
necessary quotient-address rows, not realized covers.

There is also an exact count without enumerating denominator tuples.  For
`p>=1`, let `a_p(D)` count nondecreasing `p`-tuples of divisors `d_i>1`
with lcm exactly `D`.  The number of such multisets whose entries merely
divide `e` is

```text
B_p(e)=binom(tau(e)+p-2,p)=sum_(f|e)a_p(f).
```

Divisor-poset Möbius inversion gives

```text
a_p(D)=sum_(e|D) mu(D/e) binom(tau(e)+p-2,p).        (37lM)
```

The resulting denominator-shape and body/divisor occurrence counts are

```text
k       denominator shapes               raw occurrences
0       161,535,777,082,757               1,504,842,061,942,849
1         3,095,010,121,875                  38,954,725,590,760
2            50,874,159,718                     951,545,890,235
3               694,921,995                      21,357,714,101
4                 7,483,350                         298,255,882
5                    56,419                           3,066,274
6                       171                               6,237
7                         0                                   0.
```

Here a shape is an unordered denominator multiset, and a raw occurrence is
a `(body,D,multiset)` triple.  Numerators, directions, phases, physical
speed packets, and cover realizability are not counted.  An independent
referee reconstructs all `3,003` bodies and `251,536` body/divisor rows by
an endpoint sweep and merged cyclic arcs, and recomputes `(37lM)` by
downward divisor recurrence.

### Septimal peel and exact floor/exception compression

The enormous `k=2,p=5` column has a distinguished-prime structure that the
raw multiset count hides.  First, every support-hard resolving modulus is
divisible by seven.  Indeed, put

```text
L=14 lcm(F),          M=L/D,
```

and suppose `7` does not divide `D`.  In one residue fibre
`j=r+Dt`, the cell word of a body speed `f` has period `P=L/f`.
The quotient period

```text
P/gcd(P,D)
```

is divisible by seven, and the half-open danger block has length `P/7`.
Consequently that speed removes exactly `M/7` cells from every `D`-fibre.
Six body speeds remove at most `6M/7`, so every fibre retains a safe cell:

```text
7 does not divide D  ==>  S_D=Z/DZ.                 (37lS1)
```

Because the support-transfer cutoff `887/990` is strictly below one, no
nonseptimal row is support-hard.  The exact implementation independently
checks all `96,235/96,235` nonseptimal body/divisor rows.

Now write `D=7q`.  For every denominator `d|D`,

```text
d/gcd(d,q) in {1,7}.                                (37lS2)
```

If `d` does not divide `q`, the enlarged trace is a transverse section: it
contributes at most one point to every `q`-fibre.  If `d|q`, it is vertical.
Writing

```text
d=7a+r,       w=q/d,       0<=r<7,
```

the number of hit fibres at common phase `u` has the pointwise upper
envelope

```text
Y_d(u)<=Ybar_d(u)=w(a+X_d(u)),    X_d in {0,1},
integral X_d=r/7.                                      (37lS3)
```

Equality `Y_d=Ybar_d` holds almost everywhere; when `r>0`, `X_d` can be
chosen for pointwise equality as well.  When `r=0`, a strict-mask equality
phase can lower `Y_d` from `wa` to `w(a-1)`, but never raises it above the
granted floor `wa`.  In particular `integral Y_d=q/7`, independently of
`d` and of the reduced numerator.  This is an exact seven-sheeted Fubini
law, not a probabilistic independence assumption.  Every sieve below uses
`Ybar_d`, so the seam defect only enlarges the necessary relaxation.

Suppose `c` of the five denominators are transverse and let

```text
N_c=#{b mod q: lambda_q(S_D,b)>c}.
```

For the remaining vertical denominators, sum their deterministic floors
`wa` and retain every nonconstant bit `X_d` with its weight `w` and marginal
`r/7`.  Coverage forces the corresponding open upward weighted event
throughout the compact aligned-safe carrier.  If the event is the whole
circle, its mass is already strictly above `u_2=66/91`; otherwise compact
containment in an open proper subset makes the comparison strict.  The
fractional-cover theorem above gives its exact maximum over all joint laws
with these marginals.  Its coarser expectation consequence is

```text
66 N_c < 13(5-c)q,                                 (37lS4)
```

with `N_c=0` handled separately; the implemented floor/exception test keeps
the strictly stronger exact upward-event bound.

For a positive-remainder vertical multiset `P`, zero-remainder multiplicity
`z`, and transverse multiplicity `c`, exact-lcm completion still avoids
five-tuple enumeration.  If `ell=lcm(P)`, `U_D(E)` counts transverse
divisors in the downset of `E|D`, and `Z_D(E)` counts its zero-remainder
vertical divisors, then the completion coefficient is

```text
sum_(ell|E|D) mu(D/E)
  Mult(U_D(E),c) Mult(Z_D(E),z).                    (37lS5)
```

The exact ledgers are

```text
stage                         shapes          occurrences       rows
raw support transfer      50,874,159,718   951,545,890,235     27,163
expectation screen        36,962,285,549   320,011,786,356      4,592
floor/exception screen    26,908,162,790   200,389,247,292      4,414.
                                                               (37lS6)
```

The last row occupies `2,977` bodies and `149` resolving moduli.  By number
of transverse sections `c`, its occurrence counts are

```text
c=1:      3,524,756       c=2: 46,320,209,782
c=3: 112,812,921,408      c=4: 39,762,283,189
c=5:  1,490,308,157.                                (37lS7)
```

Thus the exact screen removes `78.94%` of the raw occurrence ledger, but it
does not close `k=2`.  Before located-phase refinement its minimum survivor
is `D=168`, `F=(1,2,3,4,6,12)`, with four transverse sections and one
denominator-six vertical mask.  The relaxation still forgets the actual
common-phase locations of the optional bits, their positions inside the
seven-sheeted fibres, compatibility between quotient scales, and the
reduced numerator/unit direction.  Those are precisely the coordinates
retained by the residue-ray/common-status sidecar of THM-2941.

The Python referee reconstructs all body rows, checks `52,925` literal
exact-lcm shapes through `D<=300` and `1,880` suffix queries, and requires
byte-identical ordinary/optimized output.  Its exact C++ engine is required
to agree under `-O2/-O3`.  Canonical source/engine/output hashes are

```text
085b4e2747a48bdbc1125e894af7d4f647dfdd7be86a00cf02dea2a8667e26dc
664d0df36d104d959279605c8ea8539d61ab595b155e5157fa7d0433f1b7944c
f711376bdfa0064f70d76e42505cf1eb89dadc4c66bcacd90d985b6641c2cd75,
```

and the survivor semantic digest is
`2eec9a97f02a7b8f8e36e50f747d53186ff5b84234a14fc3ac64818b54033675`.

### Located phase closes the sole-denominator-six family

The first survivor exposes exactly which coordinate the aggregate optional
bit erased.  In the `c=4` slice, suppose the sole vertical denominator is
`d=6` and `N_4>0`.  Its reduced numerator is `s in {1,5}`.  The spike-high
event and its closed low set are

```text
V_s={u:||su||<3/7},
C_s={u:||su-1/2||<=1/14},       mu(C_s)=1/7.          (37lS11)
```

Coverage would force the two-aligned safe set `R_(a,b)` into `V_s`.  This
is impossible uniformly in the distinct aligned multipliers `a,b`.  Put
`g=gcd(a,s)`, `A=a/g`, `S=s/g`, and
`B_2(x)=x^2-x+1/6` periodically.  Exact Fourier overlap, with null
endpoints immaterial, gives

```text
mu(C_s intersect D_a)
 =1/49+
  [B_2({(S+6A)/14})-B_2({(S+8A)/14})]/(AS)
 =1/49+h/(49AS).                                     (37lS12)
```

For `S=1`, the numerator `h` by `A mod 14` is

```text
(0,-1,5,-3,3,-5,1,0,-1,5,-3,3,-5,1),
```

and for `S=5` it is

```text
(0,-5,4,-8,8,-4,5,0,-5,4,-8,8,-4,5).
```

These finite residue tables prove

```text
mu(C_s intersect D_a)<=1/14,  equality only when a=2s;
mu(C_s intersect D_a)<=1/28,  otherwise.              (37lS13)
```

Since `a` and `b` are distinct,

```text
mu(C_s intersect (D_a union D_b))<=3/28<1/7=mu(C_s).
```

Hence `C_s` contains a point of `R_(a,b)`, contradicting
`R_(a,b) subset V_s`.  The entire sole-`d=6`, `c=4`, `N_4>0` subfamily is
therefore empty without a multiplier horizon.

Divisor-Möbius completion and a residual-union pass remove exactly
`248,154,771` occurrences and `8,998,004` globally surviving shapes.
The exact floor/exception ledger becomes

```text
shapes             26,899,164,786
occurrences       200,141,092,521
rows                        4,354
bodies                      2,966
resolving moduli              147.                    (37lS14)
```

Sixty rows, eleven bodies, and two resolving moduli disappear after union
over all `c` and shapes; the minimum survivor moves to `D=336`,
`F=(1,2,3,4,8,12)`.  The referee checks `1,000` literal overlaps, `18`
small exact-lcm moduli containing `1,029` fixed-`d=6` shapes, the full
residual union, and ordinary/optimized byte identity.  Canonical
source/output hashes are

```text
9f300459b273ad1825d3fe3e9274c6afe609f2d581e9df3d2be1780d347e541b
187e4f0e48eb2c93bdeac083f09707739b86656f5e88f9e4cbdb6d5019fbd0f7.
```

This located lemma closes one exact family only.  It does not license
subtracting all `18,064,772` fixed-`d=6` shapes before the residual-union
pass, and it does not close the remaining `k=2` quotient.

### Four-drift septimal mean and exact one-spike quotient

The same seven-sheeted object gives a much smaller first exact quotient for
`k=3,p=4`.  Every one of its `26,970` support-transfer rows is septimal, so
write `D=7q`.  Grant each of `c` transverse denominators one point in every
`q`-fibre.  For the target support define

```text
N_c=#{b mod q:lambda_q(S_D,b)>c}.
```

The remaining `m=4-c` vertical masks fill whole fibres.  If `T(u)` is their
total number of filled fibres, seven-sheeted Fubini gives

```text
integral_T T(u)du=mq/7.
```

A cover would force the open event `{T>=N_c}` to contain the compact
three-aligned safe set, whose mass is at least `u_3=55/91`.  If this event
is the whole circle, then `N_c<=integral T=mq/7`, so
`55N_c<=55mq/7<13mq`.  Otherwise

```text
55/91<=mu(R)<mu({T>=N_c})<=mq/(7N_c).
```

Thus Markov gives, for `N_c>0`, the strict necessary condition

```text
55N_c<13(4-c)q,                                      (37lS8)
```

while `N_c=0` is retained separately.  No independence among the vertical
bits is assumed.

When `c=3`, there is only one vertical denominator.  Write

```text
d=7a+r,       w=q/d,       0<=r<7.
```

Its a.e. fibre count is `aw+wX`, with `integral X=r/7`, and this is a
pointwise upper envelope at strict seams.  The superlevel event can have
mass strictly larger than `55/91` exactly only when

```text
N_3<=floor(d/7)(q/d)+1_(r>=5)(q/d).                  (37lS9)
```

Divisor-Möbius inversion is performed in grouped feature space, retaining
the unique `c=3` spike denominator rather than enumerating four-tuples.
Intersecting these screens with the separately necessary support/status
screen gives

```text
stage                                      shapes       occurrences       rows
raw support transfer                 694,921,995    21,357,714,101     26,970
support/status                       694,254,050    13,280,722,299     18,599
mean screen AND support/status       400,005,870     2,934,202,044      2,120
hybrid: exact c=3, mean c!=3          398,241,574     2,548,901,482      1,904.
                                                                  (37lS10)
```

The last row occupies `1,823` bodies and `107` resolving moduli.  Its
`c=1,2,3,4` occurrence counts are respectively
`71,619,386`, `1,351,841,956`, `1,065,317,472`, and `60,122,668`: only
the `c=3` slice is replaced by `(37lS9)`, while the other slices retain the
mean screen `(37lS8)`.  Its
minimum survivor is `D=840`, `F=(1,3,4,5,6,10)`; for the unique-spike
slice the first survivor has `d=5`, remainder `5`, and allowance `24`.
The mean equality locus is empty, and the literal phase audit checks `485`
individual masks, `7,904` strict boundary cells, `441` one-spike
thresholds, and `626` four-mask Markov thresholds.

This is a necessary-state compression, not a closure of `k=3`.  The
support/status and spike predicates are intersected only as separate
necessary conditions: their shared phase bits, actual transverse-section
locations, unit directions, and cross-scale compatibility are not yet
optimized jointly.  The abandoned general all-`c` feature expansion is not
a dependency.  Canonical base/q7/primary/audit source hashes are

```text
2fcd1fa7f122517feff3d3e0b3a21a6664fefaa12588e4db008572078989d6eb
3cc07195d580c5c5c01457ea95b58837a25c2176d326d12feaccc8e0bfa28dcc
05e365a654b32e66b814dcbce9385a2d13c22a2c84a5474e0855dcab6262b055
a3a65ae7d2fd05c7efc2b9ad5338eeda65fd5ddce6a5812af378e205aeae1065,
```

and primary/audit output hashes are respectively
`aa14b0894fb368df1723351d1c8fc92a5ba0e2c2cc3808b6285cade79ab468d1`
and
`eff21ad878828f002a38d4cc258cdb02a060c80ba943335bf0f397b6de773fff`.

### Finite-ring Kakeya needles and the first three-drift sector

The preceding transversal is one face of a general exact pushforward law.
This is a finite cyclic arithmetic-needle statement, not an invocation of
planar Kakeya dimension theory.  For `d,q|D`, a unit `h mod d`, and
`0<=ell<=d`, put

```text
B_D(d;a,h,ell)
 ={x in Z/DZ:x mod d in a+h{0,...,ell-1}}.           (37m)
```

Every actual fixed-phase denominator-`d` mask is contained in such a set
with

```text
ell=ceil(d/7).
```

Let

```text
g=gcd(d,q),       H=D/lcm(d,q),       ell=Ag+r,
0<=r<g.                                               (37n)
```

Then, for every `b mod q`,

```text
#{x in B_D:x=b mod q}
 =H(A+1_(b mod g in a+h{0,...,r-1})).                (37o)
```

Indeed, the congruences

```text
x=a+jh mod d,                    x=b mod q
```

are compatible exactly when `a+jh=b mod g`, and each compatible pair has
`H` solutions modulo `D`.  Since `h` is a unit modulo `g`, the relevant
class of `j mod g` occurs `A` or `A+1` times among `0,...,ell-1`.

Thus a needle pushes forward to a constant layer plus a shorter unit needle,
repeated through `Z/qZ`.  In particular its sharp phase-free fibre cap is

```text
(D/lcm(d,q))ceil(ell/g)
 <=(D/lcm(d,q))ceil(d/(7g)).                          (37p)
```

For several masks, the sum of `(37p)` is a necessary fibre capacity.  More
precisely, `(37o)` retains the locations of the exceptional fibres and may
be iterated: inside a fixed `q`-fibre, the trace is itself a lifted unit
needle on the smaller modulus `d/g`.  This Euclidean quotient/remainder law
is the exact toothpick self-similarity that the raw cardinality quotient
forgets.

At the first still-open mixed layer, take four aligned tail combs and three
drifts.  The proved four-comb floor and sharp three-comb union theorem give

```text
u_A>=558/1183,          mu(D_a1 union D_a2 union D_a3)<=36/91.
```

Consequently every cover must satisfy

```text
|S_D|/D<=(36/91)/(558/1183)=26/31.                   (37q)
```

If a drift denominator were one, that drift would be aligned and Section 9's
two-drift closure would apply.  Hence the genuine three-drift denominator
universe is

```text
2<=d_1<=d_2<=d_3,              lcm(d_1,d_2,d_3)=D.
```

The exact literal-body census for `(37q)` has

```text
251,536 body/divisor rows,
237,758 support kills,
13,778 support-hard rows on 206 divisors.             (37r)
```

Across those rows there are `7,483,350` denominator-triple shapes and
`298,255,882` row occurrences.  Write

```text
C_i=(D/d_i)ceil(d_i/7)
```

and let `Lambda_i` be the sum of the `ceil(d_i/7)` largest `d_i`-class
loads of `S_D`.  The four nested necessary relaxations

```text
C_1+C_2+C_3,
Lambda_1+C_2+C_3,
Lambda_1+Lambda_2+C_3,
Lambda_1+Lambda_2+Lambda_3
```

leave, respectively,

```text
143,852,683; 44,573,157; 1,385,991; 544,571          (37s)
```

occurrences.  The last ledger still meets `13,577` body/divisor rows, so
this is a substantial compression rather than three-drift closure.

The diagonal denominator sector *does* close.  For

```text
(d_1,d_2,d_3)=(D,D,D),              D=7k,
```

equation `(37o)` with `q=k` says that every enlarged mask is a section of

```text
Z/DZ -> Z/kZ.
```

Three masks therefore meet any fibre in at most three points.  Among the
`2,636` diagonal capacity survivors, the exact fibre histogram kills
`2,601`; their maximum multiplicities are

```text
maximum 3: 35 rows,   4: 702,   5: 403,   6: 1,496. (37t)
```

Every one of the `35` saturated rows has `D=L_F`, contains body speed `7`
but not `14`, and has nonempty first-level fibre sizes only `2` and `3`.
They expose the promised `chi_7`/Fano coordinate, but it is not yet the
terminal one.

In fact every saturated row has `D=49m`.  Enlarge a diagonal mask to its
full `7m`-term needle and fix a residue modulo `m`.  Its trace in that
fibre is a seven-point unit arithmetic progression in `Z/49Z`, by the
iterated form of `(37o)`.  Grant the three masks independent directions
and phases in every fibre; this is again an upper relaxation.  There are
exactly

```text
1,029
```

distinct seven-point unit needles in `Z/49Z`, and `147` pass through each
point.  An exact depth-three set-cover test finds, on every one of the
`35` rows, a body slice not coverable by three such needles.  The hostile
slice sizes range from `14` to `18`, strictly below the raw three-needle
capacity `21`.  Therefore:

> **Diagonal three-drift closure.**  Four aligned tail combs and three
> arbitrary drift combs cannot cover a literal six-body carrier when the
> three drift denominators all equal their common resolving denominator.

The same proof gives a genuine septimal ladder: if `7^s|D`, a full diagonal
`D/7`-needle restricts over `Z/(D/7^s)Z` to a
`7^(s-1)`-term unit needle in `Z/7^s Z`.  At the literal LRC(14) body scale,
the `Z/49Z` test is the terminal available septimal digit.  This was the
discovery proof of the diagonal sector.  The all-divisor argument below
subsumes its emptiness conclusion, but the `Z/49Z` proof remains an
independent exact check and records the septimal toothpick self-similarity
that the coarser transport relaxation forgets.

The combined referee and an independent diagonal-only implementation both
retain all checks under optimized Python; in each case ordinary and optimized
outputs are byte-identical.  Their frozen hashes are

```text
combined script   42dc165781148c702dfcd3c6535f4d02aee516af60b5ddf602a19cb1d87695e4
combined output   2e211620ad7064ea06f7544b5fbac709d6d52d9a0e261b464ae26b595f09b669
diagonal script   d887b0c0b202b6311e9b040693af1e96152f1a2eabeada245d6c02a75d80700c
diagonal output   478e3f2640a0484e9c21f1ca2ace8ca1adf6d3066e09c937ee9196bdfd8ccfa4.
```

### All-divisor overload and common-phase activity

The phase-free fibre cap `(37p)` is useful at every divisor, not only on
the septimal diagonal.  For `q|D` put

```text
g_i=gcd(d_i,q),       H_i=D/lcm(d_i,q),       ell_i=ceil(d_i/7).
```

Every cover necessarily satisfies

```text
max_(b mod q) #{x in S_D:x=b mod q}
 <=sum_(i=1)^3 H_i ceil(ell_i/g_i).                  (37ta)
```

Testing `(37ta)` for every divisor `q` kills `125,060` of the `544,571`
all-top occurrences and leaves `419,511`.  It kills all `2,636` diagonal
occurrences, so the earlier `Z/49Z` terminal is structurally informative
but no longer logically needed for diagonal emptiness.  The full
phase-free Lorenz profile strengthens `(37ta)` by comparing the sums of
the `s` largest fibres for every breakpoint `s`; it gives `147` additional
kills.  This independent strengthening is not needed by the terminal
chain below.

The decisive next coordinate is the common aligned phase.  For `d<7`,
`gcd(c,d)=1`, and `u in R_A`, the literal mask

```text
M_(c,d)(u)={r mod D: ||c(r+u)/d||<1/14}
```

meets at most one residue class modulo `d`, and

```text
M_(c,d)(u) is nonempty
 iff ||cu||<d/14.                                    (37tb)
```

The activity set on the right has Haar measure `d/7`.  If the other two
enlarged masks have total ambient capacity smaller than `|S_D|`, a cover
would force `(37tb)` for every `u in R_A`.  This is impossible for
`d=2,3`, because

```text
mu(R_A)=u_A>=558/1183>3/7.                           (37tc)
```

This common-phase activity test kills `383,391` of the `419,511`
single-fibre survivors.  There are `6,756` further survivors with
`d_1=d_2=2`.  Reflection `r -> D-1-r` preserves `S_D` and swaps parity,
so its two parity loads are equal.  In every one of these rows each parity
load exceeds the exact `q=2` fibre cap of the third mask.  Hence the two
denominator-two masks would have to be active, on opposite parities, for
every `u in R_A`.  Their individual activity sets have mass `2/7`, again
contradicting `(37tc)`.  The exact residual is therefore

```text
544,571 -> 419,511 -> 29,364.                        (37td)
```

The activity referee explicitly reconstructs the `q=2` histogram on every
even row; even total cardinality is not used as a substitute for reflection
balance.

### Eight-state threshold transport

The scalar and Lorenz bounds allow each needle to place its exceptional
fibres independently at every rank.  The next quotient retains the joint
status word.  Fix an outer divisor `t|D`, and write

```text
M=D/t,       g_i=gcd(d_i,t),       ell_i=A_i g_i+r_i,
0<=r_i<g_i,       R_i=(t/g_i)r_i.                    (37te)
```

By the exact trace law `(37o)`, needle `i` has load `H_i A_i` in every
outer fibre and one extra layer of height

```text
H_i=D/lcm(d_i,t)
```

on exactly `R_i` of the `t` fibres.  For a status word
`E subset {1,2,3}`, let `n_E` count outer fibres on which precisely the
needles in `E` receive their extra layer.  Every literal configuration
satisfies

```text
sum_E n_E=t,        sum_(E contains i)n_E=R_i,        n_E>=0.  (37tf)
```

We relax the integral, arithmetically located counts to real counts obeying
only `(37tf)`.  If

```text
c_E=min(M,sum_i H_i(A_i+1_(i in E)))
```

is the raw union-capacity upper bound of status `E`, and

```text
sigma_b=#{x in S_D:x=b mod t},
```

then every cover must satisfy, for every positive target load `y`,

```text
#{b:sigma_b>=y}
 <=max {sum_(E:c_E>=y)n_E : (37tf)}.                 (37tg)
```

This is the missing distributional statement: a needle may have enough
total mass while possessing too few exceptional fibres to serve all heavy
body fibres.  The polytope has eight variables and four equality
constraints.  The exact referee enumerates its `58` nonsingular four-column
bases out of `binom(8,4)=70`, using rational arithmetic, and independently
matches each decisive primal value with a dual value.

Applying `(37tg)` at every nontrivial divisor `t|D` (the case `t=1` is
vacuous) kills `29,345` of the `29,364`
activity survivors.  The `19` survivors lie on only five body words:

```text
F                              D          occurrences
(1,4,7,9,11,12)               38,808          2
(1,6,7,8,9,11)                77,616          4
(1,4,7,8,9,11)                77,616          2
(1,3,7,8,10,11)              129,360          1
(1,5,7,8,9,11)               388,080         10.     (37th)
```

There is a strictly stronger all-arity form which will be useful beyond
this closure.  For target-load multiplicities

```text
h_y=#{b:sigma_b=y},       H_y=sum_(z>=y)h_z,
```

a literal cover by any number of needles induces **one** status table
`n_E`, independent of `y`, such that

```text
sum_E n_E=t,       sum_(E contains i)n_E=R_i,
sum_(E:c_E>=y)n_E>=H_y for every y>0.                (37tg+)
```

Conversely, `(37tg+)` is exactly the feasibility condition for this
fractional capacity relaxation.  Introduce real variables `x_(y,E)>=0`;
their row sums are `h_y`, their bit marginals are `R_i`,

```text
sum_y x_(y,E)=n_E,
```

and `x_(y,E)=0` when `c_E<y`.  The nested-neighbourhood Hall criterion,
or greedy matching after sorting loads and capacities, eliminates `x` and
gives `(37tg+)`.  Thus the separate maxima in `(37tg)` are valid but may
choose incompatible status tables at different thresholds.

### Upward status transport is exactly fractional cover

The one-threshold optimization in `(37tg)` has a closed form for arbitrary
drift arity.  Let `A` be an upward-closed family of subsets of `[p]`, let
`r=(r_1,...,r_p) in [0,1]^p`, and define

```text
F_A(r)=max_mu mu(A),
```

where `mu` ranges over all laws on `2^[p]` with one-marginals `r_i`.  Let
`H=min(A)` be the clutter of inclusion-minimal true sets and put

```text
C(H)={w>=0:sum_(i in H)w_i>=1 for every H in H},
tau_H(r)=min_(w in C(H)) sum_i r_i w_i.
```

Then

```text
F_A(r)=min(1,tau_H(r)).                               (37tgC)
```

For the empty event take `tau=0`; for the full event, whose minimal clutter
contains the empty set, the cover polyhedron is infeasible and
`tau=+infinity`.

The upper bound is immediate: for every `w in C(H)` and every state `E`,

```text
1_A(E)<=sum_(i in E)w_i,
```

because a true state contains a minimal true set.  Taking expectations and
also using `mu(A)<=1` proves
`F_A(r)<=min(1,tau_H(r))`.

For attainment, first suppose `q=tau_H(r)<1`.  The dual fractional-packing
problem is

```text
max sum_(H in H)lambda_H
subject to sum_(H contains i)lambda_H<=r_i,       lambda_H>=0.
```

Choose optimal `w,lambda`, put mass `lambda_H` on each minimal true state
`H`, and write

```text
q=sum_H lambda_H,       a_i=sum_(H contains i)lambda_H,
s_i=r_i-a_i,            b=1-q,       Z={i:w_i=0}.
```

Complementary slackness gives `w_i s_i=0`, so all residual demand lies on
`Z`.  No true state is contained in `Z`, since such a state would contain a
minimal `H` of cover weight zero.

If `s_i>b`, move

```text
theta_i=s_i-b
```

of the current true-event mass from states lacking `i` to the corresponding
states with `i` added.  At the moment coordinate `i` is processed, the
available mass is `q-a_i`, and

```text
theta_i=r_i-a_i-1+q<=q-a_i
```

is exactly `r_i<=1`.  Splitting atoms if necessary, these moves preserve
the event by upward closure and leave residual marginals
`s'_i=min(s_i,b)`, still supported on `Z`.  Realize the remaining mass `b`
on `2^Z` by independent inclusions with probabilities `s'_i/b`.  Every
such state is outside `A`, so the resulting law has the prescribed
marginals and event mass `q`.

Now suppose `tau_H(r)>=1`.  Consider the upward blocking polyhedron

```text
Q=conv{1_H:H in H}+R_+^p.
```

If `r` were outside `Q`, strong separation would give `u` with

```text
u.r<inf_(x in Q)u.x=min_(H in H)u(H).
```

The recession cone forces `u>=0`; the strict inequality and `r>=0` force
the right side to be positive.  Normalizing `u` by that minimum would give
a fractional cover of cost below one, a contradiction.  Hence
`r>=sum_H alpha_H 1_H` coordinatewise for a probability vector `alpha`.
Start with mass `alpha_H` on the true states `H` and add each coordinate
`i` to mass `r_i-sum_(H contains i)alpha_H` of states currently lacking it.
There is enough such mass because `r_i<=1`, and upward closure again
preserves the event.  This constructs an all-true law with marginals `r`
and proves `(37tgC)`.

In the unnormalized `t`-fibre problem, with exact activity counts `R_i`,
`(37tgC)` becomes

```text
max_(real status tables) sum_(E in A)n_E
 =min(t,min_(w>=0, w(H)>=1)sum_i R_iw_i).             (37tgU)
```

Thus the `2^p`-cell real transport LP is exactly a fractional cover with
`p` variables and one constraint per minimal heavy status.  Integer status
tables are bounded above by the floor of `(37tgU)`; equality is not claimed
without an integrality argument.

This caveat is sharp.  Let `p=7`, let `H` be the seven lines of the Fano
plane, take the upward event of containing a line, and set `t=3`, `R_i=1`
for every point.  Summing the seven cover constraints, each point occurring
on three lines, gives `3 sum_i w_i>=7`; the uniform cover `w_i=1/3`
attains equality.  Hence `(37tgU)` has value `7/3`.  But an integer status
table assigns each point to exactly one of three states.  Two true states
would contain disjoint Fano lines, impossible because every two lines
meet.  One true state is attainable, so the integer optimum is `1`, not
even `floor(7/3)=2`.

In the LRC application
`A_y={E:c_E>=y}` is upward because the raw capacity upper bound `c_E` is
monotone under adding active needles.  The formula is exact for each fixed
threshold of the **real marginal relaxation**; it does not assert that the
raw capacities, integer tables, arithmetic locations, or literal needle
covers are realized.  The single-table, all-threshold condition `(37tg+)`
remains strictly stronger.

There is a useful equality upgrade in the common aligned phase.  Suppose the
status bits are open activity conditions on `u`, with the stated
one-marginals, and a cover would force the compact aligned-safe set `R_A`
inside a **proper open** upward event `V`.  The fractional-cover theorem is
an upper relaxation of the actual joint status law, while compactness and
proper openness make containment strict in measure.  Hence

```text
u_k<=mu(R_A)<mu(V)<=min(1,tau_H(r)).                  (37tgO)
```

Thus `tau_H(r)<=u_k` already kills the branch; equality does not survive.
For a single denominator `2<=d<7`, `(37tb)` gives activity marginal `d/7`.
The denominators ruled out from being forced active throughout `R_A` by
this one-marginal test are exactly

```text
k=1: {2,3,4,5,6},   k=2: {2,3,4,5},   k=3: {2,3,4},
k=4: {2,3},         k=5: {2},          k=6: empty.   (37tgA)
```

Larger denominators are merely not excluded by this test; no containing
phase is asserted to exist.
This alphabet explains the recursive migration seen in the exact sectors:
as `k` decreases, one new small denominator becomes activity-sensitive,
while survivors with none of these denominators must be attacked in
divisor fibres or by located arithmetic-progression shape.

This distinction is already strict in the four-needle scale relevant to
the next sector.  Take `D=28`, `t=2`, denominators `(4,4,28,28)`, and
widths `(1,1,4,4)`.  The only variable status bits have marginals `(1,1)`,
and their capacities are

```text
c_00=4,       c_10=c_01=11,       c_11=14.
```

A target with two fibre loads `(14,8)` passes every separately optimized
tail inequality and the maximum total-capacity test: disjoint high bits
give capacities `(11,11)` and coincident high bits give `(14,4)`.
No single table works, since the load-`8` tail forces the disjoint table
while the load-`14` tail forces the coincident one.  For four drifts the
joint relaxation has only `16` variables, five marginal equalities, and
at most fifteen consolidated tail inequalities; it should therefore be
tested before any local depth-four needle recursion.

### The `98/99` terminal local needles

The quotient trace retains one more self-similar coordinate.  In a fibre
of `x=b mod q`, of height `m=D/q`, put

```text
g_i=gcd(d_i,q),       n_i=d_i/g_i,
w_i=ceil(ell_i/g_i)=ceil(n_i/7).                   (37ti)
```

The trace of the `i`th enlarged mask is contained in the lift to
`Z/mZ` of a `w_i`-term unit arithmetic progression in `Z/n_iZ`.
The second equality in `(37ti)` is the elementary identity
`ceil(ceil(d/7)/g)=ceil((d/g)/7)` for `g|d`: quotienting a tooth does not
merely bound its width, but returns the same canonical width rule at the
reduced denominator.  This is the precise all-scale toothpick
self-similarity.
We grant the three masks independent phases and unit directions inside
the chosen fibre.  This is an upper relaxation: the actual global clocks
couple those choices rather than enlarging them.

For `18` rows in `(37th)`, take `m=99`, `q=D/99`, and `b=0`.  The only
trace-family types that occur are

```text
(m,n,w)       number of distinct lifted needles
(99,11,2)                     55
(99,33,5)                    330
(99,99,15)                 2,970.
```

For the remaining row take `m=98`, `q=1,320`, and `b=0`; its types are

```text
(98,49,7)                  1,029
(98,98,14)                 2,058.                    (37tj)
```

These family sizes are analytic, not unexplained census constants.  For
`2<=w<n-1`, a cyclic `w`-term unit progression in `Z/nZ` has exactly
`n phi(n)/2` distinct translates and directions: among the `n phi(n)`
pairs `(a,h)`, only reversal

```text
(a,h) -> (a+(w-1)h,-h)
```

identifies two descriptions.  Thus the five entries are respectively
`11*10/2`, `33*20/2`, `99*60/2`, `49*42/2`, and `98*42/2`.

The five distinct body slices have sizes `28,30,32,34,34`, and their
body/family data collapse to nine terminal profiles.  The four `m=99`
slices are invariant under `x -> -x`; only the `m=98` endpoint slice
breaks that exact reflection.  An exact
grouped set-cover recursion chooses the first uncovered target point and
tries every unused-family needle through it.  None of the `19` slices is
coverable by one needle from each of its three families.  A separate
literal-residue construction and incidence-bitset pair solver reproduces
all `19` failures; it does not use the bit-replication implementation of
the primary referee.

The audit also separates mechanism from search cost.  Five rows already
fail the sum of the three exact single-family target-intersection maxima;
ten more fail the exact union maximum of their repeated full-period pair
plus the third-family maximum.  Only four all-`(99,99,15)` rows remain.
Their exact three-needle target-coverage maxima are `27<32` on the
`(1,6,7,8,9,11)` slice and `24<28` on the
`(1,5,7,8,9,11)` slice.  Moreover, all `18` modulus-`99` targets omit
row zero modulo `9` and column zero modulo `11`, reflecting the common
body speeds `9` and `11`.  More generally, if `M=ab`, `gcd(a,b)=1`, and
a target avoids one row modulo `a` and one column modulo `b`, then a
`w<M` unit progression meets it in at most

```text
w-floor(w/a)-floor(w/b)+1:
```

the forced row and column hits can overlap only at their one CRT
intersection.  The terminal is therefore a finite-ring Kakeya orbit-cover
gap on a CRT grid, not a knife-edge phase search.

Therefore:

> **Uniform three-drift closure.**  Four aligned tail combs and three
> arbitrary drift combs cannot cover the safe carrier of six distinct
> literal body combs.

The proof is numerator-free because all unit directions and phases were
granted at the terminal fibre.  Its preserved coordinate is not merely
fibre mass but the finite unit-needle shape.  The two terminal moduli show
that the septimal `Z/49Z` stalk is one branch of a broader composite
toothpick recursion: the final obstruction may live on `2*7^2` or on
`3^2*11`, as selected by the body support.

The primary activity and terminal referees retain all guards under
optimized Python, and their ordinary and optimized outputs are
byte-identical.  The frozen hashes are

```text
activity script   067424a0edb126ad8e15f9f56ad77489bfdbb9a5fd2dd42e8a3eb4599438dfd9
activity output   cab8e7b4a63177e89f19979c3c8129360dc486786b5a594bc871b25c4bf5c723
Lorenz script     8d9ef3a750b7d09eb247090c04d3825b48be355c5fc653e57935fea7d98f49a5
Lorenz output     b1b5801c53411e3adf661e52765e040db29769076ab02b67168b5e9b087239c1
terminal script   13e524e728736480798d52acc736afe0cbd7b651a487aeef5d771d2b5dfa1338
terminal output   435c34b249255c8659e62e89e3a22e42b2fa8dd67c4cbd4b188d5e9806206f52
audit script      4b2f85d59f018bda6d75f52cc399ec6e3282b3015d8dba748eac8f41208c2019
audit output      0cd0ebb33a9cba143abf64ba9f5f18aa5da5dc63a5f74435ffb69698cf863af7
support script    8db781fb3e7dc8fdc4df2bf3c6d83869a9ffe52f41c7d70c25bbd0a9b0122bea
support output    808ec922a881e1d6d9541539ee51ff520e44c4fa7c98208b315c82f91df59e81
support audit     417830cff7a767227d93bbcee42ad57b75adf2b335dc5fa8fe50e85a972bb792
audit output      43a73d9daa2beafb69541db6dc9bf205d9f4d9e4ac0ebf83a268c868551923f4
cover audit       708465d1f6b154a4fb8d477bd793a2dcd7f8f652716fc1f9d8baf86d4d0ed9db
cover output      0eb35e2f990ecca043a18c08987bab92fed42968265bbb8ca6f7144f3787115c.
```

## 10. The critical affine-profile residual is empty

THM-2184 leaves one honest critical possibility: a fixed seven-tail
two-torus profile could in principle have zero safe mass.  The
divisor-minimal Fourier coordinate rules it out on every positive grid
carrier.

Let `C` be a nonempty union of `1/L` cells and fix

```text
c_i in Z_(>0),             r_i in Z,                 1<=i<=7,
W_i(N)=NLc_i+r_i.                                      (39)
```

For all sufficiently large `N`, the displayed speeds are positive; an LRC
packet also requires them to be distinct.  Put

```text
d(y)=1_(||y||<1/14),

Phi_(c,r)(t)
 =integral_T product_(i=1)^7(1-d(c_i x+r_i t)) dx,

P_(C;c,r)=integral_C Phi_(c,r)(t)dt.                  (40)
```

THM-2184's grid transfer gives

```text
|mu(C intersection intersection_i D_(W_i(N))^c)
       -P_(C;c,r)|
 <=5||r||_1/(2NL).                                    (41)
```

The key point is that its integral on every positive carrier is always
positive, even though the pointwise profile can vanish at isolated slow
phases.  Fix `t` and let

```text
m_t(x)=sum_(i=1)^7 d(c_i x+r_i t).
```

Every summand has mass `1/7`, so `integral m_t=1`.  The negative part of
`m_t-1` is precisely the uncovered set and has integral `Phi_(c,r)(t)`;
the mean-zero identity makes the positive part have the same integral.
Therefore

```text
||m_t-1||_(L1)=2Phi_(c,r)(t).                         (42)
```

Let

```text
c_0=min_i c_i,              I_0={i:c_i=c_0}.
```

At `x`-Fourier frequency `c_0`, no comb with larger slope contributes.
Since

```text
hat d(1)=sin(pi/7)/pi>0,
```

equation `(42)` gives the quantitative divisor-minimal obstruction

```text
Phi_(c,r)(t)
 >=sin(pi/7)/(2pi)
   |sum_(i in I_0) exp(2pi i r_i t)|.                 (43)
```

If `P_(C;c,r)=0`, then `(43)` makes the exponential polynomial on the
right vanish almost everywhere on one open carrier cell.  It must therefore
vanish identically.  After equal exponents `r_i` are grouped, however, all
its coefficients are positive integers.  It is not the zero polynomial.
This contradiction proves

```text
P_(C;c,r)>0                                           (44)
```

for every fixed positive slope vector and integer residue vector.  Equations
`(41)` and `(44)` give the effective terminal

```text
N>5||r||_1/(2L P_(C;c,r)).                            (45)
```

In particular, if the minimum slope `c_0` is unique, `(43)` has the
pointwise floor

```text
Phi_(c,r)(t)>=sin(pi/7)/(2pi).                        (46)
```

Thus a profile approaching zero must first repeat its divisor-minimal
frequency.  This is the shifted two-scale analogue of the
Mirsky--Newman divisor-minimal obstruction, now with the slow phase retained.

There is a rational quantitative form when all seven leading slopes agree.
Absorb their common value into `N`, and suppose the fixed residues
`r_1,...,r_7` are distinct.  If

```text
g_1(t),...,g_7(t)
```

are the cyclic gaps between the seven centers `-r_i t`, then

```text
Phi_(1,r)(t)
 =sum_j max(g_j(t)-1/7,0)
 =(1/2)sum_j |g_j(t)-1/7|.                            (47)
```

Consequently the profile vanishes exactly when the centers are a coset of
`(1/7)Z/Z`.  For any nonzero residue difference `d=r_i-r_j`, summing the
gap errors along the arc between those two centers gives

```text
Phi_(1,r)(t)
 >=(1/2)dist(dt,(1/7)Z/Z)
 =||7dt||/14.                                         (48)
```

On an interval of length `ell`, put

```text
n=7|d|,       n ell=q+s,       q=floor(n ell),       0<=s<1.
```

The exact translated-interval minimum is

```text
inf_(|I|=ell) integral_I ||nt||dt
 =(q+s^2)/(4n).                                       (49)
```

When `ell<=1/7`, the facts `n=7|d|` and `n ell<=|d|`
give from `(49)`

```text
integral_I Phi_(1,r)(t)dt>=ell^2/8.                  (50)
```

An independent exact census of all `3,003` literal six-body carriers finds

```text
min_F longest_component(G_F)=23/1092
```

at `F=(1,6,7,8,10,13)`.  Taking an interval of that length in `(50)` gives
the body-uniform profile floor

```text
P_(G_F;1,r)>=529/9539712.                             (51)
```

This yields a concrete uniform finite sector.  Let

```text
R_14=14 lcm(1,...,14)=5,045,040
```

and suppose seven distinct tails lie in one common quotient block,

```text
w_i=kR_14+b_i,             0<=b_i<R_14.
```

Here `sum_i b_i<7R_14`, so `(41)` and `(51)` give

```text
safe mass
 >529/9539712-35/(2k)>0
```

for every

```text
k>=315,586.                                            (52)
```

Thus the common-ruler block index is uniformly finite.  More generally,
`(44)` closes every fixed seven-tail affine ray, with its exact rational
profile providing the ray-specific terminal.

The dependency-free referee checks the equal-slope gap and coset laws,
`57,015` pair-distance inequalities, `12,480` exact interval formulas,
the all-body component minimum, the arithmetic in `(51)`--`(52)`, and
general minimum-frequency hostile controls.  Its checks remain active under
optimized Python; ordinary and optimized outputs are byte-identical.  The
frozen hashes are

```text
script  b12f0927312c1ea56b2e7fce1937b82cd49e7978c0d820c9c7621a9be838f13a
output  418e871e5a26b805ce3e161d85c5b500753e3a43090ae111a494caa4a9e93ddc.
```

This section empties the fixed affine-ray zero-profile residual, not the
union over changing slopes and residues.  A genuinely escaping packet must
therefore change its normalized affine data, or move between scale charts,
rather than travel along one fixed two-torus ray.

## 11. Scope and information audit

The proved transport is

```text
source:       labelled carrier cells and danger-tooth incidences;
map:          t=(j+u)/L;
preserved:    the complete Boolean word when every residue b is zero;
destroyed off the fixed locus:
              slope detuning b/L and cell phase bj/L;
necessary mixed-residue sidecar:
              J, bJ mod L, component centres/radii, and gcd fibers. (38)
```

For an arbitrary grid carrier, equations `(27)`--`(35)` remain necessary
rather than sufficient and need not eliminate every one-drift row.  The
double use of the six-comb floor and reflection is what closes the literal
six-body one-drift branch.  Section 9 now closes the literal six-body,
five-aligned/two-drift branch uniformly: its decisive quotient retains the
body address multiplicities but deliberately forgets the aligned shape,
clock phases, and numerator sizes.  Its finite-ring pushforward also closes
the entire four-aligned/three-drift branch: the threshold transport retains
the joint exceptional-fibre word, and the terminal `98/99` quotient retains
the local unit-needle shape.  The theorem does **not** close an arbitrary
mixed-residue seven-wall, branches with four or more drift speeds, or
LRC(14).  The next aligned sectors have two or three aligned combs and are
  finite by THM-2941, but their four/five-drift censuses have not been run;
  the zero/one-aligned sectors remain the not-yet-finitized address frontier.
  The support-transfer census shows why direct denominator enumeration is the
  wrong object: before numerators or phases it has `21,357,714,101`
  four-drift and `951,545,890,235` five-drift occurrences.  The septimal
  peel and fractional-blocker compression now reduce the latter to
  `200,389,247,292`; the sole-`d=6` located lemma reduces it further to
  `200,141,092,521`.  Mean and one-spike screens reduce the four-drift
  ledger to `2,548,901,482`.  Both residuals still deliberately forget
  most optional-bit locations.  The next exact object is the common
  seven-sheeted support
  pattern—base load together with located transverse sections and vertical
  exceptional fibres—joined to the THM-2941 unit-ray/common-status sidecar
  and, where necessary, a carrier-local multi-endpoint current.
A tournament on the seven labels forgets the metric widths, cell phases,
gcd fibres, quotient remainders, status multiplicities, minimal-heavy
clutter, and located endpoint current, and is therefore not an equivalent
quotient.
