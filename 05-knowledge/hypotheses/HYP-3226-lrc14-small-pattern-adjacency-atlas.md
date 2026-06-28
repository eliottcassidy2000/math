---
id: HYP-3226
title: LRC14 small-pattern adjacency atlas and motif payload ledger
status: EVIDENCE / executed atlas; not an LRC14 proof
source: codex-2026-06-28
tangent: T1324
technique: LTI-324
tournament_technique: LTT-224
script: 04-computation/lrc14_small_pattern_adjacency_atlas_codex_20260628.py
result: 05-knowledge/results/lrc14_small_pattern_adjacency_atlas_codex_20260628.out
reflection: 07-reflections/lrc14-small-pattern-adjacency-atlas-codex-20260628.md
related:
  - HYP-3235
  - HYP-3234
  - HYP-3233
  - HYP-3232
  - HYP-3231
  - HYP-3230
  - HYP-3218
  - HYP-3229
  - HYP-3228
  - HYP-3227
  - HYP-3217
  - HYP-3216
  - HYP-3215
  - HYP-3214
  - HYP-3225
  - HYP-3224
  - HYP-3223
  - HYP-3222
  - HYP-3221
  - HYP-3205
  - HYP-3204
  - HYP-3203
  - HYP-3202
  - HYP-3201
  - HYP-3200
  - HYP-3163
  - HYP-3154
  - HYP-3150
  - HYP-3138
  - OPEN-Q-108
---

# HYP-3226: LRC14 Small-Pattern Adjacency Atlas

## Claim After The Scout

Small patterns near the current LRC14 frontier are useful only when they are
typed payload atoms.  A motif may enter the proof workspace if it says which
LRC coordinate it preserves, which coordinate it destroys, and which sidecar
repairs the quotient loss.  It should not enter as raw numerology or as the
name of a famous problem.

The executed atlas supports this discipline.  The best motifs were not the
most exotic analogies; they were the motifs already attached to the current
normal-fan / trap-discharge payloads:

```text
comb-overlap Gram kernel
shell L_y magic quartic
normal-cone dual slack
multi-chart proof split
AP self-dual Fejer equidistribution certificate
three-gap Stern-Brocot cap-kernel recursion
consecutive plus doubled AP
modulus-covariance apex break
Toeplitz lambda-min margin
certificate-Helly separation
single-arc peeling recursion
ordered-tail exchange-rate ratio
D1/D2/D3 covariance layer split
totally-real cap field conductor packet
cyclotomic factor grading
```

## Method

The scout defines 94 motifs across 93 families, scans 8370 repo-local files,
and ranks each motif by a payload-retention score:

```text
score =
  sum(payload weights)
  + capped keyword-hit support
  - terminal-risk penalty
```

The hit count is only repo-local triage support.  It is capped in the score
and is not proof strength; generic numeric strings are especially treated as
weak evidence.  The actual proof-facing information is the typed payload set,
the destroyed coordinate, the repair sidecar, and the risk label.

## Result Snapshot

```text
repo_files_scanned=8370
motifs=94
families=93
risk_hist={'analogy': 19, 'direct': 33, 'raw': 3, 'sidecar': 39}
directed_3cycles=0
hamiltonian_path_count=1
```

Payload coverage in the atlas:

```text
SIDE_CARRIER        64
ANALYTIC_EQ         40
QUOTIENT_LEGALITY   30
EDGE_PACKET         19
CHEBYSHEV           18
TOEPLITZ            15
GEOMETRY            15
AP_NORMAL           14
TRAP_BOUNDARY       13
PGF_ROOT            13
COV_LAYER            9
HB_PERRON            9
SELBERG              9
CIRCUIT              9
P_ADIC               8
GREEN_LORENTZIAN     6
ORDERED_TAIL         5
```

The priority path begins:

```text
M073 comb-overlap Gram kernel
M080 shell L_y magic quartic
M067 normal-cone dual slack
M068 multi-chart proof split
M093 AP self-dual Fejer equidistribution certificate
M085 three-gap Stern-Brocot cap-kernel recursion
M001 consecutive plus doubled AP
M089 modulus-covariance apex break
M074 single-arc peeling recursion
M005 Toeplitz lambda-min margin
M066 certificate-Helly separation
M003 HYP-3204 exchange-rate ratio
M006 D1/D2/D3 layer split
M094 totally-real cap field conductor packet
M002 11 non-AP exchange traps plus AP
M029 Fejer-Riesz square
M009 Chebyshev V7 double root
M004 AP support projection
M020 Lorentzian exchange chamber
M072 conductance/Fiedler trap graph
M091 cyclotomic factor grading
M019 Green-current bottleneck
```

## Small Signals Worth Keeping

The atlas keeps these small patterns, but only with their payload labels:

```text
2-row skyline       AP and doubled AP are the all-bank equality face.
11/12 traps         11 non-AP exchange traps plus AP; Toeplitz slack discharges all 11.
3432/3431           anchored bounded k=8 rows / primitive rows.
39766/540225        AP support projection value.
12882/17161         worst central exchange-rate ratio.
6237419/8643600     total covariance Sigma_kappa2 for consecutive k=8.
6237419/25930800    ideal C6 Perron quotient lambda0.
0.042304730706      Toeplitz lambda-min margin at AP.
1/7                 apex obstruction, not an exact associator law.
3.149364            odd Worpitzky contribution dominates the even fold.
1,4,1               K3/Worpitzky Eulerian row.
1,11,11,1           next Worpitzky row.
2,6,12,20,30,42     doubled triangular numbers; only a sidecar unless bridged to a partial cube/simplex payload.
E/O legs            E=x^2+5x+4 and O=x^2+4x+1 strictly interlace.
V7 double root      Chebyshev/de Moivre cubic squared as magic-function dual.
F7 Fejer kernel     positive-definite Delsarte sidecar with weights (7-|n|)_+.
lambda2 pair        no-Toeplitz and green-only trap graphs stay connected at 2.537866286 and 1.208613477.
181 M-defect        precision M-matrix defect has many primitive beaters, so it is a guardrail, not a scalar objective.
K(1,q)=1/(7q)      speed-1 comb overlap has exact least resonance at q=13.
peeling recursion  cap(P)=cap(P\{1})-(1/7)(1-1/min(P\{1})) for speed-1 rows <=13.
-37/1092,-61/588   order-3 overlap residues are the binding j=4 and j=5 constants.
10q0+q3+10q6       finite shell L_y magic dual; quartic contact at shells 1,2,4,5.
rho>=18.019        cyclic PSD positivity overprices the central repair; guardrail only.
E7 Eisenstein      E_7=(7E_2(7tau)-E_2(tau))/6 is a level-7 coefficient sidecar.
B7 cubic           B^3-5B^2+6B-1=0 is a Beraha/Mahler height gauge, not an inequality.
7-|n|              AP autocorrelation is the Fejer/subshift rank-one transfer state.
K(a,b)=g/(7ab)     three-gap/Stern-Brocot recursion for the order-2 cap kernel.
K(a,13)            antipode column K(a,13)=(2a-1)/(91a) is verified for a=1..12.
cap_k ladder       cap_k=cap_{k-1}+k/C(2p,2) is the Faulhaber moment-order ladder.
depth=(p+1)/2      LRC(2p) apex moment depth matches the cyclotomic-degree law for p=3,5,7.
scale-normal       primitive projective shape plus first surviving coordinate is the route recursion.
K fold x1/2        modulus-covariance gives K^(2n)/K^(n)=1/2 until the apex break.
chi3 cosets        cubic mode cosets {1,6}/{2,5}/{3,4} are the de Moivre angles.
Phi_d grading      recursion mode=(x-1)^depth*Phi_d separates moment depth from character.
A..G local         signed recurrences require chart addresses before cancellation is legal.
|g(chi7)|^2=7      Gauss-sum modulus equals the Lee-Yang apex zero and Fejer reserve.
disc=49            Q(cos2pi/7) puts the binding cap rows on the 7^2 conductor.
2 heads            n=14 has a 7-cap head and a 3^3 witness head, both depth 3.
1/23 -> 1/14       Chen-Cusick supplies a floor-to-target lift; the 23/M=2/23 link is only bounded-bank coincidence.
```

The number 12 remains useful, but as chart/fiber bookkeeping: local maxima,
source witnesses, and duodecimal sectors are not interchangeable unless their
preserved predicate is named.

## Guardrail

The executed HYP-3225 trap-fingerprint scout, earlier spectral-regularization
motifs, p-adic tau valuation, Lindelof, Skewes, Collatz, Pell/cannonball,
Markov/Hurwitz, and Moser-de Bruijn/fibbinary threads are not discarded.  They
remain sidecar candidates until they attach to one of the live LRC payloads:

```text
AP normal-fan face
Toeplitz / moment-cone curvature
covariance distance layers D1,D2,D3
ordered-tail q0+q6 versus q3 pricing
HYP-3202 finite exchange-trap discharge
HYP-3214 Fejer/Delsarte F_7 positive-definite sidecar
HYP-3222 Hermite-Biehler / Perron gluing sidecar
HYP-3223 Green-current / Lorentzian trap classifiers
HYP-3225 Green-current / Lorentzian trap fingerprints
HYP-3227 conductance graph / Fiedler defect island
S75 comb-overlap Gram kernel / single-arc peeling sidecar
HYP-3230 three-gap / Stern-Brocot cap-kernel recursion
HYP-3231 scale-normal recursion ledger
HYP-3216 LRC(2p) moment-order ladder and 2-adic fold
HYP-3232 modulus-covariance apex break
HYP-3217 cyclotomic subfield / character-mode lattice
HYP-3233 cyclotomic factor grading
HYP-3234 signed address chart-change sidecar
HYP-3218 AP self-dual Fejer equidistribution certificate
HYP-3235 totally-real cap field / conductor packet
HYP-3215 induction-base and 23/27/14 modulus route
```

Incoming S283 added an additive-resonance coordination layer around Skewes,
Helfgott-Ruzsa/additive energy, Collatz 2-adic anomalies, and the PFR
landscape.  HYP-3226 absorbs that as a sidecar cluster: M023, M024, M025,
M069, and M071.  The usable part is not the famous-problem name; it is the
possibility of a Freiman-model endpoint packet, additive-energy residual, or
2-adic valuation sidecar that preserves one of the live payload coordinates.

Post-rebase HYP-3225 executes the first version of the trap-fingerprint table
that HYP-3226 asked for.  All 12 HYP-3202 arbitrary-swap local maxima select
`Toeplitz_lambda_min` as the first dictionary discharge coordinate; the 11
non-AP traps split into `rank2_pair_plucker_bottleneck` (6),
`green_low_connectivity_bottleneck` (2),
`mixed_green_lorentzian_sidecar` (1), and
`AFM_frustrated_high_rayleigh_debt` (2).  This turns M002/M019/M020 from a
next-hook into live sidecar evidence.

Post-rebase HYP-3214 verifies the cyclotomic magic-function sidecar:
`F_7=(de Moivre cubic)^2=sin^2(7t/2)/sin^2(t/2)` with Fejer weights
`F_7hat(n)=(7-|n|)_+`.  This upgrades M029 and M009 from suggestive
Chebyshev/Fejer language into a concrete positive-definite Delsarte sidecar.

Post-rebase HYP-3227 adds the full-bank conductance-graph evidence.  AP and
doubled AP tie for the main Green-current capacity coordinates; the no-Toeplitz
and green-only trap-discharge graphs remain connected with
`lambda2=2.537866286` and `lambda2=1.208613477`; and precision min-degree /
M-matrix defect appear as guardrail coordinates, not terminal scalars.  This
promotes the Green-current and conductance/Fiedler motifs to live sidecar
evidence for the finite trap discharge.

Incoming mac-mini S75 adds the comb-overlap Gram-kernel build.  The useful
payload is M073/M074/M075: `K(p,q)=meas(D_p cap D_q)` is PSD/Bochner as a
Gram matrix, `K(1,q)=1/(7q)`, speed-1 rows satisfy the peeling recursion, and
the binding rows expose order-3 residues `-37/1092` and `-61/588`.  The
hypothesis index currently also uses `HYP-3227` for this S75 magic-function
entry, while the detail file `HYP-3227-lrc14-green-current-conductance-graphs`
is the conductance scout.  Until that namespace collision is normalized,
HYP-3226 cites the Gram-kernel build by source tag S75 and motif ids M073-M075.

Incoming HYP-3215 adds the proof-route audit motifs M076-M079: verify whether
the induction base is really `n<=13` or only the published `n<=10`; start from
the Chen-Cusick `1/23` floor and lift toward the `1/14` apex target while
treating the earlier `23 = M=2/23` resonance as a bounded-bank coincidence;
recast the bounded core through polyhedron flatness or zonotope covering
radius; and import Rosenfeld-style exponential sums for Node-3.

Incoming HYP-3228/HYP-3229 add M080-M084.  M080 is the finite shell
`L_y` magic quartic with `10q0+q3+10q6`; M081 is the Gamma0(7) Eisenstein
coefficient sidecar; M082 is the Beraha/Mahler height gauge; M083 is the
subshift-transfer Perron-defect sidecar; and M084 is the Dirichlet-L/Stark
denominator guardrail.  These sharpen the atlas by separating exact shell
magic from false cyclic PSD positivity and by keeping modular arithmetic as a
coefficient engine, not a proof shortcut.

Incoming HYP-3230/HYP-3231/HYP-3216 add M085-M088.  M085 turns the
comb-overlap kernel into a three-gap/Stern-Brocot recursion
`K(a,b)=g(a,b)/(7ab)` with continued-fraction breakpoints; M086 names the
scale-normal packet tower; M087 records the verified `LRC(2p)` moment-order
depth law for `p=3,5,7`; and M088 isolates the 2-adic reflection fold that
halves the degree-4 k=8 packet.  The guardrail is sharp: the order-2 kernel
recursion is not the full cap proof until order-3 overlap and
inclusion-exclusion debt are discharged.

Incoming HYP-3232/HYP-3217 add M089-M090.  M089 names the modulus-covariance
law `K^(2n)/K^(n)=1/2` and its apex break in the antipode half; M090 names the
cyclotomic subfield/character-mode lattice, including the cubic de Moivre
`chi_3` mode.  These make the recursion layer sharper: the difficulty
concentrates at the apex fold, and the character modes are usable only after
they emit signed recursion packets or L-value sidecars.

Incoming HYP-3233/HYP-3234/HYP-3218/HYP-3235 add M091-M094.  M091 grades
recursion modes as `(x-1)^depth*Phi_d`, separating moment depth from character
factor.  M092 makes the A..G signed recurrences local chart addresses, not
global letters.  M093 identifies AP autocorrelation as the self-dual Fejer /
Vaaler equidistribution certificate with Gauss-sum reserve.  M094 records the
totally-real cap field `Q(cos2pi/7)` and the binding-row conductor debt
`7^1,7^2`.  These are direct motifs only because each names the sidecar needed
before the cyclotomic object can be used as a proof row.

The raw-famous-problem magnet is explicitly the sink motif: it destroys all
LRC payload unless it is retyped into one of those coordinates.

## Tournament Analysis

Vertices are motif families and payload atoms, not runners or arcs.  The
pairwise observable is proof-payload retention after quotienting.  The switch
orients toward the motif that preserves more of the HYP-3224/HYP-3205 normal
fan and trap-discharge coordinate set while retaining HYP-3221's analytic
equidistribution obligation.  The induced tournament is transitive in this
run (`directed_3cycles=0`, one Hamiltonian path), so the ranking is a ledger
order rather than evidence for a cyclic obstruction.

## Assumption Challenge

Alternate vertex sets considered:

```text
motif families
payload atoms
edge packets / tip-tail observers
fixed circle sections
sector boundaries
wall-crossing events
residue packets
cover arcs
Fourier modes
matroid circuits
proof obligations
```

The chosen quotient preserves the LRC predicate "this motif supplies or
repairs a named proof payload."  It destroys raw incidence multiplicity, exact
source identity, and sometimes the distinction between a scalar coincidence
and a structural certificate.  The challenged assumption is that a small
pattern is useful because it is recognizable.  The atlas instead demands that
recognizable patterns pay the no-free-slider tax.

## Next Proof Target

The next high-value step is now the symbolic version of the HYP-3225 table.
For each trap class, prove or derive:

```text
first failed HYP-3205 dictionary coordinate
comb-overlap Gram-kernel PSD coordinate
speed-1 peeling-recursion status
order-3 triple-overlap residue
Toeplitz lambda-min slack
Green-current bottleneck type
conductance-graph lambda2 / Fiedler-cut subcase
Lorentzian / valuated-exchange defect
precision M-matrix defect / Schur-complement edit
Worpitzky or Hermite-Biehler sidecar debt
Fejer/Delsarte dual slack against the F_7 kernel
shell L_y magic quartic slack against 10q0+q3+10q6
Gamma0(7) coefficient-row compatibility
Beraha/Mahler/subshift/Dirichlet-L sidecar status
three-gap Stern-Brocot kernel-recursion status
scale-normal packet / omega_Q exactness status
LRC(2p) moment-order ladder and 2-adic fold status
modulus-covariance apex-break status
cyclotomic subfield / chi_3 mode status
cyclotomic factor grading and Phi_d character status
signed address chart-change legality status
AP self-dual Fejer/Vaaler tail status
totally-real cap field conductor / trace status
induction-base dependency and Chen-Cusick floor-to-1/14 lift status
```

If those columns collapse to exact identities or finite inequalities, HYP-3226
turns the small-pattern atlas into a finite discharge recipe rather than a
loose list of analogies.
