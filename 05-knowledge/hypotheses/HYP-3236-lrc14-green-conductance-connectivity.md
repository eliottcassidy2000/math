---
id: HYP-3236
title: LRC14 Green conductance and algebraic connectivity certificate
status: EVIDENCE / exact bounded-bank scout; not an LRC14 proof
source: codex-2026-06-28
tangent: T1336
technique: LTI-336
tournament_technique: LTT-236
script: 04-computation/lrc14_green_conductance_connectivity_codex_20260628.py
result: 05-knowledge/results/lrc14_green_conductance_connectivity_codex_20260628.out
reflection: 07-reflections/lrc14-green-conductance-connectivity-codex-20260628.md
related:
  - HYP-3235
  - HYP-3234
  - HYP-3233
  - HYP-3232
  - HYP-3231
  - HYP-3230
  - HYP-3229
  - HYP-3228
  - HYP-3227
  - HYP-3226
  - HYP-3225
  - HYP-3224
  - HYP-3223
  - HYP-3222
  - HYP-3221
  - HYP-3218
  - HYP-3217
  - HYP-3214
  - HYP-3216
  - HYP-3213
  - HYP-3212
  - HYP-3211
  - HYP-3210
  - HYP-3205
  - HYP-3204
  - HYP-3203
  - HYP-3202
  - HYP-3201
  - HYP-3200
  - HYP-3163
  - HYP-3162
  - HYP-3161
  - HYP-3160
  - HYP-3154
  - HYP-3153
  - OPEN-Q-108
---

# HYP-3236: LRC14 Green Conductance And Algebraic Connectivity Certificate

## Claim

The electrical half of HYP-3223 is not only a metaphor.  On the exact bounded
k=8 bank, the positive covariance conductance graph gives a new AP-tight
coordinate:

```text
C_E(i,j) = Cov(1_{sector i empty}, 1_{sector j empty})
G_+(E)   = weighted graph with conductance max(C_E(i,j),0)
L_E      = graph Laplacian of G_+(E)
L_E^+    = Green kernel / Moore-Penrose inverse
```

Consecutive speeds and the nonprimitive doubled AP dilation are the only rows
that maximize algebraic connectivity `lambda2(L_E)` and total positive
conductance, and the only rows that minimize Kirchhoff index, mean effective
resistance, maximum effective resistance, and each cyclic-distance effective
resistance channel.

This is evidence for a new proof coordinate:

```text
AP/consec = maximum-conductance / maximum-spectral-gap network
HYP-3202 traps = finite bottleneck networks
negative covariances = leakage sidecar, not disposable noise
```

## Exact Scout

The script recomputes the anchored bounded bank

```text
E = {0} union A, A subset {1,...,14}, |A|=7
rows_all = 3432
rows_primitive = 3431
```

using the six inner sectors as vertices.  Edge conductance is the positive
part of pair covariance.  The positive-part operation is deliberately lossy:
negative covariance edges are retained as `negative_edges` and
`negative_mass`.

The AP profile is:

```text
total_positive_conductance = 0.721622819196
negative_edges = 0
negative_mass = 0
lambda2_algebraic_connectivity = 0.192033074001
laplacian_eigvals =
  0, 0.192033074001, 0.278704158990,
  0.291269950895, 0.327480748451, 0.353757706055
kirchhoff_index = 108.654718079151
mean_effective_resistance = 7.243647871943
max_effective_resistance = 9.713313375596
distance_R1 = 6.710085829826
distance_R2 = 7.345147613359
distance_R3 = 7.675710172645
bottleneck_pair = (1,6)
bottleneck_current_entropy_bits = 3.460234326390
```

Across the all-bank ranks:

```text
lambda2_algebraic_connectivity_MAX:
  all_beaters=0 primitive_beaters=0 all_ties=2 primitive_ties=1
total_positive_conductance_MAX:
  all_beaters=0 primitive_beaters=0 all_ties=2 primitive_ties=1
kirchhoff_index_MIN:
  all_beaters=0 primitive_beaters=0 all_ties=2 primitive_ties=1
mean_effective_resistance_MIN:
  all_beaters=0 primitive_beaters=0 all_ties=2 primitive_ties=1
max_effective_resistance_MIN:
  all_beaters=0 primitive_beaters=0 all_ties=2 primitive_ties=1
distance_R1_MIN, distance_R2_MIN, distance_R3_MIN:
  all_beaters=0 primitive_beaters=0 all_ties=2 primitive_ties=1
```

The only all-bank equality rows are AP and the doubled AP dilation; primitive
normal form leaves AP unique.  This matches HYP-3205 and HYP-3224: the new
Green coordinate exposes the same two-row face.

## Decoys And Traps

The closest non-AP Green decoy by algebraic connectivity is
`(0,2,3,4,5,6,7,8)`:

```text
lambda2 = 0.144321509290
kirchhoff = 136.424340938360
maxR = 12.847664881121
negative_edges = 1
negative_mass = 0.004453237077
```

Other high-connectivity decoys already show the split that a proof must
respect.  Some have no negative covariance edges but still lose the Green
certificate through lower `lambda2` and larger resistance.  Others carry
visible negative leakage.  Therefore "all positive covariances" remains too
weak, just as HYP-3202 warned.

For the HYP-3202 arbitrary one-point exchange traps, the scout finds `12`
primitive local maxima: AP plus `11` non-AP traps.  Every non-AP trap has a
strict effective-resistance defect relative to AP.  The primary discharge
histogram is:

```text
non_AP_trap_primary_green_discharge =
  {'kirchhoff_excess': 3, 'maxR_excess': 8}
```

Examples:

```text
(0,2,4,6,7,8,10,12):
  lambda2=0.113820536045
  kirchhoff=173.325358101925
  maxR=15.730609041900
  negative_mass=0.010485792957
  primary=maxR_excess

(0,1,3,5,7,9,11,13):
  lambda2=0.110831526196
  kirchhoff=195.383171490455
  maxR=16.111973010442
  negative_mass=0
  primary=kirchhoff_excess

(0,2,5,7,9,10,12,14):
  lambda2=0.016336734005
  kirchhoff=606.795921265124
  maxR=87.637153930345
  negative_mass=0.080333862948
  primary=maxR_excess
```

This makes the HYP-3202 finite trap manifold more geometric: the traps are
not only local exchange exceptions; they are networks with named bottlenecks.

## Rebase Integration: HYP-3231 Scale Ledger, HYP-3232 Interlocking Recursions, And HYP-3216 Family Law

Incoming mainline renumbered the universal scale-invariance recursion ledger
to HYP-3231, kept HYP-3230 as the three-gap cap-kernel recursion thread, added
HYP-3232 for the Mobius/Eisenstein/Legendre interlocking-recursion audit, and
sharpened HYP-3216 into the verified two-recursion family law for the LRC(2p)
route:

```text
HYP-3230: three-gap / Stern-Brocot recursion in the cap kernel K(a,b)
HYP-3231: universal scale-normal quotient ledger
HYP-3232: Mobius/Eisenstein/Legendre recursion interlock at apex 7
HYP-3216: moment-order ladder x 2-adic reflection fold across n=2p
HYP-3236: finite Green conductance / algebraic-connectivity face
```

This is a real constraint on the Green packet.  A raw spectral scalar such as
`lambda2` cannot remember the Farey or three-gap address of the cap kernel.
Nor can it see by itself where HYP-3232's modulus-covariant regime breaks at
the apex half `n/2=7`, precisely where the antipode speeds become the hard
LRC14 half.
The lawful scale-compatible object has to be a function packet:

```text
row -> sector gaps / cap-kernel address
    -> covariance conductance edges
    -> cut/Rayleigh/Thomson tests
    -> Green resistance sidecars
```

The conjectural bridge is that AP is not merely the best-connected finite
network; it is the fixed point where the pair-Pascal cap mass, the Fejer
autocorrelation, and the Green Dirichlet form all become scale-normal.  A proof
therefore needs either a renormalization law for conductance cuts or a
certificate that the Green slack is annihilated by the HYP-3231/HYP-3232/
HYP-3216 recursion coordinates.  Otherwise algebraic connectivity is only a
successful finite diagnostic, not yet a route proof.

Incoming HYP-3217 then identifies the Mobius/Eisenstein/Legendre/CUBIC/sextic
stack with the subfield lattice of `Q(zeta_7)` times the 2-adic fold, with the
cubic Gaussian-period mode giving the de Moivre angles.  For this Green
packet, that says the Fiedler vectors, bottleneck currents, and
distance-resistance channels should be projected onto the cubic-period cosets
before they are allowed to count as proof data.  If the Green slack is
scale-normal but not cubic-mode-aware, it has still compressed away a live
cyclotomic coordinate.

Incoming HYP-3233 then recasts the recursion modes themselves as cyclotomic
factors `(x-1)^depth * Phi_d`, with the moment-order depth carried by the
unipotent `(x-1)` multiplicity and the hard apex mode carried by the real
`Phi_7`/de Moivre cubic factor.  For HYP-3236 this adds a spectral algebra
test: the Laplacian spectrum, Green resolvent, Fiedler mode, and effective
resistance channels should be tagged by which cyclotomic factor they preserve
or forget.  A conductance graph that sees a large gap but erases the `Phi_7`
mode is still an illegal compression for the proof route.

The concurrent HYP-3234 signed-address chart-change sheaf gives the same
warning in coordinate language.  Full, even, and odd recurrences are local
signed charts with slot bases and cancellation debts, not global letter
formulas.  Green reductions should therefore retain
`signed_address_chart`, `local_slot_basis`, `chart_change_map`, and
`cancelled_same_size_slots` when Fiedler/current data is transported across
recursion charts.  Otherwise a clean electrical scalar may be hiding an
illegal chart-change compression.

Incoming HYP-3235 and HYP-3218 then push the harmonic side closer to proof:
the cap lives in the totally real cyclotomic field `Q(cos(2*pi/7))`, the
Fejer magic function is a totally positive cyclotomic square, and the
equidistribution finish is recast as an explicit cyclotomic Fejer certificate.
For HYP-3236 this suggests a direct comparison: the Green kernel and effective
resistance slack should be tested against the Fejer square / Gauss-sum margin
instead of only against the finite bounded-bank rankings.  If the Green face is
real proof data, it should either factor through that totally-real positive
certificate or name the residual electrical sidecar.

## Rebase Integration: HYP-3227 Conductance Graphs

Incoming mainline claimed HYP-3227/T1325/LTI-325/LTT-225 for the broader
Green-current conductance-graph scout while this packet was local.  The two
packets should not collide:

```text
HYP-3227: sector graph + precision graph + trap-discharge conductance graph
HYP-3236: positive-covariance Green resistance extremality coordinate
```

HYP-3227 tests all-ones Green energy, positive/precision graph capacities,
M-matrix defects, and trap-certificate algebraic connectivity.  HYP-3236
keeps the positive covariance graph as an AP-tight Dirichlet object and pushes
harder on `lambda2`, `L^+`, Kirchhoff, distance effective resistances,
bottleneck unit currents, current entropy, and the negative-leakage sidecar.
Thus HYP-3227 supplies the richer graph family; HYP-3236 supplies a sharper
Green-resistance extremality readout for one legal face of that family.

## Rebase Integration: HYP-3228 Shell Magic

Incoming mainline also claimed HYP-3228/T1326/LTI-326/LTT-226 for the
cyclotomic Delsarte shell-magic packet.  That packet identifies the exact k=8
`L_y` dual as the finite shell polynomial

```text
f(n) = ((n-1)(n-2)(n-4)(n-5))/4
E[f(N)] = 10q0 + q3 + 10q6 = 10L_y
```

So the current neighborhood now has several distinct faces:

```text
HYP-3232: Mobius/Eisenstein/Legendre interlocking recursion audit
HYP-3231: scale-normal quotient ledger
HYP-3230: three-gap/Farey cap-kernel recursion thread
HYP-3227: conductance/precision/trap-discharge graph family
HYP-3228: Delsarte shell-magic / L_y coefficient face
HYP-3236: positive-covariance Green-resistance extremality face
```

The proof task is to compare their slacks without collapsing them.  HYP-3228
prices the `q3` center through the Delsarte shell; HYP-3236 prices bottleneck
resistance through the Green kernel.  If they are the same normal-fan face,
that has to be proved by a packet map, not asserted from shared AP equality.

## Rebase Integration: HYP-3229 Modular Magic Sidecar

Incoming mainline then claimed HYP-3229 for a modular/arithmetic sidecar audit
around the Fejer route: Gamma0(7), Dirichlet-L, Beraha/Mahler, subshift, and
the comb-overlap Gram kernel.  That packet is an arithmetic coefficient-engine
audit, not the Green-resistance certificate itself.  HYP-3236 should use it
only after naming which coefficient sidecar is being imported into the
conductance/Dirichlet-form packet.

## Rebase Integration: HYP-3225 Trap Fingerprints

Concurrent mainline work claimed HYP-3225 for the trap-local Green/Lorentzian
fingerprint scout.  That packet evaluates the `12` HYP-3202 arbitrary-swap
local maxima plus one-swap neighbors and finds that Toeplitz `lambda_min` is
still the universal first dictionary discharge.  Its residual sidecar classes
are:

```text
rank2_pair_plucker_bottleneck: 6
green_low_connectivity_bottleneck: 2
AFM_frustrated_high_rayleigh_debt: 2
mixed_green_lorentzian_sidecar: 1
```

HYP-3236 is the complementary all-bank Green coordinate.  It proves, at the
bounded-bank evidence level, that AP/doubled AP are the only rows tight for
algebraic connectivity and Green resistance profiles across all `3432` rows.
HYP-3225 then refines the finite trap boundary by explaining which local
sidecar mechanism lives under the universal Toeplitz discharge.  The packets
should be read as:

```text
HYP-3236: global conductance/connectivity extremal face
HYP-3235: totally-real cap / Fejer-square cyclotomic certificate
HYP-3234: signed-address chart-change sheaf
HYP-3233: cyclotomic-factor grading
HYP-3232: interlocking recursion / apex-break audit
HYP-3231: scale-normal quotient ledger
HYP-3230: three-gap/Farey cap-kernel recursion thread
HYP-3217: subfield-lattice / cubic de Moivre mode
HYP-3229: modular/arithmetic magic sidecar audit
HYP-3227: broader conductance/precision/trap-discharge graph family
HYP-3228: Delsarte shell-magic / L_y coefficient face
HYP-3225: trap-local typed sidecar taxonomy
HYP-3226: reserved motif atlas that should ingest both
```

## Rebase Integration: HYP-3214 Fejer Magic Function

Incoming HYP-3214 identifies the cyclotomic Delsarte/Beurling-Selberg magic
function with the order-7 Fejer kernel:

```text
F_7(theta) = sin^2(7 theta / 2) / sin^2(theta / 2)
           = (de Moivre cubic)^2
hat F_7(n) = (7 - |n|)_+
```

This is the harmonic sibling of the Green packet.  HYP-3214 says AP is sharp
because its 7-sector autocorrelation is the positive-definite Fejer kernel
with double zeros at the de Moivre angles.  HYP-3236 says the same AP row is
sharp after the finite covariance packet is clipped to a conductance graph and
read as a Dirichlet form via `L`, `L^+`, `lambda2`, and effective resistance.

The shared proof target is therefore not a scalar invariant but a commuting
diagram of functions:

```text
AP orbit / autocorrelation -> Fejer positive-definite kernel -> harmonic dual
sector covariance packet   -> conductance Dirichlet form      -> Green dual
cap kernel address         -> scale-normal recursion          -> family law
```

The warning from HYP-3214 is essential: `F_7` is not the 14-clock cap
functional.  The cap lives in the Johnson/pair-Pascal magic function, while
`F_7` lives on the 7-sector coverage side.  HYP-3236 must therefore keep the
Green conductance graph, negative leakage, odd sidecars, and pair-Pascal cap
as distinct coordinates until a dual certificate explicitly glues them.

## Compression And Information View

The construction is a chain of functions:

```text
row E
  -> covariance matrix C_E
  -> positive conductance graph G_+(E)
  -> Laplacian L_E
  -> Green kernel L_E^+
  -> lambda2, resistance, bottleneck-current entropy
```

Each arrow is a compression.  The legal quotient test from HYP-3201 applies:
the proof predicate may pass through the compression only if destroyed
information is zero-defect, reconstructible, dual-annihilated, or named as
residual debt.

The most important destroyed coordinate is the negative part of `C_E`.
Discarding it turns a signed covariance kernel into a conductance graph.  That
is useful because Laplacian tools become available, but it is not lawful by
itself.  HYP-3236 therefore treats

```text
(G_+(E), negative_edges(E), negative_mass(E), odd/Worpitzky sidecars)
```

as the actual packet.  Plain `G_+(E)` is a shadow.

There is also a noncommutative/nonassociative warning.  Network reduction by
Schur complement is associative only when the boundary terminals and eliminated
coordinates are retained.  If one compresses first to a scalar such as
`lambda2`, then star-mesh edits, exchange moves, and trap discharge no longer
compose.  This is the same failure-of-compression pattern as HYP-3201's
commutativity/associativity examples: a law that is true in a lifted packet can
fail after quotienting away the bracket, order, or boundary sidecar.

## Relation To HYP-3224 And HYP-3222

HYP-3224 says the live proof object is a normal-fan payload cube.  HYP-3236
adds a new visible face:

```text
root-locus / AP support face
moment-cone / Toeplitz lambda_min face
coefficient / ordered-tail exchange face
ferromagnetic / covariance-layer face
spectral / Perron-HB-Joukowski face
scale-recursion / three-gap cap-kernel face
Green-conductance / algebraic-connectivity face
```

The Green face is compatible with the Perron side of HYP-3222.  A conductance
Laplacian has a constant ground state and a spectral gap `lambda2`; the HYP-3222
ideal C6 Perron quotient has the all-ones mode as top covariance mode.  These
are dual readings of the same coherent phase: covariance wants the all-ones
current, while the conductance Laplacian wants no bottleneck in the orthogonal
directions.

The Hermite-Biehler/Joukowski side still matters because conductance is an
even-degree shadow.  It does not see the odd Worpitzky associator by itself.
A proof must glue the Green face to the HYP-3222 even/odd interlacing sidecar
and the HYP-3204 ordered-tail pricing of central `q3`.

## Proof Frontier

The theorem-facing target is a Rayleigh/Thomson extremality statement:

```text
For every primitive legal k=8 row E,
  either G_+(E) has lower algebraic connectivity than AP,
  or its Green resistance profile has a named bottleneck excess,
  or the missing negative/odd sidecar discharges the row through another chart.
Equality forces AP; all-bank equality also allows the doubled AP dilation.
```

Possible proof routes:

1. Prove a conductance rearrangement lemma: AP maximizes every positive
   covariance conductance cut or a sufficient family of cuts, then use
   Rayleigh monotonicity / Cheeger / Poincare to force `lambda2`.
2. Prove a Thomson dual: every non-AP row has a unit-current demand whose
   energy is larger than AP, with HYP-3202 traps handled by the finite
   bottleneck list.
3. Lift the Green packet through the HYP-3231/HYP-3232 scale-recursion stack:
   express the relevant conductance cuts as functions of the three-gap/Farey
   cap-kernel address, then test whether the HYP-3216 moment-order ladder and
   the Eisenstein apex-break preserve the same Rayleigh slack across LRC(2p).
3. Glue the Green coordinate to the HYP-3224 normal fan: show Green slack is
   controlled by AP support, Toeplitz slack, or distance-layer slack.
4. Lift Schur-complement network edits to the HYP-3223 Lorentzian/valuated
   exchange language; a legal exchange move should reduce a Green bottleneck
   or emit a named circuit sidecar.

## Tournament Analysis

Vertices are conductance certificates and proof obligations, not runners,
gaps, or raw sectors:

```text
algebraic_connectivity_lambda2
kirchhoff_green_resistance_profile
distance_effective_resistance_channels
schur_complement_trap_discharge
positive_covariance_conductance_graph
negative_covariance_leakage_sidecar
spectral_dictionary_compatibility
raw_total_positive_conductance
raw_covariance_scalar
plain_positive_association
runner_gap_graph
```

Pairwise observable: which carrier preserves Green current plus leakage
sidecar payload.

Switch/gauge: orient from lower payload retention toward higher payload
retention; ties follow the displayed priority path.

Fingerprint:

```text
score_hist = {0:1,1:1,2:1,3:1,4:1,5:1,6:1,7:1,8:1,9:1,10:1}
directed_3cycles = 0
sccs = singletons
hamiltonian_path_count = 1
priority_path =
algebraic_connectivity_lambda2
-> kirchhoff_green_resistance_profile
-> distance_effective_resistance_channels
-> schur_complement_trap_discharge
-> positive_covariance_conductance_graph
-> negative_covariance_leakage_sidecar
-> spectral_dictionary_compatibility
-> raw_total_positive_conductance
-> raw_covariance_scalar
-> plain_positive_association
-> runner_gap_graph
```

## Assumption Challenge

Alternate vertices considered: runners, gaps, sector boundaries, wall
crossings, residues, covariance edges, conductance graphs, current paths,
Fiedler modes, Schur-complement edits, trap rows, and proof obligations.

Chosen vertices: conductance graphs and Green-current certificate carriers.

Preserved predicate: AP covariance extremality, distance-layer payload, trap
identity, effective-resistance bottlenecks, algebraic-connectivity margins,
and the AP/doubled-AP equality face.

Destroyed information: negative covariance signs unless leakage sidecar is
retained; raw runner/gap identities; odd Worpitzky debt; exchange chamber
sidecars.

Challenged assumption: positive association or a raw conductance scalar is
enough.  The graph must retain where current flows, where the Fiedler
bottleneck sits, and where negative covariance was clipped.
