# Session139 verification and checkpoint lineage

**Independently accepted scope through this checkpoint:** exact weighted stability at n=3..9 and n>=48;
all-modulus contexts/clocks/content; finite complete prime states and weighted
simulation at3,5,7,11; integral/weighted cocycle incidence; the exact stated
regular isospectral-circle range. The first weighted spectral boundary is n5.
**10..47 is unresolved by this reviewed work.** No RH/GRH or integer-spectrum closure is claimed.

| Checkpoint | Scope | Evidence |
|---|---|---|
| 1, published `400dcd629598c30762eb9adb8bb85bea581698a7` | n>=748, weighted n5 spectral loss, corrected historical recovery | [Frozen scope](CHECKPOINT1.md), [original replay](REPLAY_CHECKPOINT1.json), [original source integrity](SOURCE_INTEGRITY_CHECKPOINT1.json); 42 published owned blobs verified |
| 2, published `9c7eb0b13e48e68b8fd076feccbf56d04cdd54aa` | n>=526; all-modulus functions/clocks/padding; exact finite prime quotients and actual5/11 projectors | [Frozen scope](CHECKPOINT2.md), [replay](REPLAY_CHECKPOINT2.json);86 cumulative owned blobs verified |
| 3, published `b68db2afda029c83ad7087f09fdb9297037f8414` | Integral/weighted cocycle, regular critical-circle obstruction, context content/simulation and narrow historical repair | [Frozen scope](CHECKPOINT3.md), [fresh replay](REPLAY_CHECKPOINT3.json), separately pinned independent audits;153 cumulative owned blobs verified |
| 4, published `8da4684eef0e8f35041889ae4f0b5dbea37e1a40` | Tensor product margin and independent arithmetic proof; n>=140; prime padding orbits; separate all-modulus tensor-coset appendix; narrow order-eight naming repair | [Frozen scope](CHECKPOINT4.md), [fresh replay](REPLAY_CHECKPOINT4.json), [final handoff](FINAL_HANDOFF.md);246 cumulative owned blobs verified |
| 5 | Independently audited n68 alternate; independent acceptance of concurrent n56 and slot-aware n48 frontiers | [Scope](CHECKPOINT5.md), [proof](concurrent68/PROOF.md), [fresh replay](REPLAY_CHECKPOINT5.json) |

All earlier scientific lane files, audit reports and replay receipts remain
byte-unchanged. Root README/SYNTHESIS/BOARD/VERIFICATION/FINAL_HANDOFF and SHA256SUMS are
maintained summaries/manifest. `CHECKPOINT1_SHA256SUMS` preserves the original
whole-packet manifest **at checkpoint-one's commit**; `CHECKPOINT2_SHA256SUMS`
does the same at checkpoint two, and `CHECKPOINT3_SHA256SUMS` at checkpoint
three; `CHECKPOINT4_SHA256SUMS` records the published fourth checkpoint.
These are not commands to
validate later rewritten summaries against old hashes. The original
`PUBLIC_FILES.json` is likewise the historical checkpoint-one selector.
The root `SHA256SUMS` validates the cumulative current packet.

## Reviewed snapshot reaches n>=48 without rewriting upstream sources

[Concurrent source](../stability_748_guarded_encodings_20261102/STABILITY_48.md) · [Independent audit](independent_concurrent48/AUDIT.md) · [Immutable pins](CONCURRENT48_SOURCE_PINS.json)

**ACCEPTED after independent argument and equality review.** This packet's
reviewed range is3..9 and n>=48; later global routers may point to newer
external claims without those claims being automatically endorsed here.
The independent68 alternate and reviewed56 predecessor remain unchanged.

The48 proof credits the accepted tensor reserve and stronger degree bound,
then adds exact matching slots and the maximum of two valid group invoices,
`max(w*S,w*S/2+S^2/4)`. It does not add those overlapping invoice bounds.
The two-regime scalar envelope and polynomial identities are checked exactly;
the slot constraint forces12 heavy edges with strict endpoint slack1/100.
Uniform tail monotonicity and the full equality transport are reviewed.
Any optional general scalar iff statement is interpreted for d>0; the actual
theorem fixes d=5/6. No source correction is needed.

The independent implementation has **47 active checks and760 exact graph
controls**, with continuum-polynomial reasoning rather than extrapolation
from a grid. Its47 endpoint fails the first bootstrap advance. A false24/5
scalar ratio gives the exact negative minimum -1/288; this is a local
counterexample to that stronger invoice, not to the stability theorem at47.

```bash
# From this packet in a checkout containing the pinned upstream source:
python3 replay_concurrent48.py --output /tmp/session139-concurrent48-replay.json
```

The fresh wrong-working-directory replay runs **16 executions** normally and
with optimization: four scientific positives, four source-integrity positives
and eight expected rejections. Independent outputs and both source verifiers
match frozen receipts; same-size audit corruption is rejected. The unmodified
pinned producer runs only in a temporary directory, reproducing25,351 local
checks,729 whole-graph controls and277 finite-R controls. These are checking
evidence, not the all-order proof. [Replay](REPLAY_CONCURRENT48.json).

Across the credited combined work,1452 of the original1490 unresolved orders
are now certified in this snapshot; **10..47 (38 orders) remains unresolved**.
No threshold-optimality or phase-transition claim is made. The global router
bytes from the newer parent are preserved; all29 non-owned incoming changed
entries and all other parent entries are inherited exactly, without ledger
queries. [Latest reconciliation](CHECKPOINT5_RECONCILIATION_48.md).

## Reviewed n56 predecessor: independently audited n>=56

[Concurrent source](../stability_priority_spectrum_20261102/STABILITY_56.md) · [Independent audit](independent_concurrent56/AUDIT.md) · [Immutable pins](CONCURRENT56_SOURCE_PINS.json)

**ACCEPTED after independent review.** The original source's same-agent audit
wording is preserved as history, not silently rewritten. This release adds
an independently written proof review and40 active endpoint checks. The
source combines Session139's accepted tensor reserve with the credited
stronger discarded-degree inequality and flexible retained-baseline dual.
All contributions, earlier proofs and their independent ranges remain clear;
the n68 proof is a valid alternate rather than the current best cutoff.

The independent bootstrap uses different Taylor degrees4,5,5 and the weaker
penalty1/10000 instead of copying the producer's degree12 certificate. The
local dual is checked coefficient-wise, its remainder is nonnegative over
the stated domain, uniform tail monotonicity is explicit, and the complete
boundary/equality chain is reviewed. Threshold55 fails the first reserve
advance. The bad-dual control keeps the square identity but violates the
remainder sign, producing a genuine local counterexample. Neither control
is a counterexample to the stability theorem at55.

```bash
# From the packet in a checkout containing the pinned upstream source:
python3 replay_concurrent56.py --output /tmp/session139-concurrent56-replay.json
```

The fresh wrong-working-directory replay has **16 executions** in normal
and optimized Python: four scientific positives, four source-integrity
positives and eight expected rejections. The independent outputs match the
frozen receipts. The unmodified pinned producer runs only in a temporary
directory and its4351-gate output matches both frozen producer receipts;
no upstream file is written. Both source verifiers reject same-size audit
corruption. [Replay](REPLAY_CONCURRENT56.json). The finite producer scan is
verification evidence, not a substitute for the uniform proof.

This predecessor certifies **3..9 and n>=56**, covering1444 of the initially1490
unresolved orders. The independently reviewed48 source above extends that
range; the56 source is not rewritten to claim48.
The proof does not claim threshold optimality. All32 newly concurrent paths
are inherited; the29 non-owned changed entries remain exact parent entries.
The n56 source, THM-1780 correction and opaque ledger are not rewritten.
[Final reconciliation](CHECKPOINT5_RECONCILIATION.md).

## Fifth checkpoint alternate: a credited proof combination reaches n>=68

[Combined proof](concurrent68/PROOF.md) · [Independent audit](independent_concurrent68/AUDIT.md) · [Final source review](independent_concurrent68/FINAL_SOURCE_REVIEW.md)

**ACCEPTED, including the stated equality.** The factor-four stronger
`W>=2D^2/n` degree cancellation and retained `-RB` light budget are explicitly
credited to the concurrent stability-priority-moduli packet. The combination
reproves those used lemmas and imports Session139's independently accepted
tensor product reserve. This is not a new claim of independent verification
of the whole concurrent n276 theorem.

Each discarded edge is paid by its assigned greedy selected/discarded wedge,
canceling the diagonal subtraction in the discarded-degree square identity.
The assigned and discarded/discarded wedge sets are disjoint. Stronger square
completion then bootstraps the reserve `6/25 -> 63/100 -> 49/50 -> 99/100`
for R>=34. After retaining `-RB`, the exact light-budget margins are both
`43/71400`; the geometric-tail endpoint has slack
`1854611611/766651662336`. Explicit monotonicity supplies all larger orders.
The proof preserves the perfect-matching-only correction and binary baseline
once, and handles zero variance, zero discarded mass and the equality locus.

The independent exact verifier has **38 active conditions** and exhausts760
weighted simple graphs through order4 with edge weights in `{0,1/8,1/4}`.
These test the cancellation identity and joint budget separately; finite
controls are not substituted for the all-real-weight proof. Normal and
optimized outputs match. The false67 certificate fails its first bootstrap
mass check; omitting selected-edge payment fails a three-vertex path. These
are rejected proof shortcuts, not counterexamples to stability at67.

Final source pins and the reviewed definition/status changes are recorded
without modifying the original audit. Earlier748/526/140 proof bytes, all
concurrent proof packets and guarded-FRACTRAN precision work remain intact.
No canonical file or central ledger changes in this checkpoint.

```bash
# From this packet directory:
sha256sum -c SHA256SUMS
python3 replay_checkpoint5.py --output /tmp/session139-checkpoint5-replay.json
```

The fresh wrong-working-directory replay performs **14 executions** under
normal and optimized Python: two scientific positives, four source-integrity
positives and eight expected rejections (67, missing selected payment and
same-size proof corruption rejected by both source verifiers in both modes). Output bytes match the frozen
receipts; source-integrity comparisons are active with optimization.
[Replay receipt](REPLAY_CHECKPOINT5.json) · [Source pins](SOURCE_PINS_CHECKPOINT5.json).

This independently audited alternate certifies **3..9 and n>=68**, covering
1432 of the1490 initially unresolved orders. The later concurrent n56/n48 audits
above supply stronger reviewed ranges; this proof is not rewritten to
claim56. Neither cutoff is asserted optimal.

## Fourth checkpoint: the genuine dimension-improving dependency

[Primary tensor proof](tensor_matching/PROOF.md) · [Independent review](stability/INDEPENDENT_TENSOR_MATCHING_REVIEW.md) · [Independent arithmetic proof](stability/EISENSTEIN_BERNSTEIN_TENSOR_PROOF.md) · [Exact audit](tensor_matching_audit/AUDIT.md)

**ACCEPTED.** Odd constant/even other coefficients imply the sharp tensor
norm floor on the real unit cube. Grouping all even principal Pfaffians
and observer counterparts compensates cancellations; it does not assert
that each individual Pfaffian is bounded below. Rounding the complement
of a chosen matching preserves its variances and yields the weight-sensitive
product margin. The matching restriction, perfect-matching penalty and
binary floor are explicit; a product over an arbitrary uncertainty graph
is not asserted. Ordinal weighted pairs attain the bound, without claiming
a complete equality classification for this standalone floor.

The independent exact audit covers12,604 exhaustive low-dimensional tensor
controls,300 higher-dimensional controls,876 rational tournaments through
order9,7,394 complete grouped transforms,20 ordinal equalities and eight
false-strengthening hostiles. The separate Eisenstein--Bernstein route
retains all coefficient slices and uses nonvanishing in F4; it is not a
fixed-polynomial evenness shortcut. Producer controls comprise762 checks,
including18 phase and484 coefficient cases; independent arithmetic-slice
checks reconstruct the transforms rather than importing producer code.

[N140 proof](stability/PROOF_N140.md) · [Independent audit](independent_threshold140/AUDIT.md) · [Final status bridge](independent_threshold140/FINAL_STATUS_REVIEW.md)

**ACCEPTED, including the stated equality.** The new product reserve gives
`B+1-exp(-M)-4^-R`, with the product penalty used only for a perfect matching.
For R>=70, exact Taylor bounds bootstrap the retained reserve from
`249999/1000000` through `7/25` and `57/100` to `99/100`. The light-deficit
dual then has exact final slack `1/448`. The46 named independent conditions
remain active under optimization. The false139 certificate and altered
perfect-matching indicator are both rejected: neither is a counterexample
to the unresolved theorem at139. The original audit pin and final proof pin
differ only in the status line; reconstruction and the final review bridge
are recorded in [source pins](SOURCE_PINS_CHECKPOINT4.json).

Thus the causal chain is **binary odd-Pfaffian parity -> tensor compensation
(with an independent arithmetic proof) -> matching product reserve ->
bootstrap -> n>=140**. The former one-quarter dual barrier limited the older
input, not all proof strategies. The RH-adjacent spectral work is a rigorous
information-loss boundary, not this dimension-improving mechanism.

[Local-factor hostiles](stability/LOCAL_FACTOR_HOSTILES.md) reject the tempting
conditional factors: the five-vertex fractional-cross residual is `-7/8192`,
and the eight-vertex uncertain-pair/two-triangle residual is `-q`. Independent
literal checks agree. The old hostile note's pending absolute-floor question
is now answered by the new tensor proof; [scope reconciliation](stability/CHECKPOINT4_PRODUCER_SCOPE.md)
preserves that chronology instead of rewriting the old note.

## Fourth checkpoint: padding orbits and the separate coset extension

[Prime padding proof](all_moduli/PRIME_PADDING_ORBITS.md) · [Independent review](independent/PRIME_PADDING_ORBITS_AUDIT.md)

**ACCEPTED.** Formal odd-prime parity-function fixed counts are1 for odd
iterates and `2*p^(gcd(d,p-1)/2)-1` for even iterates. Mobius inversion yields
periodic-state and orbit counts. Complete actual realization remains restricted
to3,5,7,11; formal higher-prime spaces are not relabeled as actual tournaments.
The all-modulus reduced-coordinate operator and weighted operator are distinct.
The historical ten-node metagraph/eight-profile cardinality firewall remains
explicit, with all three historical blob pins rechecked.

The verifier covers24 formal primes through97,2,068 root counts,131 matrix
moduli,3,009 basis instances and202,222 operator updates. The independent
work is a proof review plus replay of the producer's two-formula verifier,
not an independent new census or second implementation of the full audit.

[All-modulus tensor appendix](tensor_matching_moduli/APPENDIX.md) · [Independent acceptance](tensor_matching_moduli/REVIEW.json)

**ACCEPTED separately.** If the constant coefficient lies in `rho+MZ` and
all others in `MZ`, the same induction proves the real norm floor
`dist(rho,MZ)^2 * product(1+x_i^2)`, for every integer M>=2 and unit-cube point.
The zero-coset case, ties and composite moduli need no division. Exact controls
include91,212 low-dimensional examples,720 mixed-point/transform comparisons,
9,486 equality cases and four invalid-assumption hostiles. Ordinary binary
Pfaffian grouping supplies only the M=2 coset; no stronger tournament bound,
further cutoff, or all-modulus Eisenstein mechanism is inferred.

## Fourth checkpoint: assembled replay and final historical repair

```bash
# From this packet directory:
sha256sum -c SHA256SUMS
python3 replay_checkpoint4.py --output /tmp/session139-checkpoint4-replay.json
python3 verify_1475_repair.py
python3 -O verify_1475_repair.py
```

The fresh, wrong-working-directory relocation runs ten scientific kernels
and two manifest verifiers under normal and optimized Python: **44 executions,
20 scientific positives, four source-integrity positives,20 expected
rejections**. Positive pairs match the frozen receipts. Mutations test strict
instead of non-strict tensor/coset inequalities, arithmetic and conditional-
factor constants, the padding count, the140 endpoint, the perfect-matching
indicator, and same-size source corruption. All checks are explicit runtime
conditions, not removable Python assertions. All fourth-checkpoint kernels
use the standard library. [Receipt](REPLAY_CHECKPOINT4.json).

One further explicitly authorized naming correction changes
[THM-1475](../../../01-canon/theorems/THM-1475-the-pfaffian-is-the-odd-function-and-its-spectrum-has-gaps.md):
the order-eight extremizer is a skew-Hadamard sign matrix obtained by adjoining
a source to Paley7, not an impossible even-order regular tournament. A separate
independent exact reconstruction verifies `SS^T=7I`, Pfaffian `-49` and score
sequence `7,3,3,3,3,3,3,3`, normally and with optimization; a false `-48`
Pfaffian is rejected in both modes ([four-run replay](HISTORICAL_1475_REPLAY.json)). The original census
is unchanged. [Before/after pins](HISTORICAL_1475_REPAIR.json) · [Independent receipt](HISTORICAL_1475_INDEPENDENT.json).

## Third checkpoint: cocycle, second moment and regular geometry

[Cocycle proof](cocycle/PROOF.md) · [Integral audit](independent/INTEGRAL_COCYCLE_AUDIT.md) · [Weighted audit](independent/WEIGHTED_COCYCLE_AUDIT.md)

**ACCEPTED.** Integer signed conservation precedes modular reduction. The
support floor n-3 is pointwise for a nonzero binary cochain; it does not
apply to the possibly cancelling mean after weighted rounding. The square
identity is binary; weighted formulas explicitly use expectations of squares.
Disjoint-edge M22 is a necessary second-moment sidecar in the displayed
cospectral same-mean hostile. No new stability constant follows from it.

The integral audit exhausts64 labeled four-vertex and1024 five-vertex
configurations,1024 switchings,1536 relabelings and32768 modular boundary
controls. The weighted audit verifies the polynomial M22 difference and3136
exact binary roundings, including an unstructured nonzero-mean four-vertex
case. It imports only the frozen independently written bridge/integral helpers,
not the producer's implementation. The unsigned boundary and wrong second-
moment constant are rejected normally and with optimization.

[Regular proof](regular_spectral/REGULAR_CRITICAL_LINE.md) · [Independent family/slice audit](independent/REGULAR_SPECTRAL_AUDIT.md)

**ACCEPTED.** The general critical line, reciprocal pole circle and finite
functional equation follow from the real skew restriction. The full regular
rank-two3+2 circle has exact H/Delta width t^4/27, but the full isospectral
fibre is not classified. The producer's rational-polynomial checks and the
independent full rational-parameter enumeration agree. Only the all-one
marked resolvent is fixed: the audit explicitly detects individual vertex-
root difference5*t^2/18. Three producer mutations and false independent H
coefficient are rejected. A separate exact-circle coefficient mutation is
also rejected in the assembled replay. No RH/GRH or spectral-sufficiency
impossibility claim is made.

## Third checkpoint: content and weighted/binary context simulation

[Content](all_moduli/CONTEXT_CONTENT.md) · [Simulation](all_moduli/WEIGHTED_BINARY_SIMULATION.md) · [Independent audit](independent/CONTENT_SIMULATION_AUDIT.md)

**ACCEPTED after incorporated scope clarifications.** The finite gcd formula
uses equality of integer ideals, not a claim that signed linear combinations
are actual contexts. The zero state exists exactly at odd moduli; odd-prime
minimum order isp, and the named composite minima15/33/55 are5/11/11.
The general composite construction is only an upper bound. Content33 for
cyclic R11 is exact; Paley11 uses55 dividing its content, not exact equality.

Weighted profile multilinearity extends the homogeneous-context kernel.
Good-prime binary simulation at3,5,7,11 preserves modular contextual H only,
not exact H, order, spectrum or real inequalities. The particular rational
pair cannot reduce at2 or3. Its all-context spectral propagation uses an
orthogonal conjugacy fixing both row and column all-one vectors; a one-sided
similarity is explicitly rejected. The first-channel response difference is
H(Q)/3456, and the mixed reflector rectangle is1/1728.

The independent binary profiles use disjoint odd-cycle collections followed
by a Stirling transform, unlike the producer's path/subset partitions. Weighted
profiles use an independent endpoint/partition DP. Thirty-one literal pointed
contexts (largest substituted order9), full R11 profile, all120 R5 permutations,
Paley11 affine subgroup closure and explicit mod11 positive words are checked.
False content and weighted-profile mutations remain active under `-O`.

The additional [R11/Paley11 bridge](independent/REGULAR_CONTENT_BRIDGE_AUDIT.md)
checks degree5, cyclic triangles55 and the regular critical-line identities for
both graphs. Their fourth traces880/660 explicitly reject cospectrality;
coarse spectral data do not determine arithmetic context content. No exact
Paley11-content value is claimed.

## Third checkpoint: fresh assembled replay and historical repair

```bash
# From this packet directory:
sha256sum -c SHA256SUMS
python3 replay_checkpoint3.py --output /tmp/session139-checkpoint3-replay.json
```

The fresh relocation runs nine kernels in both normal and optimized Python:
**40 executions,18 positive runs and22 expected rejected controls**. All nine
positive output pairs agree exactly and match the frozen receipts. The shared
independent helper and the pinned content receipt are copied explicitly; there
are no producer imports in the independent audits. All new kernels use the
standard library. [Receipt](REPLAY_CHECKPOINT3.json) · [Reviewed pins](SOURCE_PINS_CHECKPOINT3.json).

One explicitly authorized historical file is repaired:
[THM-1455 / why-eleven-spectra](../../../01-canon/theorems/THM-1455-why-eleven-spectra.md).
It now says Paley7 has14 of35 odd foursets, explains absence of k=2 by the
sharp n-3 support floor, and records that odd-order evenness was already
proved. The dated correction links the new proof and retains the finite-census
scope. [Exact before/after hash record](cocycle/CANON_REPAIR.json).
The distinct same-ID sibling and central mistakes ledger are unchanged;
old historical blob pins continue to identify the old source version honestly.

## Second checkpoint: light-mass dual closure at n>=526

[Proof](stability/PROOF_N526.md) · [Independent endpoint/full audit](independent_threshold526/AUDIT.md) · [Second structural audit](independent/N526_STRUCTURAL_AUDIT.md)

**ACCEPTED, including equality.** The proof retains actual light mass U,
uses a nonempty heavy matching, and pays the light-deficit budget with an
exact nonnegative quadratic dual. Endpoint-convexity and monotonicity—not
sampled orders—supply the uniform geometric tail bound. All divisions,
strict signs, parities and matching-coefficient directions were reviewed.
The final exact slack is `81/8416000`.

The producer and independent integer/Fraction certificates match their pinned
receipts in normal and optimized Python. The independent code imports no
producer implementation. A false extension to525 is rejected, as is a changed
dual multiplier `lambda=d/63`, in both modes. These are controls of the proof
certificate, not counterexamples at525. The older748 proof and its tests are
not rewritten to claim526. Accepted Session138 inputs/equality are inherited,
not all re-proved in this release.

## Second checkpoint: all-modulus operational arithmetic

[Proof](all_moduli/PROOF.md) · [Core/orbit audit](independent/ALL_MODULI_AUDIT.md) · [Padding audit](independent/ALL_MODULI_PADDING_AUDIT.md)

**ACCEPTED.** The all-M equivalence uses integer Newton coefficients, not
division by factorials modulo a composite. Exact transitive periods and
transients use prime-power valuation arguments plus CRT. Universal source
padding has sharp index A(M)-1 and period L(M); the individual delta-dependent
onset is a sufficient bound only. The actual image of Psi is a finite family,
not an asserted ring or a realization of all arbitrary functions.

The exact class counts5,49,685,322101 at primes3,5,7,11 combine all-order
upper bounds and complete finite constructor-quotient reachability. The
producer uses parity-tagged half-vectors; the independent orbit audit uses
full p-entry functions. Its raw unary counts are5,48,685,322100, not the
producer's tagged counts6,48,686,322100. Their final contextual counts agree
after identifying zero tags and adjoining a disjoint-support join at5,11.
All-depth unary zero exclusion at5,11 follows from invertibility, not a
bounded search. These results do not extrapolate to every prime/composite.

Independent executable coverage:

- 209 moduli for exact first-repeat clocks via Newton power sums rather than
  the producer's Stirling recurrence; literal order-nine H55 profile.
- 7,776 signed degree-five profile/modulus pairs, plus full-function orbit
  exhaustion at3,5,7,11.
- All1,099 labeled nonempty tournaments through order5 via literal vertex
  permutations and cuts, rather than producer Hamilton-subset/partition DP;
  70,336 actual tournament/modulus comparisons, 89,402 padding comparisons
  over132 moduli, and the actual E4 response117/172 modulo275.
- Producer separately checks519 clock moduli,70,588 period comparisons,
  35,168 Psi comparisons,225 actual operator/profile checks and its full
  constructor quotient orbits.

All five inherited source blob pins match the declared source commit. The
base theorem's proof/verifier/receipts/pins are verbatim. Two producer summary
files only lose references to the excluded later CONTEXT_CONTENT addendum;
[acceptance and curated selector](all_moduli/ACCEPTANCE.md) document that cut.

## Fresh second-checkpoint relocation and failure controls

```bash
# From this packet directory:
sha256sum -c SHA256SUMS
python3 replay_checkpoint2.py --output /tmp/session139-checkpoint2-replay.json
python3 replay.py --output /tmp/session139-checkpoint1-replay.json
```

The second replay copies six standalone kernels to a new temporary directory.
It performs **24 executions:12 positive runs and12 expected rejections**.
Each normal/optimized output pair must agree and match the frozen receipt.
False threshold, dual multiplier, dyadic period, orbit cardinality and actual
E4-response mutations are rejected by explicit runtime checks, including
under `python3 -O`. All second-checkpoint kernels use the standard library.
Source/output hashes are recorded in [the receipt](REPLAY_CHECKPOINT2.json).

The original first replay remains unchanged:22 executions, eight positive
runs and fourteen expected rejections. Its historical polynomial cross-check
requires SymPy; the two stability kernels and symbolic bridge audit use the
standard library. The archived `independent/SHA256SUMS` has original absolute
scratch paths; use the cumulative root manifest for portable validation.

## First checkpoint retained: n>=748, spectral loss and history

The [n748 audit](independent_threshold/AUDIT.md) accepts the original shared
wedge-budget theorem and equality;20 producer and22 independent exact checks
remain recorded. False747 extension and producer h=0 are rejected. The proof's
original review-pending sentence is its drafting snapshot, not current status.

The [bridge audit](independent/AUDIT.md) independently reconstructs the full
bivariate rational conjugacy, determinant and Hamiltonian polynomial with
standard-library arithmetic. The all-one-root conjugacy proves the full
rooted response is preserved. M23 repairs exact H at n5; binary first-blindness
remains6. Auxiliary binary cocycle controls were retained without promoting
the separate later integral/weighted cocycle continuation.

The [historical checks](riemann_five_eleven/verify_history.py) reproduce the
corrected transfer/tetranacci/pentanacci discriminants, Bernoulli denominator
repeat, Paley11 character sums/ranking33 versus35, and a separate Held–Karp
rational witness. Twenty historical source pins and six inherited input hashes
are preserved in [the first pin record](SOURCE_PINS.json). The exploratory
probe is non-load-bearing; no global priority, RH/GRH or spectral-proof
impossibility is inferred.

## First concurrency reconciliation (checkpoint4)

Publication was guarded and stopped when main advanced. The final local tree
is based on `c16226249f1f27a335cf7c84afdeaa90879b9580`; all concurrent unrelated
entries are preserved. Four overlapping summary/router paths merge current140
with the separate276/314/466 proof credits and guarded-FRACTRAN precision work.
Scientific proofs and all44+4 final replay outcomes remain unchanged.
[Reconciliation](RECONCILIATION.md) · [Tree-entry evidence](RECONCILIATION.json).

## Documentation, scope and publication checks

Native `agents/check_docs.py` needs baseline documents absent from this sparse
checkout. Its unmodified same-version checker passes with read-only overlays
on the clean materialized Session137 ancestor, with every intervening path
accounted for. No baseline mutation or hydration is used. The pre-reconciliation receipt
is [checkpoint-four documentation verification](DOC_OVERLAY_VERIFICATION_CHECKPOINT4.json);
the fourth-checkpoint startup was179966/180000 bytes and frontier450/450 lines. New links,
source hashes, frozen prior lanes, conflict markers and explicit staging scope
are independently checked. After concurrent main changes, a [scoped check](DOC_RECONCILIATION_CHECK.json)
revalidates startup budgets and merged links. The full checker is not rerun
against the opaque concurrent mistakes-ledger change; its remote blob is
preserved exactly instead.

Only the explicitly listed THM-1455 and THM-1475 historical repairs change
older canon; the central mistakes ledger and the other1455 slug remain
untouched. Exploratory probes and unsupported extensions are excluded.
The fifth checkpoint inherits that same opaque ledger through parent `dde756b2`
without reading, querying or changing it, and uses
a [bounded documentation check](DOC_VERIFICATION_CHECKPOINT5.json):179964/180000
startup bytes and450/450 frontier lines, not an unreported full-ledger scan.
Prepared local commits are not remote publication: the publishing agent
separately verifies the actual main commit and every owned blob.
