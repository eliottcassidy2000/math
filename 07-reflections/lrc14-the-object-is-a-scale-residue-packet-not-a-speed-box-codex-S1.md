# LRC(14): the object is a scale-residue packet, not a speed box

**Session:** codex-2026-07-14-S1  
**Status:** frontier audit plus exact correction; LRC(14) remains open in this repository.

## 1. The honest frontier

For a 13-set of distinct positive integer speeds `S`, write

`M(S)=max_t min_{v in S} ||vt||`.

The target is `M(S)>=1/14`.  The strongest stable chain at this checkout is:

1. The non-covering sieve reduces the new content to divisor-complete sets,
   conditional on the repository's imported `LRC(<=13)` result (THM-366/523).
2. THM-738 proves every family with at least ten speeds in `{1,...,14}` is
   lonely.  Consequently the far-count class `#{v>14}<=3` is closed.
3. THM-731/732 turn a chosen peel into an exact autocorrelation discrepancy,
   expressible by rational Bernoulli `B2` endpoint-pair sums.
4. THM-755 proves the capped Fourier envelope: for a fixed core `P`, every
   peeled speed `v>r_P/(pi|G'_P|)` certifies positive lonely measure.
5. Named coherent packs, near-AP bodies, aligned carriers, and hierarchical
   clusters have additional exact routes.

What is not proved is completeness of these routes on the `f>=4` covering
class.  The latest raw finite-band claim does not survive scale covariance.

## 2. How the object sharpened over time

The repository's long path is not a list of unrelated analogies.  Each phase
identified data that the previous quotient had forgotten.

### Phase A: rational witnesses and divisibility

The first durable reduction was arithmetic.  If a modulus `q<=14` is missed,
the time `1/q` gives an explicit route through a lower-runner instance.  This
isolated covering/divisor-complete sets as the only new class.

What it preserved: an explicit witness and the exact threshold.  
What it destroyed: where witnesses lie inside the remaining covering class.

### Phase B: the torus line and exact circle arrangement

The central-box view changed the problem from runner folklore to intersection
of a one-parameter torus orbit with a closed box.  Cutting the circle at all
danger-arc endpoints made the safe set an exact rational interval complex.
Endpoint owners, wall crossings, and equality points separated `M>=1/14` from
the stronger statement that the safe set has positive measure.

What it preserved: all metric and boundary information.  
What it destroyed: almost nothing, but the raw arrangement grows with scale.

### Phase C: gaps, Farey cells, and Ostrowski/three-distance structure

Gap and rational-grid work explained why pair sums and differences are the
candidate denominators of exact maxima.  Farey and continued-fraction cells
organize when endpoint order changes.  The deep well `{1,...,12,182}` and the
`14/183` covering-floor candidate arose naturally here.

What it preserved: exact event order and denominator geometry.  
What it destroyed: a global accounting of many simultaneous blockers.

### Phase D: covering depth, Bonferroni trees, and exact finite bodies

Coverage multiplicity turned union-of-danger-arcs geometry into factorial
moments and exact Bonferroni trees.  THM-738 is the strongest reproducible
finite theorem in the current endgame: 1001 ten-element bodies and 4,677,712
exact bottom checks.  This is why `f<=3` is genuinely closed.

What it preserved: exact interval widths and covering obligations.  
What it destroyed: scale coherence outside the fixed small body.

THM-741 would extend this to nine-element bodies and hence `f<=4`, but its
metadata is still `CLAIMED`; its 2002-body run and output are absent.

### Phase E: Fourier, relation lattices, and singular series

Fourier expansion exposed additive relations and explained AP extremality, but
absolute coefficient bounds diverged at the threshold.  The useful lesson was
not "Fourier fails".  It was that the correct spectral object must retain signs,
localization, and the exact interval endpoints.  Unsigned energy and truncated
Bonferroni bounds cannot see the required cancellation.

What it preserved: arithmetic resonance across all scales.  
What it destroyed when used naively: location, signs, and boundary ownership.

### Phase F: autocorrelation, Bernoulli `B2`, and the capped envelope

THM-731 kept the whole leave-one-out good set intact and bounded its covariance
with one runner by positive autocorrelation discrepancy.  THM-732 converted
that discrepancy to an exact `B2` edge-pair sum.  THM-755 then spliced two
views of every Fourier coefficient: the origin cap `|c_m|<=|G'|` and the spoke
envelope `|c_m|<=r/(pi m)`.

What it preserved: exact metric slack and endpoint phases.  
What it destroyed: nothing per fixed core; the mistake was treating a
per-core cutoff as a global raw-speed cutoff.

### Phase G: safe peels and recursive reduction

A runner safe at a complement's tight point can be removed without changing
`M`.  This reveals that much of the covering class is a lower-runner problem in
disguise.  The exact lemma is sound; the claim that every irreducible residual
lands in an existing tile remains empirical.

What it preserved: exact `M` along a reduction chain.  
What it destroyed: a canonical choice of tight point and a proof that all
irreducibles have a listed normal form.

### Phase H: coherent packs, affine clusters, and hierarchical clocks

Pack-clock and cluster-clock work showed directly that absolute diameter is the
wrong complexity parameter.  Large speed sets may be a dilation or translation
of a tiny shape.  Their arithmetic is controlled by one shared clock plus a
small number of detuned offsets.

What it preserved: scale coherence and exact residue fibers.  
What it destroyed if compressed too far: primitivity and divisor covering.

This phase was the clearest warning against the later raw finite-band claim.

### Phase I: tournaments, observer changes, and leader handoffs

Raw tournaments on runners were useful telemetry for parity, symmetry, and
wall changes, but no unweighted orientation determines lonely clearance.  The
leader ledger sharpened the natural vertices from runners to wall-handoff
events and pair-sum obligations.  Tournament data becomes faithful only with
metric sidecars such as endpoint signs, exact phases, and `B2` weights.

What it preserved: ordering, cycles, switch sensitivity, and observer changes.  
What it destroyed without sidecars: distances, equality, and scale.

### Phase J: topology, sheaves, normal fans, and controlled forgetting

Cech, barcode, oriented-matroid, cocycle, and sheaf language did not produce a
standalone proof.  It did produce the correct audit rule: a quotient is safe
only if every forgotten label is fiber-constant, reconstructed, a coboundary,
dual-annihilated, or sent to a named residual.

What it preserved: compatibility and handoff obligations.  
What it destroyed when left unlabelled: the numerical LRC predicate itself.

### Phase K: formalization and executable certificates

Lean now checks substantial algebraic, spectral, Raabe, and grid-deficit
pieces.  Exact Python interval and Bernoulli engines produce reviewable
witnesses.  This exposed a second important distinction: a theorem schema can
be formally sound while the finite coverage claim that invokes it is sampled.
The current formal composition also retains `LRCUpTo13` as an external input.

## 3. The scale correction

For any core `P` with positive-measure good set and integer `c>=1`, the
degree-`c` circle cover gives

`G'_{cP}=T_c^{-1}(G'_P)`,

so

`|G'_{cP}|=|G'_P|`, `r_{cP}=c r_P`, and `v*(cP)=c v*(P)`.

This is HYP-6780.  A finite residual interval for each fixed core is therefore
not a uniformly finite collection of integer speed sets.

The family

`V_c={c,2c,...,12c,13c+1}`

is primitive.  Under `13|c`, it is covering exactly when
`gcd(c,14)>1` or `c=1 mod 14`.  Its first covering scale is `c=26`.
For every such `c`, it lies below the top-peel cutoff and has exact
`M(V_c)=1/13`.

Three current claims are thereby corrected:

- the global terminal cutoff near `500` was sampled, not proved;
- `f>=4` does not imply `M>=0.097` (these examples have `f=13` and `M=1/13`);
- a covering near-dilate does not require `32760|c` and can occur inside the
  reported `(220,475]` band.

None of this refutes LRC(14).  It refutes the proposed compactness variable.

## 4. The underlying object we were missing

The live object is an arithmetic bundle over projective packet shapes:

`normalized cluster tree`
`+ coherent scale parameters`
`+ divisor-covering residue fibers`
`+ detuned offsets`
`+ exact endpoint/wall complex`
`+ certificate or blocker sidecar`.

The base records scale-normal geometry.  The fiber records precisely the data
that covering and primitivity need.  The endpoint complex supplies metric
clearance.  The certificate sidecar records which route actually proves the
packet lonely.

This reconciles several views that previously seemed to compete:

- the torus arrangement is the faithful geometric total space;
- CRT and p-adic data describe arithmetic fibers;
- pack/cluster trees describe the projective base;
- `B2` autocorrelation measures metric interaction of boundary events;
- the endogenous pair-sum blocker complex (HYP-6785) is a finite per-family
  incidence presentation;
- safe peels and certificate handoffs are morphisms between packets.

## 5. Recursive research program

The correct recursion should act on a normalized cluster tree, not on
`max(S)` alone.

1. Extract a maximal coherent packet `cH` or translated packet `a+cH`.
2. Retain offsets and residues modulo the divisor obligations `2,...,14`.
3. Dispatch a large packet with few detuned speeds by THM-668/737.
4. Dispatch translated or two-level clusters by THM-739/740.
5. On the incoherent remainder, choose the peel maximizing exact THM-731 slack
   and apply THM-755/732.
6. If no certificate fires, refine the cluster tree at its largest normalized
   gap or expose the first nonzero blocker/cocycle class.
7. Prove that the chosen complexity decreases: fewer packets, smaller offset
   alphabet, fewer uncovered residue obligations, or smaller endpoint SCC.

This would turn the infinite raw family into finitely many normalized base
types with controlled arithmetic fibers.

## 6. Tournament Analysis after the assumption challenge

Candidate vertex sets considered: runners, gaps, fixed circle sections,
section boundaries, wall events, residues, cover arcs, Fourier modes, matroid
circuits, blocker edges, and proof obligations.

For the scale audit, quotient carriers are the most faithful vertices.  The
pairwise observable compares preservation of the LRC predicate, covering,
primitivity, capped ratio, witness, scale orbit, and finite state.  Predicate-
first and compression-first gauges produce two edge flips but agree that
`shape+residue` is first.  Both tournaments are transitive, have no directed
3-cycles, singleton SCCs, score histogram `{0:1,1:1,2:1,3:1,4:1}`, and one
Hamiltonian path.

For the autocorrelation route, a better future tournament uses signed endpoint
wall events as vertices, phase half-circle orientation as the switch, cyclic
endpoint order as the tie path, and
`sigma_e sigma_f B2({v(x_f-x_e)})` as the mandatory metric sidecar.  The
unweighted tournament preserves interlacing but destroys the discrepancy.

## 7. Ranked next contributions

1. Prove a scale-normal uncapped-family trichotomy: coherent pack, translated/
   hierarchical cluster, or genuinely incoherent compact atlas.
2. Implement and complete the missing resumable THM-741 2002-body exact run.
   This would extend the unconditional finite closure from `f<=3` to `f<=4`.
3. Turn the sampled `q<=25` good-period observation into a theorem on the
   normalized incoherent class, or find its first counterexample.
4. Complete the Lean interval-autocorrelation-to-`B2` single-pair structural
   identity.  The surrounding Fourier/Raabe chain is already formalized.
5. Audit the external `LRC(<=13)` citation and make the formal dependency
   explicit in every claimed top-level assembly.
6. Reconcile THM-724/726/757 status language: the named rigidity lanes and
   exact banks are strong, but the global multi-killer floor is still described
   as conjectural in THM-757.

## 8. Session contribution

- Proved and recorded the scale law and unbounded covering ray in HYP-6780.
- Corrected THM-757's covering condition, minimum scale, band statement, and
  enumeration scope.
- Corrected THM-758's global finite-band and `0.097` claims.
- Added an exact scale audit with quotient-carrier Tournament Analysis.
- Repaired the shared `M_exact` implementation and removed floating-point
  certification from the THM-755 protocol.

The frontier did move, but by removing a false compactness claim and replacing
it with the projective arithmetic object the remaining proof must actually
control.
