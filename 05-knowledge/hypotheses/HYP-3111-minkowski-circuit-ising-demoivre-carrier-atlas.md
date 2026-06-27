---
id: HYP-3111
title: LRC14 Minkowski, circuit-complexity, Ising, and De Moivre carrier atlas
status: EVIDENCE / executable carrier scout; not a proof
source: codex-2026-06-27-S264
tangent: T1187
technique: LTI-248
tournament_technique: LTT-146
related:
  - HYP-3110
  - HYP-3109
  - HYP-3108
  - HYP-3107
  - HYP-3106
  - HYP-3103
  - HYP-3097
  - HYP-3096
  - HYP-3089
  - HYP-3088
  - THM-575
  - THM-573
  - OPEN-Q-108
external_sources:
  - https://en.wikipedia.org/wiki/Minkowski%27s_theorem
  - https://en.wikipedia.org/wiki/Circuit_complexity
  - https://en.wikipedia.org/wiki/Lee-Yang_theory
  - https://en.wikipedia.org/wiki/Ising_model
  - https://en.wikipedia.org/wiki/Quintic_function#De_Moivre's_quintic
---

# HYP-3111: LRC14 Minkowski, Circuit-Complexity, Ising, And De Moivre Carrier Atlas

## Claim

This lane merges the user's requested Minkowski theorem, circuit complexity,
Ising model, and De Moivre quintic prompts into the live HYP-3107--HYP-3110
proof frontier.

It is not a proof.  The intended contribution is a controlled-forgetting
interface:

```text
minkowski_lattice_body
  -> convex symmetric finite-packet forcing sidecar;
proof_circuit_complexity
  -> size/depth/uniformity ledger for proof-state gates;
ising_partition_zero
  -> finite spin partition-function zero geometry refining Lee-Yang signals;
demoivre_quintic_fold
  -> exact finite-depth algebraic cancellation/fold witness;
lee_yang_root_curve
  -> HYP-3109 whole-PGF-root sidecar;
observer_gluing_certificate
  -> HYP-3107/HYP-3097 terminal proof interface;
finite_address_packet
  -> CRT/dyadic/level-7 address carrier for the residual row.
```

The preserved LRC predicate is: a residual primitive row must retain enough
coverage, CRT/dyadic lift, observer-gluing, root-locus, endpoint-owner, and
finite-address data to imply `LRC14Statement` through the HYP-3107 frontier.
The destroyed coordinates are raw runner labels, raw time, branch choices in
algebraic folds, and any scalar-only root/moment/lattice/circuit count without
its sidecar.

## Remaining Work

1. Replace the toy finite Ising carrier graph with the actual HYP-3109
   zero-real one-swap component and compare partition-zero angles before and
   after crossing the `#real=0 -> #real=2` wall.
2. Search for an LRC-native q-body whose inequalities are component count,
   segment clearance, endpoint-owner debt, and finite-address status rather
   than arbitrary Euclidean radius.
3. Turn the De Moivre Laurent identity into a formal lemma and test whether
   HYP-3110's branch/orbifold packets need a fifth-root orbit sidecar.
4. Use the proof-circuit ledger as a regression check for future shortcuts:
   any shortcut must supply `q_witness_gate`, or both `level7_sieve` and
   `dyadic_lift`, or route through `observer_gluing_certificate` plus
   `finite_address_packet`.

## Assumption Challenge

Do not assume the useful vertices are runners, arcs, roots, graphs, or algebra
systems.  The candidate vertices are proof obligations and carriers:
Minkowski convex bodies, circuit gates, Ising spin packets, De Moivre algebraic
folds, Lee-Yang root curves, observer-gluing certificates, and finite-address
packets.  The quotient preserves progress toward a finite LRC14 certificate;
it destroys raw labels and branch coordinates unless a sidecar explicitly
stores them.

## S264 Scout Readout

Artifacts:

- `04-computation/lrc_minkowski_circuit_ising_demoivre_carrier_codex_s264.py`
- `05-knowledge/results/lrc_minkowski_circuit_ising_demoivre_carrier_codex_s264.out`
- `07-reflections/minkowski-circuit-ising-demoivre-carrier-atlas-codex-s264.md`

Main measurements:

1. **Minkowski q-lattice body.**  The named HYP-3108/HYP-3109 packets span
   affine q-rank `6`.  The selected independent q-difference basis has
   covolume proxy `6.795578624e-12`.  For the Euclidean ball in this generated
   lattice, the Minkowski threshold radius is `0.020934`, just below the
   shortest nonzero named q-difference `E_star_k12 <-> gw12` with norm
   `0.021715`; the ratio `vol(B_R)/(2^rank*covolume)` at that shortest named
   radius is `1.245847`.
2. **Proof-state circuit complexity.**  A monotone proof-frontier circuit with
   ten input obligations has gate-size `8`, depth `4`, maximum fan-in `3`, and
   all ten inputs essential.  Minimal certificates are:
   `q_witness_gate`; `level7_sieve AND dyadic_lift`;
   `demoivre_fold AND observer_gluing_certificate AND finite_address_packet`;
   and two root/ear observer exits of size `5`.
3. **Ising/Lee-Yang finite packets.**  Three ferromagnetic toy packets
   (`path6`, `cycle6`, and the seven-carrier proof graph) have partition zeros
   numerically on the unit circle with maximum radial error about `1.11e-16`.
   This supports the HYP-3109 discipline: retain the whole zero set, not a
   single scalar moment.
4. **De Moivre quintic fold.**  The script verifies exactly, as a Laurent
   identity,
   `x = u - a/u` implies
   `x^5 + 5*a*x^3 + 5*a^2*x = u^5 - a^5*u^-5`.
   If `y=u^5` and `y^2+b*y-a^5=0`, this gives
   `x^5 + 5*a*x^3 + 5*a^2*x + b = 0`.  This is the cleanest
   Lean-facing algebraic fold in the session, but it needs a finite-address
   sidecar to touch LRC14.
5. **Tournament Analysis.**  Vertices are proof carriers/obligations, not
   runners or arcs.  The sidecar tournament has score histogram
   `{0:1,1:1,2:1,3:1,4:1,5:1,6:1,7:1,8:1}`, no directed 3-cycles, singleton
   SCCs, one Hamiltonian path, and priority path:
   `finite_address_packet -> observer_gluing_certificate ->
   proof_circuit_complexity -> zero_real_ear_map -> lee_yang_root_curve ->
   demoivre_quintic_fold -> minkowski_lattice_body -> ising_partition_zero ->
   raw_scalar_p0`.

## Proof-Frontier Consequence

HYP-3111 does not close LRC14.  It sharpens the frontier into the following
proof obligation:

```text
residual row
  -> finite_address_packet
  -> observer_gluing_certificate
  -> one of:
       root/ear/Lee-Yang sidecar with retained zero geometry,
       algebraic De Moivre-style finite fold plus address branch,
       q-lattice convex-body sidecar whose body is declared and matched
          to a legal LRC predicate.
```

The circuit readout says the route is shallow as a formal DAG, but not cheap:
the observer exit still needs both the finite-address packet and a carrier
certificate.  The Minkowski and Ising imports are therefore not substitutes for
the witness route; they are gauges for whether a quotient has retained enough
coordinate data to make observer gluing legal.

## Incoming S262 Root-Lattice Integration

Post-rebase, S262 added a bounded-bank root-lattice reachability supplement to
HYP-3108.  It strengthens the constraint on this hypothesis: in the
`{0}+7` anchored bank, high `p0` correlates positively with nearest-root radius
and residue entropy, but negatively with Bravais peak:

```text
corr(p0,nearest_root)    = +0.899
corr(p0,real_roots)      = -0.483
corr(p0,Bravais_peak)    = -0.430
corr(p0,residue_entropy) = +0.541
```

It also finds strict one-swap descent traps at the root-lattice level, with
Savitch midpoint depth `3` for the reachable cases.  Therefore HYP-3111 should
not use Minkowski as a "large Bragg peak" crystallinity heuristic.  The
Minkowski body must be reciprocal-flat and proof-native: inequalities should
encode root stratum, segment clearance, entropy/flatness, finite-address
status, and observer debt rather than raw lattice concentration.
