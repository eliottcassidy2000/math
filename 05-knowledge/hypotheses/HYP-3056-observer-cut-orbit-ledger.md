---
id: HYP-3056
title: Observer-cut orbit ledger for controlled-forgetting payloads
status: SYNTHESIS / formal addendum to HYP-3054; not a proof
source: codex-2026-06-26-S220
tangent: T1138
related:
  - HYP-3054
  - HYP-3055
  - HYP-3053
  - HYP-3052
  - HYP-3051
  - HYP-3050
  - HYP-3049
  - HYP-3048
  - HYP-3047
  - HYP-3043
  - HYP-3039
  - HYP-3037
  - HYP-3034
  - HYP-3032
  - HYP-3031
  - HYP-3024
  - HYP-3018
  - HYP-2997
  - HYP-2995
  - HYP-2963
  - THM-572
  - OPEN-Q-108
---

# HYP-3056: Observer-Cut Orbit Ledger For Controlled-Forgetting Payloads

## Claim

HYP-3054 says that every proof quotient has an observer-extension/cut
payload.  HYP-3055 gives the concrete duodecimal first-defect instance.  The
next exact object is not a scalar payload but an orbit:

```text
C_q(x,o) = orbit_Aut_q(x)(boundary(o,x), incidence(o,boundary), q_E(E_o x)).
```

Here `q:X->Y` is the visible quotient, `o` is the next observer or proof
operation, `E_o(x)` is the extended object, and `Aut_q(x)` is the automorphism
group of the currently visible fiber.  The controlled-forgetting question is
whether `C_q` is invisible because it is safe, or invisible because the
quotient has thrown away the coordinate that decides the next proof step.

In short:

```text
controlled forgetting = quotient + next observer + cut-payload orbit + discharge
```

For LRC14, the necessary condition is that every route/status-changing pair in
one coarse fiber is separated, reconstructed, exact, dual-annihilated,
descended, boundary-stopped, or named as residual debt by this orbit ledger.

## Formal Addendum To HYP-3054

Let:

```text
q : X -> Y
P : X -> {boundary, strict-open, route labels, proof exits}
Obs_q(x) = admissible observers over x
```

Examples of observers are new tournament vertices, directed edges, endpoint
owners, deleted runners, danger-cover endpoints, primitive denominator clocks,
Haar squares, analytic smoothing kernels, automaton refinements, residual
capacitor cuts, proof pages, and matrix observability columns.

For an observer `o`, write:

```text
B_o(x)       = the old boundary/fiber slice seen by o
I_o(x)       = incidence or comparison data between o and B_o(x)
q_E(E_o(x))  = the visible shadow after extension
G_q(x)       = Aut_q(x), the visible-fiber symmetry group
```

The payload is the orbit:

```text
C_q(x,o) = [B_o(x), I_o(x), q_E(E_o(x))] / G_q(x).
```

A forgotten observer-cut defect is a pair:

```text
q(x)=q(x'), but C_q(x,o) != C_q(x',o')
```

which changes boundary/open status, route, owner current, topology, certificate
availability, or residual name after the next operation.

The quotient `q` is controlled-forgetting safe only if every such defect has a
recorded discharge:

```text
fiber_constant
reconstructed_from_sidecars
exact_coboundary_or_cycle_boundary
dual_annihilated
descended_familywise
boundary_stopped_at_AP_GW
named_residual_debt
```

The orbit language matters because the missing coordinate is usually not a
raw label.  It is a label modulo the symmetries that the quotient still
respects.  This is the bridge between A000568 orbit counts, HYP-3048 sidecar
observability matrices, HYP-3039 hidden-coordinate ledgers, and HYP-2995
cocycle records.

## Ledger Schema

Every future HYP-2963-style packet quotient should be able to emit rows of:

```text
base_quotient
fiber_id
observer_kind
extension_object
visible_automorphism_group
cut_boundary_slice
cut_payload_word
cut_payload_orbit_id
destroyed_coordinate
changed_lrc_predicate
separating_sidecar
discharge_mode
residual_debt_name
audit_source
```

The `changed_lrc_predicate` field is the guardrail.  If two packets have
different payload orbits but no LRC predicate changes, the cut is harmless for
that proof.  If a predicate changes, the sidecar or discharge must be named
before the quotient can be reused.

## Cross-Problem Ledger Entries

| Thread | Coarse quotient | Observer | Orbit payload | Necessary condition |
|---|---|---|---|---|
| A000568 node perspectives | rooted node cache | new vertex / ordered pair | incident word orbit, sector deck, cross-sector orientation | The first `R(5)=48` to `U(6)=56` defect is not deeper node memory; it is observer-cut orbit data. |
| AP/Goddyn-Wong | residue/automatic/status shadow | endpoint completion and owner deletion | closed-arc `H1` owner support, owner-deletion persistence | Any zero-open non-AP/GW packet must either carry the same owner-essential boundary orbit or name THM-572/F7 debt. |
| q=23 drop/add square | exact `M=2/23` diagonal | off-diagonal Haar swap | exact-M zeta plus endpoint-owner strip orbit | Equal exact magnitude is unsafe until the two-coordinate cut and owner current are accounted for. |
| Residual capacitors | protected status fiber | two-route proof obligation | `residual_capacitor_id`, first cut, zeta/owner/topology exit | Mixed route pairs are finite cut obligations, not anonymous residual mass. |
| K33/state lift | coarse topology or exact scale | proof-page route lift | state-lift id, cross-handoff exit, exact-M/topology cut | K33 is a route-lift payload, distinct from AP/GW boundary homology. |
| Pair-good decoys | false-switch count | blocker generator | blocker tooth, active owner, barcode/normal-fan relation | Decoy count is secondary to the generator orbit. |
| Automata/lacunary shadows | finite word, carry state, residue | exact LRC packet observer | magnitude cocycle, safe-component stalk, endpoint/topology handoff | Finite-state quotients are telemetry until the exact-packet orbit is retained or discharged. |
| Analytic sieve clocks | `mu`, `mu/n`, `mu^2/phi`, smoothing scalar | blindness observer | squarefree blindness report, primitive-period deck, endpoint-owner repair | Analytic capacity meters need a packet-family blindness certificate. |
| Diagonal layers | raw `k^2+k` line observations | layer bridge | rectangle/hourglass cycle orbit, line-potential word | Line counts descend only after cycle residues vanish or become named sidecars. |
| Matrix atlas | spectrum/rank/norm | sidecar column or Schur probe | observability column separating route/status pairs | Matrix invariants are scouts until the column orbit separates proof-critical fibers. |

## Necessary Conditions For A Proof Quotient

1. **Orbit invariance.**  The proposed payload must be invariant under
   `Aut_q(x)`.  If it depends on a destroyed label, it is not yet a quotient
   field.
2. **Observer completeness.**  The next operation has to be named.  A quotient
   can be safe for endpoint completion and unsafe for route handoff.
3. **Predicate relevance.**  A payload is mandatory only when changing it can
   change boundary/open status, route, owner current, topology, certificate
   availability, or residual name.
4. **Discharge explicitness.**  "Probably redundant" is not a discharge mode.
   The row must say reconstruction, exactness, dual annihilation, descent,
   boundary stop, or residual debt.
5. **Tournamentizable comparison.**  If several payloads compete, the pairwise
   relation must be declared: which one separates more proof-critical pairs,
   which costs less, and what tie Hamiltonian path is used.
6. **No scalar replacement.**  Counts such as decoy mass, raw line count,
   exact `M`, automatic word, or spectrum can enter only after their forgotten
   observer-cut orbit has been audited.

## Tournament Analysis

Vertices should be ledger columns and discharge mechanisms, not runners:

```text
cut_payload_orbit
observability_column
cross_sector_orientation
owner_current_orbit
closed_H1_owner_support
residual_capacitor_cut
rectangle_hourglass_residue
primitive_period_deck
analytic_blindness_report
automaton_exact_packet_handoff
raw_scalar_shadow
named_residual_debt
```

Pairwise observable:

```text
number of route/status-changing coarse-fiber pairs separated,
boundary/open preservation,
availability of exact/coboundary discharge,
availability of dual annihilation,
familywise descent,
proof cost,
residual debt introduced.
```

Switch:

```text
A -> B
```

when `A` separates every proof-critical pair separated by `B`, separates at
least one more, or has an earlier certified discharge at equal separation.
The tie Hamiltonian path is:

```text
cut_payload_orbit >
observability_column >
cross_sector_orientation >
owner_current_orbit >
closed_H1_owner_support >
residual_capacitor_cut >
rectangle_hourglass_residue >
primitive_period_deck >
analytic_blindness_report >
automaton_exact_packet_handoff >
raw_scalar_shadow >
named_residual_debt
```

Directed cycles in this tournament are not failures.  They are warnings that
two discharges do not commute, usually meaning the proof needs a bicomplex or
fiber product rather than a linear sidecar order.  That is the abstract form of
the Haar rectangle/cocycle lesson.

## Breakthrough Target

The missing bridge from controlled-forgetting philosophy to proof machinery is
an executable observer-cut ledger:

```text
coarse fiber pair
  -> all admissible next observers
  -> payload orbit under visible automorphisms
  -> sidecar column
  -> discharge mode
```

If the ledger has no mixed boundary/open defects and every mixed route defect
has a certified discharge, then the quotient becomes a theorem carrier rather
than a research metaphor.  This is the concrete way to combine HYP-3054's
observer-extension calculus with HYP-3039's hidden-coordinate ladder.
