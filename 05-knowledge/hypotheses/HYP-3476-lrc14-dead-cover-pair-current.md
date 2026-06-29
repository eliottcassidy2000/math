# HYP-3476: LRC14 Dead-Cover Pair-Current Exception Audit

**Status:** RESERVED STUB / computation in progress; not a proof.
**Created:** codex-2026-06-29
**Owner:** codex
**Artifacts claimed:**
- `04-computation/lrc14_dead_cover_pair_current_codex_20260629.py`
- `05-knowledge/results/lrc14_dead_cover_pair_current_codex_20260629.out`
- `07-reflections/lrc14-dead-cover-pair-current-codex-20260629.md`

## Question

HYP-3472 proves on the current `135`-row bank that every dead-component row
has a low-rank E/branch gate touching the dead-cover projection, but a single
gate gives a projection-edge cut only on `123/130` dead rows and a separating
current only on `121/130`.

This hypothesis asks whether the named HYP-3472 exceptions are an artifact of
using one gate at a time.  For the edge-cut exception rows

```text
random_covering_001
random_covering_031
random_covering_039
random_covering_062
random_covering_074
random_covering_086
random_covering_101
```

and the two additional separating-current exceptions

```text
covering_AP_with_84
ap_omit_12_tail_84x01
```

test whether a pair of low-rank E/branch gates, or a mirror-orbit-aware pair,
removes enough adjacent blocker labels from the HYP-3451 dead-cover projection
to create a projection-edge cut or a separating branch current.

## Planned Method

Reuse the HYP-3472 projection and gate-adjacency code.  For each target row,
enumerate unordered pairs of low-rank E/branch survivor gates.  For each pair,
remove the union of adjacent blocker labels from the dead-cover projection and
measure:

- remaining edge count and component count;
- largest-component drop;
- whether the pair is a projection-edge cut;
- whether the resulting components are branch-separated;
- whether the pair is the exact mirror partner of HYP-3475 or belongs to the
  same mirror-orbit family.

## Assumption Challenge

Tournament vertices will be proof obligations and finite current carriers, not
runners or arcs.  Candidate vertex sets considered before fixing this quotient:
runners, gaps, fixed circle sections, section boundaries, wall-crossing
events, residues, cover arcs, Fourier modes, matroid circuits, survivor gates,
mirror orbits, gate pairs, blocker labels, graph cuts, and formal proof
obligations.

The quotient preserves the predicate:

```text
dead-cover obstruction on a target row has a <=2-gate boundary-current carrier.
```

It destroys individual interval geometry, ordered branch orientation inside a
mirror pair, and full row identity unless the script retains typed sidecars.
HYP-3474 is the guardrail: count or color profiles are telemetry unless they
reconstruct the projection/cut predicate.

## Expected Proof Pull

If all HYP-3472 exceptions close under two gates, HYP-3473 can expose a smaller
terminal packet: universal touch, then at most two adjacent gate labels give a
projection-current certificate, with AP84 and random rows handled by finite
pair clauses.  If not, the remaining exception rows become a sharper debt list
for HYP-3455, HYP-3451 conductance, HYP-3460 phase-branch bypass, two-adic
descent, owner-current, signed-SPEC/Rprime, or state-lift routes.

## Links

HYP-3476 depends on HYP-3475, HYP-3474, HYP-3473, HYP-3472, HYP-3471,
HYP-3462, HYP-3470, HYP-3461, HYP-3460, HYP-3459, HYP-3458, HYP-3455,
HYP-3453, HYP-3451, HYP-3450, HYP-3438, HYP-3436, HYP-2595, THM-523,
T1436, LTI-436, LTT-336, and OPEN-Q-108.
