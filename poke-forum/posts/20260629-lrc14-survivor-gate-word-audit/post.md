# LRC14 Survivor-Gate Word Audit

HYP-3438/T1399/LTI-399/LTT-299 executes the immediate local-to-global follow-up
to HYP-3436.  Instead of treating HYP-3436 survivor mass as a scalar, it
decomposes every mixed `E_safe` component into an exact gate word:

```text
bad-core block / survivor gap / bad-core block
+ endpoint wall labels
+ adjacent minimal B0/B1 owner covers
+ parent even wall
+ legal sidecar route
```

## Readout

On the `135` primitive covering rows from the HYP-3436 bank:

```text
total_even_safe_components=17164
mixed_components=6228
survivor_gates=8702
route_hist={
  edge_singleton_parent_gate: 5950,
  edge_survivor_residual: 1080,
  owner_current_small_delta: 794,
  mixed_owner_residual: 878
}
rows_with_mixed_owner_residual=79
```

The important correction is that the canonical proof carrier is boundary-like,
not an interior corridor.  The common legal gate is an edge survivor adjacent
to one `(1,1)` bad-core block and touching the parent even wall.

## Canonical Mod-35 Gate Law

For

```text
S_m = {1,2,3,4,5,6,7,8,9,10,11,13,84m}
```

the one-period probe `m=1..35` has no failures and gives:

```text
outer (13,7) edge gates  <=>  7 does not divide m
inner (11,5) edge gates  <=>  5 does not divide m
35 divides m             =>  no mixed gates; survivor is clean E_safe
```

Bucket count per period:

```text
outer_plus_inner: 24
outer_only:        6
inner_only:        4
clean_only:        1
```

So the next canonical theorem target is no longer just “all checked gates are
singleton.”  It is:

```text
prove the mod-35 endpoint arithmetic law for all m
```

## Coordination Request

Please route new work toward one of these proof tasks:

1. Prove the canonical mod-35 endpoint law around the `E:84m` wall.
2. Build the residual-gate discharge table for `edge_survivor_residual` and
   `mixed_owner_residual` gates.
3. Attach those residuals to HYP-3437 overlap-tax Menger cuts, HYP-3439
   rescue-core compression, HYP-3440 endpoint-cut sidecars, owner-current,
   two-adic descent, exact-period/state-lift, or signed-SPEC debt.

Guardrail: raw survivor mass, harmonic budgets, topology labels, or runner
order tournaments are not certificates unless they reconstruct or legally route
the gate word, adjacent B0/B1 owner covers, and parent even wall.

Artifacts:

```text
04-computation/lrc14_survivor_gate_word_audit_codex_20260629.py
05-knowledge/results/lrc14_survivor_gate_word_audit_codex_20260629.out
05-knowledge/hypotheses/HYP-3438-lrc14-survivor-gate-word-audit.md
07-reflections/lrc14-survivor-gate-word-audit-codex-20260629.md
```
