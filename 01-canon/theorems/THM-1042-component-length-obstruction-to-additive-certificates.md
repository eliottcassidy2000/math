---
id: THM-1042
title: The component-length obstruction — G_{1..k}(1/14) has longest component L_max(k) with 1/L_max(k) > k+1 for every k ≥ 3, so an additive certificate can NEVER absorb the next consecutive speed. This closes the additive route on every family with a small-integer core, and explains uniformly why THM-1015 needed large killers.
status: >
  PROVED + VERIFIED-EXACT + CORRECTED + INDEPENDENTLY HOSTILE-AUDITED.
  Positive component lengths and the additive obstruction are exact.  The
  original AP11 table cell `#comps=20` was false: the implemented
  positive-length convention gives 14 arcs, while the full closed safe set
  has 18 components after four isolated equality walls are included.  The
  measure, L_max=1/77, threshold 77, and every method-level consequence
  survive.  This is a NEGATIVE result about a method, not about LRC(14).
source: klein-2026-07-18-S327
audit: >
  CORRECTION AUDIT PASS (MISTAKE-464, 2026-08-23).  The canonical executable
  returns 14 positive arcs.  An import-independent exact wall graph finds
  128 equality walls, 32 closed-safe walls, 14 strict-safe cells, 14
  nondegenerate closed arcs and the four isolated safe points
  3/14,5/14,9/14,11/14, for 18 closed topological components.  It independently
  reproduces measure 10931/194040 and L_max=1/77 in 165 active gates; normal
  and optimized runs byte-match.  No proved downstream theorem consumes the
  false multiplicity 20.
depends_on:
  - THM-1004  # the interval-survival tail whose hypothesis this measures
  - THM-1015  # the clustered closure, whose large-killer hypothesis this explains
related:
  - THM-1026  # opus's five-slot ledger (the alternating route, killed separately in klein-S325)
script: 04-computation/lrc_complength_klein_S327.py
output: 05-knowledge/results/lrc_complength_klein_S327.out
script_sha256: 50628156cb00d61da5ce5326f2e329694a85e7cfda45ba2c8ad1f8b01f082bb9
output_sha256: 9772bfe2ff3285b03bbc83ae46351b34c7add00170d7400f0ec140c7474b0a91
independent_audit_script: 04-computation/lrc_thm1042_ap11_component_correction_independent_audit.py
independent_audit_output: 05-knowledge/results/lrc_thm1042_ap11_component_correction_independent_audit.out
independent_audit_script_sha256: 7279ebbac013d92e92445538073989cde465aff45753f4e378791c3e3731bf93
independent_audit_output_sha256: 5e2899c96cb172b37523e1c5bce1921374bde179d23420f61e537049ffb9f413
independent_audit_semantic_sha256: c82d19fc42d241f1abd00e475ac58462816a0c23f1a35d60f29e3a4177c87b22
hash_basis: raw LF bytes
---

# THM-1042 — the component-length obstruction

The additive certificates (THM-1004's interval-survival tail, and the measure-recursion variant) absorb a
new speed `w` into a good set by charging it a *proportional* loss. That accounting is only valid when the
good set's components are **longer than `w`'s arc period `1/w`**; otherwise a single period spans a whole
component and the component is clipped or swallowed outright, with no proportional regime. So

```text
base B admits speed w   only if   w > 1/L_max(B),        L_max(B) := longest component of G_B(1/14).
```

This theorem computes `L_max` and shows the criterion is never satisfiable by a consecutive speed.

## The exact data

Write `G_B^>` for the strict safe set and `G_B^>=` for the closed safe set.
Every positive-length component is an arc with endpoints
`(14k ± 1)/(14v)`, `v ∈ B`, so its length is a rational with denominator
dividing `14·v·v'` for some `v,v' ∈ B` (e.g. `1/588`,
`588 = 14·6·7`).  Let `r_+(B)` be the number of positive-length arcs.
The closed set can additionally contain isolated equality-wall components;
they have zero measure and do not affect `L_max`.  Computed exactly:

| `B` | `μ = |G_B|` | `r_+` | `L_max` | `1/L_max` | next speed `k+1` |
|---|---|---|---|---|---|
| {1..3} | 0.69048 | 4 | 5/21 | 4.2 | 4 |
| {1..4} | 0.61905 | 6 | 9/56 | 6.2 | 5 |
| {1..5} | 0.50476 | 10 | 4/35 | 8.8 | 6 |
| {1..6} | 0.45714 | 12 | 1/12 | 12.0 | 7 |
| {1..7} | 0.33469 | 18 | 3/49 | 16.3 | 8 |
| {1..8} | 0.26582 | 20 | 5/112 | 22.4 | 9 |
| {1..9} | 0.18107 | 20 | 2/63 | 31.5 | 10 |
| {1..10} | 0.13798 | 20 | 3/140 | 46.7 | 11 |
| {1..11} | 0.05633 | 14 | 1/77 | 77.0 | 12 |

For `B={1,...,11}`, the full closed set `G_B^>=` has 18 topological
components: the 14 nondegenerate arcs counted by `r_+`, together with the
four isolated safe points

```text
3/14, 5/14, 9/14, 11/14.
```

The strict safe set has 14 open components.  The point `3/14` is the minimal
endpoint hostile: speed 5 enters danger immediately to its left and speed 9
immediately to its right.  This is the first failed implication in the old
count, namely that every closed-safe component contains a positive safe wall
cell.

**`1/L_max(k) > k+1` at every row, and the gap widens**: the ratio `(1/L_max)/(k+1)` runs
`1.05, 1.24, 1.47, 1.71, 2.04, 2.49, 3.15, 4.25, 6.42`. `1/L_max` grows superlinearly while the next
available speed grows linearly.

## Consequence

**An additive certificate can never absorb a consecutive speed.** Hence it fails on every family whose
speeds contain a run of small consecutive integers — the AP, the deep well, the GW family, and every
covering family with a small-integer core. Checked against the deep well `{1,…,12,182}`: *every* initial
split is blocked, and the blocking speeds are exactly the consecutive ones.

```text
base {1..7}  needs w > 16.3 ; killers 8,9,10,11,12,182 ; blocked by 8,9,10,11,12
base {1..11} needs w > 77.0 ; killers 12,182           ; blocked by 12
```

## What this explains, uniformly

Three separate failures now have one cause:

1. **THM-1015 required large killers** (thresholds 65.7 … 347.5). Not a convenience — those thresholds are
   `1/L_max` of the respective bases. The clustered stratum closed *because* its killers were large.
2. **The measure-recursion (klein-S326) died at `w = 8`** with boundary `2δN/w = 0.321` exceeding the
   surviving measure. `8 < 16.3 = 1/L_max({1..7})`.
3. **The radius-3 fragmentation wall** (klein-S314): the same short components, seen from the
   Hamming-radius side.

Changing the state variable — largest interval → (measure, positive-length
component count) — removes the `r < 1/(2δ) = 7`
cap but not this, because both formulations price a speed against the component scale.

## Scope, stated plainly

This is a theorem about the **method**, not about LRC(14). It says the additive/measure family of
certificates cannot reach the small-speed regime, and gives the exact threshold at which they can. It does
**not** bound `M` for any family, and nothing here contradicts LRC(14). Combined with klein-S325 (the
alternating truncation `B₅ < 0` on real families, first clear at `B₁₁`) and klein-S324 (no pairwise-only
invariant characterizes tightness), the three certificate families the fleet has been using are now each
delimited by an exact, checkable criterion.

The honest reading: the small-speed / compact regime is not awaiting a sharper constant. It is outside the
reach of proportional-loss accounting, because there the good set has no component longer than the arcs
being added.

Reproduce the corrected table and its independent endpoint audit with

```bash
python3 04-computation/lrc_complength_klein_S327.py
python3 -O 04-computation/lrc_complength_klein_S327.py
python3 04-computation/lrc_thm1042_ap11_component_correction_independent_audit.py
python3 -O 04-computation/lrc_thm1042_ap11_component_correction_independent_audit.py
```

Each normal/optimized pair must byte-match its corresponding frozen output
in `05-knowledge/results/`.
