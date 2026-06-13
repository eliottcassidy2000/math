---
source: codex-2026-06-03-S593
status: proved cap lemma + computation + synthesis for the LRC n=14 Cprime residual
tags: [LRC, n14, Cprime, pigeonhole, cap-face, endpoint-cover, THM-398, HYP-2140, tournament-analysis]
---

# The dual pigeonhole cap face

The user prompt was to consider the "dual pigeonhole cap face" while cycling
between exploration, computation, and formalization toward LRC `n=14`.

The concrete reading that survived is this:

- Bprime is the primal cap obstruction: one component of `G(S')` is too long for
  one `v=nw` danger cap.
- The dual cap obstruction sums inside one primary `n`-clock cell: too many
  short components in that cell exceed the total cap capacity there.

## Lemma H

For `S=S' union {v=nw}` and

```
I_r = [r/n,(r+1)/n),
```

the danger set `D_{nw}` has exactly `w` caps in each `I_r`, each of length
`2/(n^2 w)`.  Thus `D_{nw} cap I_r` has total capacity `2/n^2`.

So if

```
mu(G(S') cap I_r) > 2/n^2
```

for some `r`, then `G(S')` is not contained in `D_{nw}` and `S` is loose.  The
quantitative lower bound is

```
mu(G(S)) >= max_r (mu(G(S') cap I_r) - 2/n^2).
```

This has now been folded into THM-398 as Lemma H.

## Computation

`04-computation/lrc_dual_pigeonhole_cap_face_s593.py` implements the exact
component scan and compares the cap criterion against S581's local owner
criteria.

In deterministic multiple-of-`n` samples:

| n | rows | dual cell-cap routes | residual after cap+local routes |
|---|---:|---:|---:|
| 6 | 2492 | 726 | 145 |
| 8 | 2499 | 1676 | 29 |
| 10 | 2500 | 2205 | 8 |
| 12 | 2500 | 2298 | 0 |
| 14 | 2500 | 2460 | 0 |

The cap criterion is not merely duplicate language.  It gives cap-only exits
before any S581 local owner-descent component exit in smaller rows:

```
n=6: 41
n=8: 26
n=10: 19
```

At `n=14`, it is dominant in the deterministic sample: `2460/2500` rows are
cell-cap overloads, and the remaining `40` are local owner-descent exits.

## Rebase integration: S595

The later monad-compute S595 push strengthens the background Cprime picture:
an exact open-safe-measure audit finds `19600/19600` sampled/systematic
multiple-of-`n` configs loose and `0` tight across controls `n=12,13,14` and
new rows `n=15,16,17,18`, with observed positive-measure margins about
`0.016..0.023`.

This connects directly to the cap face.  Lemma H is not a proof of that whole
S595 phenomenon, but it is an aggregate certificate for the same multiple
branch.  The under-cap residual should therefore be tracked as the common
owner-congruence/large-residue residual across the `n=14..18` frontier, not as a
quirk of the fourteen-runner row alone.

## Named calibration rows

- `unit_shift_AP_n14 = {2,...,14}` is loose by dual cell-cap overload:
  cell `0` has surplus `5/196`.
- `apex_AP_replace_13_n14 = {1,...,12,14}` routes through Lemma C, with no cell
  surplus.
- `near_AP_double_14_n14 = {1,...,11,13,28}` also routes through Lemma C, with no
  cell surplus.

So the cell-cap face separates an aggregate overload decoy (`unit_shift_AP`) from
near-AP multiple rows whose looseness is owner-congruence rather than mass
capacity.

## Tournament Analysis

S593 uses two vertex choices.

For the cap lemma, vertices are the `n` primary clock cells, not runners,
components, or cap centers.  The pairwise observable is `(cell safe load, clock
order)`, the switch points toward larger load/surplus, and the tie Hamiltonian
path is decreasing load with clock order as tie-break.  On the representative
`n=14` row `(4,5,6,7,8,9,10,11,12,13,14,17,18)`, this tournament is transitive:
score histogram `{0:1,...,13:1}`, no directed 3-cycles, singleton SCCs, and one
Hamiltonian path.

For proof routing, vertices are proof criteria.  The cap route sits beside
Bprime and the owner-descent criteria as an independent aggregate exit.

Assumption challenged: the cover proof does not need one component to be too
long.  A cell can fail by having too many individually short components.

## Handoff

The remaining proof shape is now layered:

1. use Lemma H to discharge aggregate cap overloads;
2. use origin-bisection upper/lower cap overloads when small-owner congruences
   pin components to cap centers;
3. on rows with every cell under capacity, use endpoint-owner congruences plus
   the exact `Phi` ramp functional to force a positive gap.

This looks like the right companion to HYP-2105/HYP-2110: owner descent handles
local pins, while the dual cap face handles aggregate mass before the full
large-owner CRT residual is needed.

The S598 follow-up makes the middle line precise.  A small-left endpoint pinned
to a `v=nw` center cannot use the lower half of that danger cap, so a whole
primary cell has only `1/n^2` upper-half capacity available for such components;
right endpoints give the lower-half dual.  In the S598 deterministic samples,
this one-sided certificate routes `23/5000` n=14 rows after the total-cell cap
test fails, and it certifies `apex_AP_replace_13_n14` by an upper surplus
`1/1176`.

Artifacts:

- `04-computation/lrc_dual_pigeonhole_cap_face_s593.py`
- `05-knowledge/results/lrc_dual_pigeonhole_cap_face_s593.out`
- `04-computation/lrc_origin_bisection_upper_caps_s598.py`
- `05-knowledge/results/lrc_origin_bisection_upper_caps_s598.out`
- `05-knowledge/hypotheses/HYP-2140-lrc-dual-pigeonhole-cap-face.md`
