---
source: oracle-2026-06-01-S553o
status: progress (a provable LRC-adjacent theorem: "at most one runner is near"; honest n=14 counterexample-impossibility locus)
tags:
  - lonely-runner
  - first-moment
  - almost-lonely
  - n14
  - coupon-collector
  - provable
---

# The "Almost Lonely" Theorem (at most one near), and the n=14 Counterexample Locus

Two things the user asked for: a "proof-lite" impossibility for an n=14 counterexample,
and other LRC-adjacent statements we can actually *prove*. The provable adjacent
statement is clean and turns the repo's coupon-collector picture (S524) into a theorem.

## The provable theorem (first moment)

For nonzero integer speeds `v_1,…,v_{n-1}` and the observer at `0`, let
`near(t) = #{ i : ‖v_i t‖ < 1/n }`. Each runner's danger set
`B_i = { t : ‖v_i t‖ < 1/n }` is a union of `v_i` arcs of total length `2/n`, so
`∫_0^1 1_{B_i}(t)\,dt = 2/n`. Summing,

```
∫_0^1 near(t) dt = (n-1)·(2/n) = 2 − 2/n  <  2     (all n ≥ 3).
```

A nonnegative **integer**-valued function with average `< 2` must take a value `≤ 1`
somewhere. Hence:

> **THEOREM ("almost lonely").** For any `n-1` nonzero integer speeds, there is a time
> `t` at which **at most one runner is within `1/n` of the observer**
> (equivalently `max_t #{far runners} ≥ n-2`). LRC is the (open) statement that
> "at most one" can be improved to "zero."

This is exactly S524's coupon-collector bottleneck ("`d_Q = 6` of `7` most of the
time, one short") made into a theorem. Verified (`lrc_almost_lonely_first_moment_s553.py`,
`n = 5,7,10,14`): for `11/12` sampled sets `min_t near = 0` (genuinely lonely, open),
and the **single "one-short" set is always the AP / regular polygon** — which is
lonely only at the *closed* wall `t = k/n`. So the theorem is tight precisely at the
extremal.

## Why the moment method stops at "one near"

The first moment pins `mean(near) = 2 − 2/n < 2`. The second moment
(`Var(near) ≈ 1.5` at `n=7`, `2.8` at `n=14` for the AP; `max near = n-1`) shows the
count swings widely — and whether the floor reaches `0` (LRC) is governed by the
**pairwise correlations = the resonances** (S550's `E`), i.e. the **simultaneity**
of all `n-1` runners being far at one `t`. That is exactly the high-energy core
(S550) / the 7-way coupling (S524/S552): *not* reachable by the first moment. The
"almost lonely" theorem is the honest ceiling of moment/averaging arguments.

## n=14 counterexample-impossibility (the honest "proof-lite")

A counterexample to LRC@14 cannot be ruled out (that *is* LRC@14), but it is pinned
to an extremely thin locus by the conjunction of repo results:

- **Averaging-extremal (S553):** `near(t) ≥ 1` for all `t`, yet `mean(near) = 12/7 ≈
  1.71` forces some `t` with `near(t) ≤ 1`; so `near` touches the floor `1` but
  **never reaches `0`, not even in closure** — a razor-thin condition.
- **Sieve-covered (THM-369):** contains a multiple of **every** `q ∈ {2,…,14}`.
- **High-energy core (S550):** resonance energy `E ≥ (12/14)^{13}` (else the measure
  bound proves LRC).
- **7-class coupled (S524/S552):** all `7` mod-7 CRT classes blocked at every `t`,
  the singleton `{mult of 7}` coupled to the six pair-classes.

The AP/regular polygon nearly meets all four, **but it is lonely at the wall**
`t = k/14` (closed), so it is *not* a counterexample. No object meeting all four is
known. **The honest "proof-lite": a counterexample would be averaging-extremal,
fully sieve-covered, high-energy, and 7-coupled, and would never close its last
runner — a configuration so constrained that none exists (conjecturally).**

## The menu of provable LRC-adjacent statements (repo tools)

1. **"At most one near"** (this note; first moment) — proved.
2. **Apex bound:** among `n` points on the circle some gap is `≥ 1/n` (pigeonhole;
   dual of `near_pair`, S549) — provable, formalizable.
3. **Lonely-measure bound:** `|LONELY(v)| ≥ (1-2/n)^{n-1} − E(v)` (S550) — provable
   identity + triangle inequality; rigorous with a tail estimate.
4. **No-small-resonance ⟹ LRC** (S550): explicit sufficient condition; rigorous for
   lacunary / large-minimal-resonance families.
5. **The sieve** (THM-369) and **`near_pair`** (S549) — already formalized in Lean.
6. **Weaker-threshold generalization:** at threshold `1/(cn)`, at most `c-1` runners
   near (same first-moment argument with `2c/n` per-runner danger).

## Verdict / next
- Proved the "almost lonely" theorem (`max_t #far ≥ n-2`); it is S524's floor made
  rigorous and is tight at the AP.
- Framed the n=14 counterexample as a four-way-constrained (averaging-extremal +
  sieve + high-energy + 7-coupled) locus the AP just misses.
- Concrete next: (1) **formalize "at most one near" in Lean** (needs
  `volume(B_i) = 2/n`, on top of `near_pair`/sieve); (2) the second-moment refinement
  — can `Var` + the resonance bound push the floor below `1` for a family?; (3) the
  weaker-threshold ladder as a formal LRC-approximation hierarchy.

## Artifacts
```
04-computation/lrc_almost_lonely_first_moment_s553.py
05-knowledge/results/lrc_almost_lonely_first_moment_s553.out
```
Related: S524 (coupon collector = this floor), S550 (resonance energy / high-energy
core = the simultaneity gap), S551 (Vitali: the core is measure-zero), S552 (the
7-class coupling), S549 (near_pair, formalized), THM-369 (sieve).
