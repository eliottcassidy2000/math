---
source: oracle-2026-06-02-S557o
status: progress (case-wise gap lower bounds; the sieve gives the FULL conjecture off a thin core; honest positioning vs Tao's general bound)
tags:
  - lonely-runner
  - gap-bound
  - tao
  - sieve
  - first-moment
  - case-analysis
---

# Improving the Lonely-Runner Gap Bound *in Cases*: the Sieve Gives the Full Conjecture off a Thin Core

User: Tao improved the **general** gap lower bound from the trivial
`g(v) := sup_t min_i ‖v_i t‖ ≥ 1/(2k)` (k = # nonzero speeds) to
`1/(2k) + c·log k / (k² (log log k)²)` — a tiny gain. *Can we improve it further, and
what about in certain cases?*

## The honest frame: the trivial bound is the first moment

At threshold `θ` each runner's danger arc `{t : ‖v_i t‖ < θ}` has measure `2θ`, so the
**expected near-count** is `E[#near] = 2kθ`. `E[#near] < 1 ⟺ θ < 1/(2k)`, and a
nonnegative integer with mean `< 1` is `0` somewhere — *that is the trivial bound*
`g ≥ 1/(2k)`. Tao pushes `θ` a hair above `1/(2k)` by lower-bounding the **overlaps**
of the danger arcs (a second-order term). Beating his constant in general is a delicate
optimization I do not claim to win. **But "in certain cases" the repo's sieve gives the
full conjecture `g ≥ 1/(k+1)`, and structure gives far more.**

## The case hierarchy (rigorous lower bounds on `g`)

| case | bound | note |
|------|-------|------|
| (C0) first moment | `g ≥ 1/(2k)` | always; = the trivial bound |
| (C0′) Tao | `g ≥ 1/(2k) + c log k/(k²(loglog k)²)` | always; tiny gain |
| **(C1) sieve at `q*`** | **`g ≥ 1/q*`**, `q* =` least `q∈{2..k+1}` dividing **no** speed | `∞` iff sieve-covered |
| (C2) all speeds odd | `g ≥ 1/2` | `q*=2` |
| (C3) no speed `≡0 (mod k+1)` | `g ≥ 1/(k+1)` exactly | `q*=k+1` |
| (C4) no mult of 2 or 3 | `g ≥ 1/2` | smallest missed `q=2` |

**(C1) is the lever.** At `t = 1/q*`, every speed sits at `‖v_i/q*‖ ≥ 1/q*` (since
`q* ∤ v_i`), so `g ≥ 1/q*`; and if the set is *not sieve-covered* then `q* ≤ k+1`, hence

> **`g ≥ 1/(k+1)` — the full conjecture — for every non-sieve-covered speed set.**

This *beats Tao outright* (the full `1/(k+1)`, not `1/(2k)+ε`) for the entire bulk.

## What the computation shows (`lrc_gap_lower_bounds_by_case_s557.py`, n=14)

- **38/40** random primitive 13-sets are **not sieve-covered** → `g ≥ 1/14` proven by
  the sieve. Actual gaps `0.15–0.29`, i.e. `~2–4×` Tao's `1/26 ≈ 0.0385` and all `≥ 1/14`.
- **Structured:** all-odd → `g = 1/2`; no-mult-of-14 → `g = 1/14` exactly; no-mult-of-2-or-3
  → `g = 1/2`. Rigorous and far above the conjecture.
- **The wall (AP `1..13`):** sieve-covered, `g = 1/14` exactly — it *saturates* the
  conjecture, so **no general bound can ever exceed `1/(k+1)`.** Tao's bound necessarily
  lives strictly below the wall.
- **Key nuance:** the two sieve-covered random sets still have `g = 0.16, 0.23` — far
  above `1/14`. **Sieve-covered does NOT mean near the wall.** The genuinely hard locus
  (gap → `1/(k+1)`) is the even thinner *near-AP* slice inside the sieve-covered sets.

## The verdict: "improve in certain cases" is proven, and Tao's regime is tiny

- For the **bulk** (non-sieve-covered), the sieve proves the **full conjecture**
  `g ≥ 1/(k+1)` — a large improvement over Tao for those sets.
- **Sieve-covered** sets (a multiple of every `q ≤ k+1`) are the only locus where Tao's
  general bound is the operative one — and even most of *those* have `g ≫ 1/(k+1)`. The
  truly extremal locus is the near-AP slice. So Tao's hard-won general gain is needed
  only on a very thin, near-extremal set of configurations.
- On that thin core, the repo's complementary tools push toward `1/(k+1)`:
  the **measure bound (S550)** gives `g ≥ θ` for any `θ < 1/n` with resonance energy
  `E(v) < (1-2θ)^{k}` (a `θ`-deformation of the loneliness measure), and the **local LP
  (S556)** clears the `≤1` near runner off the AP. Neither closes the AP itself, where
  `g = 1/n` is tight.

## Why this reframes the problem

Tao's improvement is usually read as "the general bound is `1/(2k)+ε`." The sieve says
the *general* bound is misleading: **almost every speed set already satisfies the full
conjecture**, by a one-line argument (`t = 1/q*`). The `1/(2k)+ε` regime is the
sieve-covered, near-AP corner — a measure-tiny, arithmetically rigid family. The right
target is not "raise `1/(2k)` for all `v`" but "close the near-AP core," which is exactly
where S550/S553/S556 already point.

## Verdict / next
- **In cases (the bulk): `g ≥ 1/(k+1)` proven by the sieve; structured cases reach `1/2`.**
  This is a genuine, large improvement over Tao *for those sets*.
- General constant: not beaten (Tao's overlap optimization stands).
- Concrete next: (1) the `θ`-deformed measure bound — find the largest `θ < 1/n` with
  `E(v) < (1-2θ)^k` over the sieve-covered core, giving an explicit `g ≥ θ` there;
  (2) quantify "near-AP": which sieve-covered sets actually have small gap, and bound
  their distance from the AP below; (3) combine the sieve (`q ≤ n`) with one fine pinch
  (`q > n`, S555/S556) to cover the sieve-covered non-AP sets.

## Artifacts
```
04-computation/lrc_gap_lower_bounds_by_case_s557.py
05-knowledge/results/lrc_gap_lower_bounds_by_case_s557.out
```
Related: THM-369 (denominator sieve = the `t=1/q*` bound), S554 (A1 sieve-covered = the
only locus where Tao's bound is operative), S550 (resonance-energy measure bound, the
`θ`-deformation), S553 (at most one near), S556 (local LP / fine pinch), S551 (the wall).
