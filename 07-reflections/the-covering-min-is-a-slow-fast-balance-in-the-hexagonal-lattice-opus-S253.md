---
source: opus-2026-07-11-S253
status: GEOMETRIC SHAPE + creative argument (a PROOF for a structured class, a direction for the rest). The
  covering-min target M >= n/Phi6(n) has a clean geometric shape: (i) M = the L-infinity CLEARANCE of the
  closed geodesic gamma(t)=(v_i t) through the integer-hyperplane arrangement on T^13; (ii) the value
  n/Phi6(n) is a HEXAGONAL-lattice norm, Phi6(n) = N(n - zeta_6) in Z[zeta_6]; (iii) it is a SLOW-FAST BALANCE
  M = M_core * v_f/(v_f+s), which PROVES M >= n/Phi6 for the interval-core single-killer class (deep well =
  unique minimizer, via smallest-killer + monotonicity + LRC(n-1)), and reduces the general lower bound to a
  joint control M_core*v_f/(v_f+s) >= n/Phi6 -- an inductive target bootstrapping from LRC(n-1).
tags:
  - lrc14
  - covering-min
  - geometry
  - hexagonal-lattice
  - eisenstein
  - slow-fast-balance
  - loop-clearance
  - inductive
---

# The covering-min is a slow-fast balance in the hexagonal lattice

**opus-2026-07-11-S253.** Owner: understand the shape of the target and find creative geometric/topological
arguments toward proofs. The target (S252) is the covering-min bound `M ≥ n/Φ₆(n)` (`= 14/183` at `n=14`).
Three geometric pictures, the third of which gives an actual proof for the structured extremizers.

## The shape: M is the L∞-clearance of the loop

`M(v) = max_t min_i ‖vᵢ t‖` is the **L∞ clearance** of the closed geodesic `γ(t) = (v₁t,…,v₁₃t)` on the torus
`T¹³` through the integer-hyperplane arrangement `⋃ᵢ{xᵢ ∈ ℤ}`: the size of the largest **safe box** the loop
threads. LRC(14) = the loop reaches a box of size `1/14`; the covering-min = the covering family whose loop
threads the **smallest** box. Each coordinate `xᵢ = vᵢt` crosses a hyperplane `vᵢ` times per period, so the
loop is chopped into `Σvᵢ` arcs, and `M` is the tallest clearance-peak over them.

## Picture 1: the hexagonal / Eisenstein lattice

The covering-min value is `n/Φ₆(n)` with `Φ₆(n) = n²−n+1 = N(n − ζ₆)`, the **norm of `n − ζ₆` in the hexagonal
lattice `ℤ[ζ₆]`** (Eisenstein integers). Verified `n=7..14`: `Φ₆(n)` is always `a²−ab+b²` with rep `(1,n)`. So
the deep well `{1..n−2, n(n−1)}` **is** the Eisenstein integer `n − ζ₆`, and the covering-min is `n` over its
hexagonal norm. The optimal phase set is the comb `{n·k/Φ₆} ∪ {killer}` with three-gap `g=2`, gaps
`{1/Φ₆, n/Φ₆}` — the **hexagonal fundamental domain**. (This is why `Φ₆` — the cyclotomic/Eisenstein norm form
— governs the covering-min; ties to kps's Eisenstein/three-distance and the `X₀(14)` genus.)

## Picture 2 → the proof: the slow-fast balance

Write a covering family as a **slow core + a resonant killer**. Take the interval core `{1..n−2}` (its LRC
optimum is `t₀ = 1/(n−1)`, value `1/(n−1)`) plus a killer `v_f` that is **resonant** at `t₀` (`(n−1)∣v_f`, so
`v_f t₀ ∈ ℤ`). Perturb `t = t₀ + δ`:

- **killer clearance rises:** `‖v_f(t₀+δ)‖ = v_f|δ|`;
- **core binding runner `v=1` falls:** `‖1·(t₀+δ)‖ = 1/(n−1) − |δ|`.

They **cross** at `|δ| = 1/((n−1)(v_f+1))`, giving

> **`M = 1/(n−1) − |δ| = v_f/((n−1)(v_f+1)).`**  (verified exactly for `v_f = 182,364,546,1820,2730`.)

**Lower bound (proved for this class).** Covering needs a multiple of `n−1` *and* of `n`; the interval core has
neither, so the killer must be a multiple of `lcm(n−1,n) = n(n−1)`, i.e. `v_f ≥ n(n−1)`. Since
`v_f/((n−1)(v_f+1))` is **increasing** in `v_f`, the minimum is at `v_f = n(n−1)`:

> **`M ≥ n(n−1)/((n−1)(n(n−1)+1)) = n/(n(n−1)+1) = n/Φ₆(n),`** equality iff `v_f = n(n−1)` — the **deep well**.

So every interval-core single-killer covering family has `M ≥ n/Φ₆`, and the deep well is the **unique
minimizer**. (Non-resonant killers give `M ≥ 1/(n−1) > n/Φ₆` — the killer is already safe at `t₀`, no balance
needed; verified `v_f=183,185,200 → M=1/13`.) This **derives `Φ₆`** (`= v_f + 1 = n(n−1)+1`), makes mac-mini
S40's "two-point equioscillation" explicit as a 1-D balance, and **bootstraps from LRC(n−1)** — the core value
`1/(n−1)` is the LRC(n−1) extremal (proved for `n−1 ≤ 13`).

## Picture 3 → the direction: the general balance is inductive

For a general core with LRC optimum value `M_core` and binding speed `s`, plus a killer `v_f` resonant at the
core optimum, the same balance gives

> **`M = M_core · v_f/(v_f + s).`**

The deep well is extremal on **all three knobs**: `M_core` minimal (interval core = LRC(n−1) extremal,
`1/(n−1)`), `s = 1` minimal, `v_f = n(n−1)` minimal. Any deviation raises a factor. So the **general covering-
min lower bound reduces to**: control `M_core · v_f/(v_f+s) ≥ n/Φ₆` jointly — an **inductive** target that
bootstraps from **LRC(n−1)** (bounding `M_core`) plus the killer-clearing balance. The one open escape it must
close is a **large-`s` trade** (a core whose binding runner is fast, shrinking `v_f/(v_f+s)`) against a larger
`v_f` — the concrete geometric obstruction, now isolated. (For multi-killer families the balance becomes a
multi-constraint version: clear several resonant runners simultaneously — the 2D `line-hits-boxes` picture
generalized, the natural next computation.)

## Net

The covering-min has a clean geometric shape — the L∞ clearance of the loop, valued at a **hexagonal-lattice
norm**, realized as a **slow-fast balance** `M = M_core·v_f/(v_f+s)`. The balance **proves** `M ≥ n/Φ₆` for the
interval-core single-killer class (deep well = unique minimizer, from smallest-killer + monotonicity +
LRC(n−1)) and reduces the general bound to an **inductive** joint control bootstrapping from LRC(n−1) — a
genuine geometric route toward the covering-min, with the sole remaining obstacle (the large-`s` trade)
pinned. This is more than a reframing: it is a proof for the extremal structured class and an inductive
strategy (LRC(14) covering-min ← LRC(13) + balance) for the rest.

→ mac-mini S38 (Ostrowski ladder), mac-mini S40 (2-point equioscillation — here made a 1-D balance), klein
S267 (14/183 covering-min, verified), kps (Eisenstein / three-distance / `X₀(14)`), THM-366 (covering ⟹
divisor-complete), THM-527 (three-gap), opus-S252 (target relocation), LRC(≤13) citation (the core bound).
Files: `lrc14_covering_min_slow_fast_balance_opus_S253.py` (+`.out`).
