---
source: mac-mini-2026-07-06-S37
status: definitive honest-state assessment — the covering route is a reframing (confirmed by klein-S150 + kps-S51); the analytic core (tight scale-uniform decorrelation) is open
tags:
  - lonely-runner
  - LRC14
  - state-of-proof
  - covering-system
  - decorrelation
  - reroute
  - correctness
---

# The honest state of the LRC(14) proof: the covering is a reframing, not a reduction

The owner asked to *deeply consider the state and reroute as needed, ensure
everything is correct, then finish the formalization*. After a rigorous audit,
the honest state is that the fleet's finite-covering route **reformulates** the
crux without reducing it, and the genuine open core is a tight, scale-uniform
**decorrelation** bound. This is now confirmed by three agents (mac-mini-S36,
klein-S150, kps-S51).

## What is actually GREEN (the wired skeleton, opus-S129)

`LRC14Target ⟸ JKReduction(cite) ⟸ (A) ⟸ (C) ⟸ CoveringComplete`, all arrows
GREEN wiring or a named open Prop. The reach atom (`reach_ge_of_covering`, S34) and
the per-`q` band floor (`loose_of_band`, kps) are GREEN: **clearing at any `q` with
`μ/q ≥ 2/25` ⟹ `M ≥ 2/25`.** The reductions `(A) ⟸ (C)` (torus) and the top-level
J-K citation are wired. **The skeleton is real.**

## The finding: CoveringComplete = (C), and it cannot be finitely discharged

`CoveringComplete := ∀ v : Fin 12 → ℤ, ¬DilatedAP v → HasCoveringWitness v`, where
`HasCoveringWitness := ∃ q c μ, 0<q ∧ … ∧ 2q≤25μ ∧ clears` — **no bound on `q`**.

- **Forward** (`CoveringComplete → CruxStatement`): GREEN (opus).
- **Reverse** (`CruxStatement → CoveringComplete`): true — `margin` attains its max
  on `[0,1]` at a rational `c/q` (piecewise-linear), and `margin(c/q) ≥ 2/25` is
  exactly a witness. So **`CoveringComplete ⟺ (C)`** — a *restatement*, not a
  reduction.

- **No finite `q ≤ Q₀` discharge exists.** For every `Q₀`, `V = {i + L·kᵢ}` with
  `L = lcm(2..Q₀)` and all `kᵢ ≥ 1` varying is: non-DilatedAP; **compressed**
  (`max/min → 2`, so it does *not* peel by kps-S49's own `max>13·min` criterion);
  `≡ AP mod L` so it **fails every `q ≤ Q₀`**; and clears **only at
  `q = nextprime(Q₀)`** (loose, `M ≥ ⌈2q/25⌉/q → 2/25⁺`). Verified at `Q₀ = 39`
  (fleet's stated bound): fails all 38 moduli, clears at 41. The covering modulus is
  **unbounded**; kps's `Q₀=25`, klein's `{2..32}`, and the 140k/650k samples are
  height-range artifacts. klein-S144's *"a compressed family cannot be `≡ AP mod`
  large `L` without a far entry"* is false — lifting **all** entries (not some)
  keeps it compressed; *"compressed"* bounds the lift **range**, not its **values**.

**Confirmed independently:** klein-S150 (`{1+L,…,11+L,12+2L}`, `L=lcm(2..37)`,
clears at 41) and kps-S51 (retracted S49/S50: "CoveringComplete == (G)").

## The reroute (klein-S150) and its honest limit

klein: the escape families are **rank-2** (generators `1` and `L`) ⟹ opus's
**(A)-branch** (decorrelation / J-K relative spectrum), not the 1-D `(C)`. Two
caveats make this a *relocation*, not a *resolution*:

1. **They are still in `(C)`'s domain as written.** `CruxStatement` and
   `CoveringComplete` quantify over *all* `Fin 12 → ℤ`. The escape families are such
   families. Routing them to (A) requires the J-K reduction to *provably* only need
   `(C)` for the non-escape subclass — a J-K structural obligation, not yet
   discharged. (LRC14 speeds are unbounded, so unbounded 12-families can arise.)
2. **(A) for them is the tight open core.** The escape families approach `2/25⁺`, so
   the decorrelation that (A) needs must be **tight in the margin** and
   **scale-uniform** (works for `L → ∞`). opus-S127 flagged decorrelation as
   "load-bearing"; it is *not yet proved* for the all-positive varying-`k` case (no
   scale gap: all speeds `~L`). So (A) inherits the hard part.

**Net:** the reroute correctly names the tool (decorrelation) but does not close
anything — it moves the irreducible core from `(C)` to `(A)`.

## Where the difficulty really is (why the covering was always going to fail)

`M ≥ 2/25` for a non-AP family means `safe(V, 2/25) = Leb{t : ‖vᵢt‖ ≥ 2/25 ∀ i} >
0`. The escape families have `safe(V, 2/25) → 0` as `L → ∞` (they approach the
boundary). So there is **no uniform positive measure**, hence no uniform finite
covering — the covering modulus must grow. The genuine statement is a **positivity
of a Poisson/theta sum over the relation lattice** (Cohn-Elkies shape, my S24) that
is *tight* against the escape families. This is analysis, not a residue check.

## Recommendation

- **Relabel** (done in the proof-map): `CoveringComplete` is `(C)` restated, an
  honest open Prop; the "finite `q ≤ 39` residue check" framing is retired.
- **The formalization is a correct conditional skeleton** `LRC14 ⟸ (C)`, with `(C)`
  the open analytic crux. That is as far as the formalization honestly goes until
  the decorrelation / theta-positivity is proved.
- **The real route** is analytic: (a) the tight scale-uniform decorrelation for the
  rank-2 escape families ((A)-branch), and/or (b) the theta-sum positivity for the
  1-D families ((C)), and/or (c) Fan-Sun's *actual* n=4 argument if it genuinely
  extends (the naive covering abstraction of it does not). None is a finite pile.

This is a correction to the fleet's "the crux is now a finite covering system"
(opus-S126) — made in the spirit of getting it right before formalizing. The
skeleton stands; the crux is the analytic core it always was.

→ HYP-4687; confirms/extends HYP-4667/S36; corroborated by klein-S150 (HYP-4651),
kps-S51 (HYP-4677); corrects klein-S144, kps-S49/S50, opus-S126 framing; the reach
atom (S34) + band floor (kps) + skeleton (opus-S129) remain GREEN and correct.
