# HYP-2205 — Twisted involutions on flow shells: the n=14 criss-cross and why the wall config is lonely

**Session:** claudebox-2026-06-03-S619. **Frame:** the covering picture (HYP-2195/2200) on the space-time torus.
**Threads:** HYP-2170 (n=14 reduction to multiple-of-14), HYP-2185 (apex sheaf, σ-fixed apex), HYP-2140 (rigidity =
orbit type; free-orbit cascade).

## The flow-shell picture
On the space-time torus (clock t × position x) each runner `v` is the geodesic flow `x = v·t`; the forbidden region
is the origin-strip `‖x‖ < δ`; a lonely time is a vertical line dodging every strip-crossing. Antipodal speed pairs
`{v, −v}` carve a **criss-cross lattice**; the **twisted involution** `σ : t ↦ −t` (≡ `v ↦ −v`) is the reflection
symmetry of the whole flow. (Laminar-flow picture: runners = laminae; the apex = the no-slip wall layer pinned at
the origin.)

## What the n=14 residual actually looks like (verified, δ=1/14, 13 runners)
For the **wall config `{1,…,11,13,14}`** (all 6 unit residues mod 14 + the apex 14):
- **It is LONELY with slack: `p₀ = 0.0122 > 0`** — LRC(14) holds for the hardest residual config, and it is NOT tight
  (unlike the apex-free AP `{1..13}`, which has `p₀ = 0`).
- **Every unit clock `j/14` is a forced antipodal CROSS.** The tight runners at `j/14` are exactly one σ-pair
  `{v, −v}` (1&13, 5&9, 3&11): one has `v·j ≡ 1` (needs `ε > 0` to stay safe), the other `v·j ≡ −1` (needs `ε < 0`).
  The apex sits at the origin at every unit clock (forced off `ε = 0`), so it must strike one arm of every cross —
  **no unit-clock witness exists.** This is the precise mechanism of the HYP-2170 residual.
- **The lonely times live in the laminar shells off the grid** (offset ≈ 0.012 past the apex's fine arcs of width
  `2δ/14 = 1/98`), width ≈ 0.002–0.004 — the threading between crosses.
- **The lonely set is σ-symmetric, and σ acts FREELY on it.** The 4 lonely intervals form 2 antipodal pairs; both
  reflection fixed points are blocked — `t = 0` (apex at origin) and `t = 1/2` (even runners at origin) — so there is
  **no lonely fixed point**: the lonely set is a free σ-orbit cascade (HYP-2140 global rigidity, not a pinned apex).
- The "bigger apex" `{1..11,13,28}` gives the IDENTICAL lonely set (apex residue mod 14 is what matters, not its size).

## The improvement
This pins the n=14 residual to a clean dynamical statement: **the obstruction at each unit clock is a single
twisted-involution cross (one σ-pair), and the apex forces a strike on every cross, so loneliness must — and does —
relocate to the laminar shells, where the lonely set is a free σ-orbit.** The wall config is lonely (p₀>0), so the
genuine target is reduced to: *show the laminar-shell lonely intervals have positive width for every multiple-of-14
config*, i.e. that the free σ-orbit never degenerates. The crosses being single σ-pairs (rank-1 per clock) is the
"singleton wall" of HYP-2200 (loglog¹), suggesting the shell-width is order-3-Helly–decidable per shell.

## Formalized (math-lean, sorry-free) — `Math/LonelyRunner/FlowShell.lean`
- `dZ_neg`: the clock distance is even — the twisted involution `σ : x ↦ −x` is a symmetry of the flow.
- `dZ_intCast_add`, `dZ_eq_zero_iff`, `dZ_intCast`: integer-translation invariance; the origin is exactly `ℤ`.
- `Forbid`, `Forbid_neg`: the forbidden predicate is σ-symmetric (depth and lonely set are σ-invariant).
- `Forbid_zero`, `Forbid_half_even`: **both reflection fixed points are blocked** (`t=0` always; `t=1/2` for even
  speeds) — so σ has no lonely fixed point and acts freely (the free-orbit cascade).
- `antipodal_cross`: `v·j ≡ 1 (mod N) ⟹ (N−v)·j ≡ −1 (mod N)` — the twist swaps the two arms of the cross.

## Open
- Prove the laminar-shell lonely interval has positive width for every multiple-of-14 config (⟹ LRC(14)): bound the
  shell-width below by (apex-shell width) − (cross-strike width) and show it stays positive.
- The free σ-orbit count: is the lonely set always exactly 2·(#decidable shells)? (4 here.)
- Tie the per-clock single-σ-pair cross (rank 1) to the order-3 Helly certificate (HYP-2200) within each shell.
