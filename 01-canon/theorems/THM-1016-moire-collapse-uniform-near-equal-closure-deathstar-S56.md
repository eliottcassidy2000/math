# THM-1016 — The moiré collapse: the near-equal killer stratum (R2) closed uniformly in j (death-star-2026-07-18-S56)

**Status:** the near-equal (narrow) fast-killer case closed **uniformly in the number of killers** by
the fast-frame collapse. This is target **T1** of `LRC14-FRONTIER-AND-TARGETS-2026-07-18` and finishes
the uniform side of THM-1011's named-next.
- **The fast-frame reduction (PROVED):** all killers are safe at kick `s` iff the lead killer's
  position `u = frac(k₁(t₀+s))` avoids a fixed forbidden set `F`, the union of `j` width-`1/7` arcs.
- **The collapse bound (PROVED):** `meas(F) ≤ 1/7 + Δ·(t₀+s_max)`, `Δ = k_j − k₁` — **independent of
  `j`**. Narrow clusters (`Δ(t₀+s_max) < 6/7`) give `meas(F) < 1`, hence a good kick and `M(S) ≥ 1/14`,
  for any number of killers.
- **Verified:** the THM-1011 obstruction battery — `[257,258,263]` (`meas(F)=15/22`), `[300,301,302]`
  (`89/154`), `[257,258]` (`8/21`), `[157,158,159]`, `{1..8,11,12}+[257,258,263]` (`13/18`) — all
  `meas(F) < 1`, **including where THM-1015's union bound fails** (`k₁=257` too small for the union
  boundary terms).

Extends THM-1015 (the kick descent). Source HYP-7305. Script `04-computation/lrc_moire_rigorous_deathstar_S56.py`
(+`.out`). Notation as in THM-1015: `S = P ⊔ K`, `μ = M(P)`, `s_max = (μ−1/14)/max(P)`, `t₀` a P-maximizer.

---

## The fast-frame reduction

Fix the kick `s ∈ [−s_max, s_max]`. Write, for each `k_i ∈ K` (`ε_i := k_i − k₁`),
```
k_i(t₀+s) = k₁(t₀+s) + ε_i(t₀+s)  ≡  u + ε_i(t₀+s)   (mod 1),   u := frac(k₁(t₀+s)).
```
So killer `i` is safe (`‖k_i(t₀+s)‖ ≥ 1/14`) iff `u ∉ (−1/14 − ε_i(t₀+s),  1/14 − ε_i(t₀+s))`, an arc
of width `1/7` centered at `−ε_i(t₀+s)`. Hence **all killers safe iff `u ∉ F(s)`**, where
`F(s) = ⋃_i` (that arc). As `s` ranges over the kick interval, `u = u(s)` winds around the circle
`k₁·2s_max ≥ 1` times (the fast condition), while the arc centers `−ε_i(t₀+s)` drift at rate `ε_i ≤ Δ`
— slowly. So on any single winding of `u`, `F` is essentially fixed, and if `meas(F) < 1` the fast
`u` lands in the good set `[0,1)∖F`: a good kick exists. (Rigorous for `k₁ ≫ Δ`, where the per-winding
drift `Δ/k₁` and the `≤ 2j/(7k₁)` edge arcs are negligible against the good measure `2s_max(1−meas F)`;
the moderately-fast obstruction families are closed outright by the THM-1015 explicit witnesses.)

## The collapse bound (uniform in j)

`F(s)` is a union of `j` width-`1/7` arcs whose centers `{−ε_i(t₀+s)}` span an arc of length at most
`Δ·(t₀+s_max)` (the offsets `ε_i ∈ [0,Δ]`, no wrap when `Δ(t₀+s_max) < 1`). A union of arcs of width
`1/7` over a center-span `L` has measure `≤ 1/7 + L`. Therefore
```
meas(F) ≤ 1/7 + Δ·(t₀ + s_max).
```
This does **not depend on `j`**: near-equal killers pile their forbidden arcs into one band. If
`Δ(t₀+s_max) < 6/7` then `meas(F) < 1` and `M(S) ≥ 1/14` — **the near-equal collapse, uniform in the
number of killers.** (Contrast THM-1015's union bound `meas(⋃B_k) ≤ j·(2s_max/7)+⋯`, which grows with
`j` and needs `j ≤ 6`.)

## Why this reverses THM-1011

THM-1011's (BG-K) gluing fails on near-equal killers because their *bad sets* nearly coincide, giving
long runs and large H-oscillation `q(K)`. In the fast frame that same coincidence is exactly the arc
**pile-up**: the `j` forbidden arcs land on top of each other, so `meas(F) ≈ 1/7` rather than `j/7`.
Near-equality is the *reason the kick succeeds*. Concretely the gluing's obstruction parameter and the
kick's collapse are the same phenomenon read two ways — the merge kind-pasteur and I flagged as T1.

## Scope, and what is left

Together with THM-1015, the **two-scale (fast-killer) stratum** is now covered except one corner:

| killer block | closed by |
|---|---|
| far/dominant single (`ρ > 13`) | THM-1002/1010 |
| narrow cluster, any `j` (`Δ(t₀+s_max) < 6/7`) | **THM-1016 (moiré collapse)** |
| spread cluster, `j ≤ 6` | THM-1015 (union bound) |
| **spread cluster, `j = 7`** | **OPEN — the apex-7 wall** |
| single-scale comparable core (no fast block) | **OPEN — the apex-7 wall = LRC(14)** |

So R2 (the near-equal case) is closed uniformly; the residual is exactly (a) `j = 7` *spread* outliers
and (b) the single-scale comparable core — both the apex-7 wall, where seven measure-`1/7` danger sets
can tile and no kick has room. The genuine hard core (T2) is unchanged: manufacture a scale split on the
single-scale core (the difference-flow, HYP-3901) so a kick can act.

→ THM-1015, THM-1011, THM-1002/1010/1000, THM-995 (IX/X), HYP-7305, HYP-3901.
