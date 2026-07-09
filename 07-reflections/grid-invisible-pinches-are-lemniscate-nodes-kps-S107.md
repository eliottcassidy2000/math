# Grid-invisible pinches are lemniscate nodes; the smooth surrogate is their elliptic desingularization

*kind-pasteur-2026-07-09-S107. Owner: push the LRC math, then formalize; connect "grid-invisible
pinches" to the lemniscate of Bernoulli; be abstract. mac-mini-S64 (MISTAKE-130) had just found that
`maxgap(x)=1/7` is hit at "measure-zero rational pinches invisible to every uniform grid" — which is
why the widest-arc pigeonhole over-merged arcs. This note names the geometry of those pinches, and it
resolves the obstruction: don't sample the singular `maxgap`; sample its **elliptic desingularization**.*

---

## The three verified facts (`lrc14_pinch_{lemniscate,k13}_kps_S107`)

1. **The good-set boundary is a set of rational pinches.** The arcs of `G* = {x : maxgap{frac(e_i x)} >
   1/7}` end exactly where `maxgap(x) = 1/7`, and these abscissae are **rationals `x = m/d`** with `d` a
   cluster difference `e_i − e_j` (verified: `5/21, 16/21, 5/28, 23/28, 1/35, 6/35, …` for a spread-35
   `k=13` cluster). Their denominators are **cluster differences, not the ruler `Vmax`** — so they are
   *invisible to every `Vmax`-ruler grid* `x_j = (j+½)/Vmax` (mac-mini's phrase, made exact).
2. **Each pinch is a node — two gap-functions crossing.** At a tooth collision `x* = m/d` two teeth
   coincide and their adjacent gaps **swap order** (verified: teeth order reverses through `x*`). So
   `maxgap` has a **corner** there (measured slopes `−12 → +9` at `x*=5/21`): two linear branches
   crossing. That is precisely the **lemniscate node**: `(x²+y²)² = x²−y²` to leading order at the
   origin is `x²−y² = 0`, i.e. the two crossing lines `y = ±x`. **The grid-invisible pinch and the
   lemniscate's self-crossing are the same local singularity.**
3. **The pinches are countable — `#arcs/spread` is dilation-invariant.** Scaling `E → sE` gives
   `#arcs = 32,64,128,256` at `s = 1,2,4,8` — `#arcs/spread ≡ 0.914`, exactly linear. This is the
   Davenport–Schinzel upper-envelope bound (`#arcs = O(spread)`) realized as a **clean dilation
   invariant** of the cluster (same axis as `D3`, the longest-AP axis).

## Why grids fail, and how the lemniscate says to fix it (MISTAKE-130 resolved)

A uniform grid samples `maxgap` at `x_j`. The good set's *boundary* lives at the rational pinches
`m/d`, a measure-zero set the grid generically misses — but the grid can also *land on* a pinch (when
`d ∣ 2Vmax`) and read the corner value ambiguously. So the grid neither counts the arcs nor measures
`meas(G*)` faithfully (MISTAKE-130: the 117k widest-arc search over-merged arcs). This is the exact
pathology of sampling a function **on the real axis through a node**: the real slice sees a corner /
crossing, not the smooth object.

The lemniscate is the paradigm of the fix. Its node at the origin is invisible/singular to the naïve
Cartesian slice, but the curve is **uniformized smoothly** by the lemniscatic sine `sl(u)` (an entire
function, CM by `ℤ[i]`), and its arc length `∫₀^r dt/√(1−t⁴)` — the elliptic integral with total
`ϖ = 2∫₀¹dt/√(1−t⁴) ≈ 2.622` — is exact *because* it runs on the desingularized (elliptic) coordinate,
not the singular Cartesian one. **The moral for LRC:** do not integrate the singular `maxgap`; integrate
its **desingularization** — the continuous, bounded-variation surrogate.

That surrogate already exists: **opus-S170's smooth `maxgap`/`W`** (`W = Σ(gap−1/7)_+`), which is `C⁰`
through every pinch (verified: `W` has a finite-slope corner where `maxgap` has its node) with **Fourier
decay `1/m²`** (α=2, kps-S97/opus-S170). `W` is to `maxgap` what `sl(u)` is to the nodal lemniscate:
the entire/smooth resolution whose *grid-average converges to the true measure*. So the equidistribution
`ρ_K → ρ*` — the sole remaining Part-A node (hembed's tight window) — should be run on `W`, not on the
sharp indicator: the smooth surrogate's `1/m²` resonant sum converges absolutely (Mertens-safe,
kps-S93), and the grid-average discrepancy is the desingularized (finite, `#arcs = O(spread)`) count of
pinches, not the invisible singular set.

## The abstract picture

| lemniscate `(x²+y²)²=x²−y²` | LRC good set `G*` |
|---|---|
| origin node `y = ±x` | tooth-collision pinch (two gaps cross), rational `x = m/d` |
| singular Cartesian slice | uniform ruler grid `x_j = (j+½)/Vmax` (blind to the node) |
| `sl(u)` entire uniformization | opus-S170 smooth surrogate `W` (`C⁰`, Fourier `1/m²`) |
| arc length `∫dt/√(1−t⁴)`, constant `ϖ` | `ρ* = meas(G*)` — the true good measure |
| `ℤ[i]` CM / doubling `r²=cos2θ` | the `×2/×7` dilation-collapse maps; `#arcs/spread` dilation-invariant |
| node desingularizes to `ℙ¹` | `#arcs = O(spread)` (Davenport–Schinzel), pinches countable |

The through-line: **a grid-invisible pinch is a node; measuring `G*` by sampling the singular `maxgap`
on a grid is slicing a node on the real axis (fails); the fix is the elliptic desingularization — the
smooth surrogate `W`, whose `1/m²` Fourier decay makes its grid-average equal the true measure.** This
converts mac-mini's MISTAKE-130 from an obstruction into a prescription, and it tells the whole fleet
that the last Part-A node (equidistribution, hembed's tight window) is a *smooth-surrogate* computation,
not a sharp-indicator one.

## Actionable

- The equidistribution `ρ_K → ρ*` / discrepancy bound should be stated and proved for the **smooth `W`**
  (opus-S170), whose Fourier `1/m²` is proven-in-spirit (continuous PL), giving an absolutely
  convergent resonant sum and a `#arcs = O(spread)`-controlled grid error. This is hembed's tight-window
  residual (kps-S106 handled the large-ruler regime).
- Formalizable next step: `#arcs/spread` is a **dilation invariant** (proved computationally, dilation
  `E → sE` gives `#arcs → s·#arcs`); and the pinch abscissae are `{m/d : d ∈ diffs(E)}` — a finite,
  a-priori set. Both are clean Lean targets that make the pinch count concrete.

*Files: `lrc14_pinch_lemniscate_kps_S107.py`, `lrc14_pinch_k13_kps_S107.py` (+ `.out`). Builds on
mac-mini-S64 (MISTAKE-130, the pinches), opus-S170 (smooth surrogate, α=2), kps-S93/S97 (Mertens-safe
resonant sum, kissing), kps-S95 (the lemniscate as the second-moment/doubling curve). See
[[triangle_foundation]].*
