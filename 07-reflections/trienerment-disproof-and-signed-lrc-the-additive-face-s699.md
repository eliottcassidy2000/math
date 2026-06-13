---
source: opus-2026-06-06-S699 (remote-control)
status: ANSWERS the user's disproof-simplification question + sharpens signed-LRC (extending S674/S674b). The trienerment '=1' layer = the NORM-1 layer; 2D lattices are kissing-capped (κ≤6 ⟹ ≤3n = Harborth); the CM disproof escapes by making '=1' a FIELD-NORM condition (unbounded). So the trienerment simplifies the CONCEPT of the disproof (kissing-cap → field-norm), not the arithmetic. Signed LRC: gauge (per-runner signs) preserves M (theorem) and exposes pair-sums = the 2n-1 worry-set shells, but is a side-channel, NOT a witness-reduction (the frame shift changes M, verified).
tags: [trienerment, three-state, unit-distance, disproof, kissing-number, norm-1, CM-field, class-field-tower, signed-lrc, gauge, pair-sum, 2n-1, eisenstein, additive-face]
---

# Trienerment & the disproof; signed LRC & the additive face

**Prompt (user):** (1) LRC with runners able to move backward (negative speed / every-other
clockwise) — tricks live here. (2) unit distance as a trienerment: `=1` = edge, `<1`/`>1` =
opposite arrows; would this simplify the recent disproof of the grid as optimal?

Both reframes were partly built by S674/S674b (signed LRC → pair-sum clocks; the `<1/=1/>1`
threshold tournament). Here I answer the **disproof question** (unaddressed) and pin the exact
boundary of the signed-LRC trick.

## 1. Does the trienerment simplify the disproof? Partly — it localizes the *concept*.

The trienerment's `=1` (tie/wall) layer **is the norm-1 layer** of the point set. For a **2-D
lattice** the per-point `=1` degree is the **kissing number** `κ`, and `κ ≤ 6` (verified
`…s699a.py`):
```
   square ℤ² (Gaussian)        κ = 4  ⟹  unit distances ≈ 2n
   triangular (Eisenstein ζ₆)  κ = 6  ⟹  unit distances ≈ 3n     ← the 2-D lattice maximum
```
> So **every 2-D lattice gives ≤ 3n unit distances** (linear; `Harborth = ⌊3n−√(12n−3)⌋` is exactly
> this cap minus the perimeter). The grid is optimal *among 2-D lattices* because `κ_max = 6`.

**The trienerment reading of why the grid is NOT globally optimal:** the `=1` layer being capped at
`κ≤6` is a **Euclidean kissing** fact. The CM construction (Sawin/OpenAI) is **not a 2-D lattice** —
it is a 2-D image of a **high-rank unit group** in which `=1` means **field-norm 1** (`β·β̄ = 1`).
By Hilbert-90 (`β = γ/γ̄`) plus infinite class-field towers (Golod–Shafarevich), the field-norm-1
layer is **unbounded** — it is *not* kissing-capped. Hence `n^{1.014}` unit distances.

> **The trienerment localizes the disproof to one sentence:** *replace the Euclidean-kissing `=1`
> layer (cap `κ≤6`) by a field-norm `=1` layer (no cap).* Equivalently, the unit-distance graph is
> the **norm-1 Cayley graph of the point group**; the grid (`ℤ²`) is a *suboptimal group* for
> norm-1 packing (`κ=6`), the CM field is a *better group* (norm-1 unbounded via the tower). The
> toy `ℚ(ζ_m)` has `2m → ∞` roots of unity on `|z|=1`, vs the 2-D cap `6`.

**Honest verdict.** The trienerment **simplifies the conceptual statement** of the disproof (a
group-with-a-denser-norm-1-layer, escaping the geometric kissing cap) and explains *why* a 2-D
lattice can't win. It does **not** simplify the *arithmetic* core (the `n^{1.014}` count is the
class-field-tower / Golod–Shafarevich estimate, untouched by the combinatorial reframe). So: yes for
the **why**, no for the **how-many**.

## 2. Signed LRC — the exact boundary of the negative-speed trick

**Gauge (per-runner signs) preserves `M` — a theorem.** For any signs `ε_i∈{±1}`,
`‖ε_i v_i t‖ = ‖v_i t‖` for every `i` and `t`, so `min_i` is identical pointwise and
`M({ε_i v_i}) = M(S)` exactly. (My `M_exact` returned a spurious "False" here — a *bug*: it builds
breakpoints `m/d` with `range(d)`, empty for negative `d`, so it silently drops witnesses on
negative speeds. The identity is one line and certain.) So the gauge is an exact `(ℤ/2)^{n−1}`
symmetry of loneliness — and it preserves which residues are `0 (mod n)`, so it **cannot remove a
multiple of `n`** (no free C2b reduction).

**The gauge exposes the pair-SUMS = the worry-set shells.** The runner-runner clock is
`ε_i v_i − ε_j v_j ∈ {v_i−v_j, v_i+v_j}`; opposite signs give **sums**. **Verified** (AP13,
`C=2n−1=27`): opposite-sign sums fill the interior residues `3..25` mod 27 — the pair-sum sieve
(THM-401), i.e. the **worry-set shell ledger made geometric as opposite-direction collisions**.
This is the real trick (S674/S674b): observer-blind, pair-visible.

**The honest boundary (the tempting false move).** "Sit on another runner / reverse the frame" is a
**frame shift** (subtract `v_k`), which **changes the config and `M`** — verified: `V=(1,2,3,7)`,
`M=1/4`; shifts give `M = 2/7, 1/3, 1/5, 4/11`. So the frame shift is *not* an `M`-preserving
witness transfer; it relates *different* configs. **The gauge is a side-channel (it classifies
configs by their pair-sum address), not a witness-finder.** The negative-speed idea animates the
additive face; it does not, by itself, prove a row loose.

**Complex/cyclotomic view (the helix).** Signed runners are complex frequencies
`e^{2πi ε_i v_i t}`; the witnesses are `(2n−1)`-th roots of unity (THM-403, cyclotomic), and the
gauge chooses which **root-products** `e^{2πi(ε_i v_i+ε_j v_j)t}` resonate. The pair-sum modulus
`2n−1` is the order of the root the opposite-direction pairs beat against.

## 3. The unification: both reframes expose the ADDITIVE face as a 3-valued lift

| | LRC (signed) | unit distance (trienerment) |
|---|---|---|
| binary object | runner danger `‖v_i t‖ ≷ 1/n` | distance `=1` or not |
| 3-valued lift | sign `ε_i` ⟹ clocks `v_i ∓ v_j` | `<1 / =1 / >1` threshold |
| additive face exposed | pair-**sums** mod `2n−1` (shells) | **norm-1** layer (`=1` ties) |
| the cap / the escape | worry-set lives at the sums | kissing `κ≤6` vs field-norm ∞ |
| cyclotomic anchor | `(2n−1)`-th roots (THM-403) | Eisenstein `ζ₆` / CM roots |

Both are the **additive resonance** (pair-sum / norm) face, lifted from a binary to a 3-valued
(signed / threshold) structure, and both anchor on **roots of unity** — the same `π/3`/cyclotomic
object as S599u (`Cl₂(π/3)`, the Eisenstein angle). The trienerment's `=1` ties are the unit-norm
roots; the signed LRC's sum-clocks are the `(2n−1)`-th roots; the worry-set and the unit-distance
optimum are both "maximal root-of-unity packings."

## 4. Honest status

- **Verified:** 2-D lattice kissing numbers (square 4, triangular 6); the `≤3n` cap = Harborth; the
  frame-shift changes `M`; the gauge pair-sum split (`3..25` mod 27).
- **Theorem (one line):** gauge invariance `M({ε_i v_i}) = M(S)` (the `M_exact` "False" was a
  negative-speed bug, noted).
- **Answer to the disproof question:** the trienerment **simplifies the concept** (kissing-cap `→`
  field-norm escape; grid = suboptimal norm-1 group) but **not the arithmetic** (`n^{1.014}` is the
  class-field-tower count, unchanged).
- **Honest limits:** the signed gauge is a side-channel (classifies configs by pair-sum address),
  **not** a C2b witness-reduction; the frame shift is not `M`-preserving. New here vs S674/S674b:
  the disproof = kissing-cap-vs-field-norm reframe, the norm-1-Cayley-group reading, the exact
  gauge/frame boundary, and the cyclotomic unification.

**Artifacts:** `04-computation/trienerment_disproof_kissing_cap_s699a.py`,
`signed_lrc_gauge_vs_frame_s699b.py` (+`.out`s). Builds on S674/S674b (signed LRC / threshold
tournament), THM-401 (`2n−1`), THM-403 (cyclotomic witnesses), S599u (Eisenstein/`Cl₂(π/3)`),
HYP-2170/2201 (unit-distance lattice), Sawin/OpenAI disproof. New: **HYP-2257**.
