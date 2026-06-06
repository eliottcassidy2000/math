---
source: opus-2026-06-06-S699 (HN ∪ LRC ∪ UD unification)
status: UNIFICATION — Hadwiger–Nelson, the unit-distance problem, and the LRC are three invariants (chromatic number χ / edge density |E| / independence density α) of the SAME forbidden-distance Cayley graph. Master quantity: the independence density α, with χ_f = 1/α. Shared technique: the Fourier transform of the forbidden-distance measure (Bessel J0 for the plane, Dirichlet kernel for LRC arcs, lattice structure factor for UD) via the Hoffman/Lovász spectral bound (verified: triangular χ≥3 tight, plane χ≥3.48). Shared structure: Eisenstein/π/3 (Cl2(π/3)). Shared hardness: incommensurate rotations = LRC Diophantine resonance (Moser spindle θ=33.56°).
tags: [hadwiger-nelson, unit-distance, LRC, distance-graph, cayley-graph, chromatic-number, independence-density, fractional-chromatic, hoffman-bound, bessel, eisenstein, incommensurate-rotation, unification]
---

# Hadwiger–Nelson, unit distance, LRC: one graph, three invariants

**Prompt (user):** pursue the Hadwiger–Nelson / LRC / unit-distance unification; how the problems
are the same underlying thing; how insights of one are keys to the other.

They are the **three classical invariants of one object** — the *forbidden-distance Cayley graph* —
and the master quantity is its **independence density**.

## 1. One graph, three invariants

Let `G_d(X)` be the **distance graph**: vertices = points of a space `X` (with a group structure),
edges = pairs at the forbidden distance `d` (a Cayley graph `Cay(X, {‖·‖=d})`). The three problems
ask the three standard invariants of `G_d`:

| problem | the invariant of `G_d` | value |
|---|---|---|
| **Hadwiger–Nelson** | **chromatic number `χ(G_d(ℝ²))`** | `∈ [5,7]` |
| **unit distance** | **edge density** (max `|E|` on `n` vertices) | `≈κn/2` (lattice) … `n^{1.014}` (CM) |
| **LRC** | **independence / covering density `α`** of `G_d(circle)`; its tournament's **dichromatic** `χ` | `α = p_0`; `χ = 2` (worry-set, THM-402) |

> They are linked by one identity: for a vertex-transitive distance graph,
> ```
>        χ_f(G_d) = 1 / α(G_d),     χ ≥ χ_f,
> ```
> so the **independence density `α` is the master quantity**: HN's `m_1` (plane 1-avoiding density,
> `≈0.229 ⟹ χ_f≈4.36`), LRC's `p_0` (the lonely-set / safe density), and UD's edge density (which
> *bounds* `α` from above via Turán/kissing) are all `α` of the same graph family in different
> ambient groups (`ℝ²`, the circle, a lattice). **The lonely set of LRC is an independent set of a
> distance graph; the chromatic number of HN is `1/`(its density). Same `α`.**

## 2. Shared technique: the Fourier transform of the forbidden-distance measure

The forbidden-distance set is a measure `μ_d` (the sphere `‖·‖=d`); its **Fourier transform** is the
graph's "symbol", and the **Hoffman/Lovász spectral bound** turns it into a `χ`/`α` bound. Verified
(`…s699g.py`):
- **Triangular lattice (UD/Eisenstein):** symbol `2(cos a+cos b+cos(a+b))`, range `[-3,6]` ⟹
  Hoffman `χ ≥ 1 − 6/(−3) = 3` — **tight** (`χ=3`).
- **Plane (HN):** the unit-circle measure transforms to the **Bessel `J_0`**; `min J_0 = −0.4028`
  (at `x≈3.83`) ⟹ spectral `χ ≥ 1 − 1/(−0.4028) = 3.48` (the same method; the true `χ≥5` needs the
  combinatorial Moser/de Grey graphs).
- **LRC (circle):** the danger arc `‖vt‖<δ` transforms to the **Dirichlet/Fejér kernel**; the
  covering-depth moments (THM-406) are its autocorrelations; `p_0 = Σ(−1)^{|S|}μ(∩D_i)` is the
  inclusion–exclusion of the same measure.

> **One method:** the chromatic/independence number of all three is read from the spectrum of the
> forbidden-distance measure — Bessel (plane), structure factor (lattice), Dirichlet kernel
> (circle). HN's `J_0`-bound, UD's Turán/kissing bound, and LRC's first-moment/`p_0` bound are the
> same Fourier–spectral computation.

## 3. Shared structure: Eisenstein ζ₆ / π/3 (and `Cl₂(π/3)`)

All three problems' **extremal/hard objects live on the Eisenstein lattice (angle π/3)**:
- **UD optimum** = the triangular = Eisenstein lattice (`κ=6`, S699a);
- **LRC worry-set** witnesses = roots of unity (THM-403); `n=14`'s prime-3 (`C=27=3³`) = the
  Eisenstein norm-3 ideal (HYP-2170), and `Cl₂(π/3)=1.0149` is the shared constant (S599u);
- **HN's chromatic-forcing graphs** (Moser spindle `χ=4`, de Grey `χ=5`) are built on the Eisenstein
  lattice. The triangular-lattice 3-coloring `(i−j) mod 3` (the Eisenstein norm-3 ideal) is exactly
  the LRC `n=14` prime-3 structure.

## 4. Shared hardness: incommensurate rotations = LRC Diophantine resonance

HN's chromatic *lower* bounds (>3, the hard direction) come from **incommensurate rotations** of the
Eisenstein lattice: the Moser spindle uses `θ = 2·arcsin(1/(2√3)) = 33.56°` — an *irrational*
multiple of `π`. De Grey's graph stacks many such rotations. **This irrationality/incommensurability
is precisely the LRC's Diophantine resonance structure** (the speeds' irrational ratios, the
two-block `2^E−3^k` gaps, linear forms in logs). So:

> **The Diophantine machinery of LRC (resonance gaps, linear forms in logarithms, the cyclotomic
> shell calculus) is the natural language for HN's incommensurate-rotation constructions, and the
> Eisenstein-lattice combinatorics of UD/HN is the natural language for the LRC worry-set.** Hard in
> one = hard in the other, for the *same arithmetic reason* (incommensurability over `ℤ[ζ₆]`).

## 5. The keys (insights of one as keys to the other)

- **HN → LRC.** The **measurable/fractional chromatic density** technology (Falconer's `m_1≥1/5`
  via Fourier/autocorrelation; the Lovász theta spectral bound) is the right tool to bound the LRC
  **lonely-set measure `p_0`** — both are the independence density `α` of a distance graph, bounded
  by the forbidden-distance measure's transform. The "Vitali wall" (LRC measure-blindness on the
  worry-set, THM-406 M2) is the HN phenomenon that `χ > χ_f` (the integrality gap between the
  *fractional* density bound and the *true* chromatic number — the combinatorial Moser/de Grey jump).
- **LRC → HN.** The **cyclotomic worry-set** (THM-403), the **doubling-rigidity** (THM-404), and the
  `2n−1` **shell calculus** are the arithmetic of `ℤ[ζ₆]`/`ℤ[ζ_n]` — the home of the Eisenstein
  constructions HN needs. The LRC's **covering-depth inclusion–exclusion** (THM-406) is the density
  method for HN's independence ratio.
- **UD → both.** The **kissing-number cap** `κ≤6` (S699a) is a *degree* bound, hence a chromatic
  *lower* bound (`χ ≥ ` clique/degree facts) and an independence *upper* bound (Turán) — UD's edge
  density directly constrains both HN's `χ` and LRC's `α`.

## 6. The one-sentence unification

> **Hadwiger–Nelson, the unit-distance problem, and the lonely-runner conjecture are the chromatic
> number, the edge density, and the independence density of one forbidden-distance Cayley graph;
> all three are read from the Fourier transform of the forbidden-distance measure, attain their
> extremes on the Eisenstein (π/3) lattice, and are hard for the same Diophantine
> (incommensurate-rotation / resonance) reason.** The integrality gap `χ > χ_f = 1/α` *is* the LRC
> "Vitali wall" — the gap between the measurable density and the true combinatorial answer.

## 7. Honest status

- **Verified:** Hoffman `χ≥3` (tight) for the triangular lattice; Bessel `J_0`-bound `χ≥3.48` for
  the plane; the Eisenstein/π/3 shared structure; the Moser incommensurate rotation `33.56°`.
- **Established (standard math, here mapped):** `χ_f = 1/α`, `χ≥χ_f`; the Hoffman/Lovász and
  Fourier-transform (`J_0`/Dirichlet) bounds; HN=`χ`, UD=`|E|`, LRC=`α` of distance graphs.
- **New (the synthesis):** the three-invariants-of-one-graph framing with `α` as the master
  quantity; the `χ>χ_f` integrality gap = the LRC Vitali wall; the incommensurate-rotation = LRC
  Diophantine transfer; the Eisenstein/π/3 anchor across all three (building on S599u, HYP-2170).
- **Directional (the keys, not all proven):** HN density methods → LRC `p_0`; LRC cyclotomic/shell
  calculus → HN Eisenstein constructions; UD kissing → HN/LRC bounds.

**Artifacts:** `04-computation/hn_lrc_ud_unification_s699g.py` (+`.out`). Builds on THM-402/403/404
(dichromatic/cyclotomic/doubling), THM-406 (covering-depth/Vitali), S599u (`Cl₂(π/3)`/π/3), S699a
(kissing), HYP-2170 (Eisenstein χ=3), Falconer/Lovász/Hoffman, de Grey/Moser. New: **HYP-2264**.
