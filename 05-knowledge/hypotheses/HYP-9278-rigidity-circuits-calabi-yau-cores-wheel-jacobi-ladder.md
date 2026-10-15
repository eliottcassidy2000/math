---
id: HYP-9278
title: "Paley pattern sums climb a modular ladder at 2D rigidity circuits (the count |E| = 2|V| - 2 makes the arrangement double cover canonically trivial; new motives of weight |V| - 2). Wheels are Greene hypergeometric values (Jacobi power sums); W_4, W_6 are the Ahlgren-Ono / Frechette-Ono-Papanikolas level-8 forms; W_8 is not of level 32. Planar Fourier duality (classical, Kramers-Wannier/Wu-Wang) fixes the canonical-class line and swaps its sides: the cube equals q times its dual octahedron."
status: >
  OPEN (the dictionary as a general statement).
  PROVED (standard): (1) wheels as power sums of Jacobi sums, Tr((DS)^n) = sum_psi (chi(-1) J(psi, chi))^n
  = (q-1) q^{n-1} nF_{n-1}(1) (Greene); (2) planar Fourier duality U(G) = chi(-1)^{|E|/2} q^{(|V| - |V*|)/2} U(G*) for
  connected plane graphs with |E| even (a specialisation of Kramers-Wannier / Wu-Wang 1976); (3) the canonical-class count
  (the double cover of P^{|V|-2} branched along |E| hyperplanes has K = pi^* O(|E|/2 - |V| + 1); "Calabi-Yau" for these
  singular covers presumes a resolution) and, for plane graphs, |E| = 2|V| - 2 iff |F| = |V|.
  CITED: T_4 (Ahlgren-Ono 2000), T_6 = (q-1)(q^2 - q a_q - b6_q) (Koike, proved by Frechette-Ono-Papanikolas, IMRN 2004);
  rechecked at 39 primes.
  FINITE-EXACT: T_8 is not in the span of S_8(Gamma_0(32)) plus twisted lower weights (39 primes, 35 functions, rank
  17 vs augmented 18); U(Q_3) = q U(K_{2,2,2}) at 32 primes (audit: 35); the circuit-reduction rows; duality at 7 graphs
  (audit: rotation-system embeddings, signs checked). AUDITED: MISTAKE-613.
  REFUTED (pre-registered in this session): "Laman-sparse pure cores are Tate (polynomial)" - the cube Q_3 is not
  polynomial; it carries its dual's W_4 form.
session: opus-2026-10-14-S35
source: 05-knowledge/results/paley_least_resistance_runs_one_modular_form_20261014.md (sections 4b-4d)
scripts:
  - 04-computation/experiments/paley_wheel_ladder_20261014.py
  - 04-computation/experiments/paley_wheel_ladder_fit_20261014.py
  - 04-computation/experiments/paley_wheel_jacobi_20261014.py
  - 04-computation/experiments/paley_rigidity_dictionary_tests_20261014.py
  - 04-computation/experiments/paley_cube_fano_core_20261014.py
  - 04-computation/experiments/paley_planar_fourier_duality_20261014.py
related: [THM-4652, THM-4653, HYP-9274, HYP-9275, HYP-9276]
---

# HYP-9278: rigidity circuits are the Calabi–Yau cores of Paley pattern sums

Notation: `U(G) = Σ_{x: V → F_q} Π_{(u→v) ∈ E} χ(x_v − x_u)`, with `χ` the quadratic character. `q ≡ 3 (mod 4)` gives Paley tournaments and `q ≡ 1 (mod 4)` Paley graphs.

## Proved ingredients

1. **Calabi–Yau count.** `U(G)` is, up to trivial terms, a point count of the double cover of `P^{|V|−2}` branched along the hyperplanes `x_u = x_v`. A double cover of `P^n` branched in degree `2d` has `K = π^*O(d − n − 1)`.
   - It is Calabi–Yau iff `|E| = 2|V| − 2`.
   - It is Fano-type below that count and of general type above it.
   - Odd `|E|` gives `U = 0`.
   - `|E| = 2|V| − 2` with all proper subgraphs Laman-sparse is a 2D rigidity circuit.
   - For plane graphs (Euler's formula) `|E| = 2|V| − 2` iff `|F| = |V|`, so every self-dual plane graph has the Calabi–Yau count.
2. **Wheels are power sums of Jacobi sums.** `U(W_n) = q·Tr((DS)^n)`, and `Tr((DS)^n) = Σ_ψ (χ(−1)J(ψ, χ))^n`.
   - The proof substitutes `t_i = x_{i+1}/x_i` and uses Fourier analysis on `F_q^×`.
   - The nonzero spectrum of `DS` is the Jacobi-sum family `{χ(−1)J(ψ, χ)}`.
3. **Planar Fourier (Kramers–Wannier) duality.** For a connected plane graph with `|E|` even and the induced dual orientation,

   `U(G) = χ(−1)^{|E|/2}·q^{(|V| − |V*|)/2}·U(G*)`.

   - The proof is a Fourier expansion: the vertex sum becomes a flow sum, and flows are face potentials modulo constants.
   - The key input is `χ̂ = χ(−1)G·χ` with `G² = χ(−1)q`: the quadratic character is a Fourier eigenfunction.
   - Duality preserves the Calabi–Yau count and exchanges Fano-type (`|E| < 2|V| − 2`) with general type.

## The ladder (FINITE-EXACT)

| graph | `|V|`, `|E|` | type | arithmetic |
|---|---|---|---|
| Y-junction | — | elliptic double cover | Legendre trace, weight 2 |
| `K_4 = W_3` (self-dual) | 4, 6 | circuit | 0 for tournaments; `η(4z)⁶`, weight 3, for Paley graphs |
| `W_4` (self-dual) | 5, 8 | circuit | `T_4 = −(q−1)(q + a_q)`, `a = η(2z)⁴η(4z)⁴`, weight 4, level 8 |
| `K_6 − P_6` (self-dual, planar) | 6, 10 | circuit | `Sym²` of the weight-3 level-32 newform (mac-mini's HYP-9276) |
| `K_(2,4) + 2K_2` | 6, 10 | circuit | `q(q−1)b_q²` for Paley graphs (HYP-9275 B) |
| `W_6` (self-dual) | 7, 12 | circuit | `T_6 = (q−1)(q² − q a_q − b6_q)`, `b6` the weight-6 level-8 newform `q + 20q³ − 74q⁵ − 24q⁷ + …` (Frechette–Ono–Papanikolas) |
| `W_8` | 9, 16 | circuit | `₈F₇(1)`: not in the span of `S_8(Γ_0(32))` plus twisted lower weights; no modular evaluation known; OPEN |
| `K_4` plus a degree-2 vertex | 5, 8 | contains `K_4` | `−b_q`, the `K_4` form |
| octahedron minus two edges at a vertex | 6, 10 | contains `W_4` | `q + a_q`, the `W_4` form |
| octahedron `K_(2,2,2)` | 6, 12 | general type | `(5 − q)a_q + 3q`, the `W_4` form twisted |
| cube `Q_3` (dual of the octahedron), one of the five connected cubic graphs on 8 vertices, the smallest pure cores in the Fano regime | 8, 12 | Fano, Laman-sparse pure core | `q·((5 − q)a_q + 3q)` = `q` × octahedron |
| `K_(3,3) − e`, prism minus a rung | 6, 8 | Fano, reducible by THM-4652 after series steps | polynomial |

## The hypothesis

For every graph `G`:
1. The non-Tate part of `U(G)` is generated by modular motives attached to rigidity circuits contained in `G`, or, when `G` is planar, in its dual `G*`. The duality case is forced by part 3 of the proved ingredients.
2. A circuit on `|V|` vertices contributes a motive of weight `|V| − 2`. For wheels this is the Ahlgren–Ono / Frechette–Ono–Papanikolas theory of `ₙF_{n−1}(1)` through `m = 3`. `W_8` leaves the level-32 trivial-character world.
3. Self-dual plane circuits (wheels, `K_6 − P_6`) are the carriers of genuinely new motives at each size.

**Refuted version.** The first formulation in this session, "Laman-sparse pure cores are Tate", predicted that the cube `Q_3` would be a polynomial. Pure cores in the Fano regime need `|V| >= 8`, and the cube is one of the five connected cubic graphs on 8 vertices. The prediction was registered before computing, and it failed. The repair is part 1's duality clause.

**Also refuted.** This file first predicted that `W_8` needs weight-8 forms of level 16 or 32. A test with more primes than functions shows `T_8` is not in the span of all trivial-character forms of level dividing 32.

## Open

- Prove the weight assignment for circuits beyond wheels.
- Identify `T_8 = ₈F₇(1)` (up to normalisation). Candidates: nebentypus `χ_8`, level 64, products with odd-weight forms.
- Test non-planar circuits: `K_(3,3) + e` vanishes for tournaments by symmetry, so test its Paley-graph side.
- (Answered by the literature, audit I1: the level-8 forms of weights 4 and 6 are Ahlgren–Ono's and Frechette–Ono–Papanikolas' evaluations of `₄F₃(1)` and `₆F₅(1)`.)
