---
source: opus-2026-06-06-S699 (alternating group graph)
status: EXTENSION of the HN/LRC/UD unification to the non-abelian alternating group A_n. The alternating group graph AG_n is a forbidden-distance Cayley graph; NEW: its spectral ratio λmin/λmax = −1/2 for ALL n (λmin=−(n−2), λmax=2(n−2)) — the SAME ratio as the triangular lattice (the UD optimum), the triangle/Eisenstein signature ⟹ Hoffman χ≥3, α≤N/3. The full 3-cycle Cayley graph's spectrum = the CHARACTER RATIOS of A_n (the non-abelian Fourier transform = the Bessel/Dirichlet analog), verified A_5={20,0,0,5,−4} from irreps {1,3,3,4,5}. A_5 ≅ icosahedral group (geometric closure with HN). A_n = the parity group of tournament Ham-paths (Rédei/determinant).
tags: [alternating-group, AG_n, cayley-graph, representation-theory, character-ratio, hoffman-bound, eisenstein, triangle, icosahedral, hadwiger-nelson, parity, redei, unification]
---

# The alternating group graph: the non-abelian member of the distance-graph unification

**Prompt (user):** consider the alternating group graph and how it changes as n varies, regarding
the recent HN/LRC/UD work; extend and explore.

The alternating group graph is the **non-abelian member** of the forbidden-distance Cayley family
(S699g), and it is doubly apt: `A_n` is also the **parity group** behind Rédei/the determinant
(the repo's home turf). Two clean new findings.

## 1. AG_n is a forbidden-distance Cayley graph; its spectral ratio is −1/2 for all n

`AG_n = Cay(A_n, {(1\,2\,i),(1\,i\,2): 3≤i≤n})` — even permutations, adjacent iff they differ by
one of the special 3-cycles. Computed (`…s699h.py`, `n=3..6`):

```
   n | |A_n| | deg=2(n−2) | triangles | λmax | λmin | ratio | Hoffman χ≥ | α≤
   3 |   3  | 2          | 1         |  2   | −1   | −1/2  | 3         | 1
   4 |  12  | 4          | 8         |  4   | −2   | −1/2  | 3         | 4
   5 |  60  | 6          | 60        |  6   | −3   | −1/2  | 3         | 20
   6 | 360  | 8          | —         |  8   | −4   | −1/2  | 3         | 120
```

> **`λmax = 2(n−2)`, `λmin = −(n−2)`, so `λmin/λmax = −1/2` for *every* `n`.** Hence Hoffman gives
> `χ(AG_n) ≥ 1 − λmax/λmin = 3` and `α(AG_n) ≤ N/3` for all `n`.

**This is exactly the triangular-lattice signature.** In S699g the triangular (Eisenstein)
lattice — the unit-distance optimum — has `λmin=−3, λmax=6`, ratio `−1/2`, Hoffman `χ≥3`. `AG_n`
has the *same ratio* because it is **triangle-rich**: the 3-cycle generators `(1\,2\,i),(1\,i\,2)`
make `e,(1\,2\,i),(1\,i\,2)` a triangle (a `K_3`), and `K_3` has eigenvalues `2,−1,−1` (ratio
`−1/2`). **The Eisenstein/π/3 structure of the HN/UD unification reappears in the alternating
group graph as the constant `−1/2` eigenvalue ratio** — the spectral fingerprint of triangle
(3-cycle / equilateral) richness, shared by the triangular lattice (geometry) and `AG_n`
(group). The growth with `n`: degree `2(n−2)`, ratio fixed at `−1/2`, `χ ∈ [3, 2n−3]`
(Hoffman ↔ greedy), `α ≤ N/3`.

## 2. The full 3-cycle Cayley graph: spectrum = representation theory (the non-abelian Fourier transform)

Take all 3-cycles as generators: `Cay(A_n, \{3\text{-cycles}\})`, `g\!\sim\!h` iff `gh^{-1}` is a
3-cycle. Its eigenvalues are the **character ratios**: for a conjugacy class `C`,
```
   eigenvalues = { |C| · χ_ρ(C)/dim(ρ) : ρ ∈ Irr(A_n) },  multiplicity dim(ρ)².
```
**Verified `A_5`** (`…s699h.py`): spectrum `{20, 0, 0, 5, −4}`. Irreps of `A_5` have dims
`1,3,3,4,5`; the 3-cycle class has `|C|=20` and character values `1,0,0,1,−1`, so the eigenvalues
are `20·\{1, 0, 0, 1/4, −1/5\} = \{20, 0, 0, 5, −4\}` ✓ (multiplicities `1,9,9,16,25`, summing to
`60`). `λmin = −4` comes from the **5-dimensional irrep** (the one with negative character ratio),
giving Hoffman `χ(A_5\text{-3cyc}) ≥ 1 − 20/(−4) = 6`. Across `n` the Hoffman bound is `3,3,6,6`
for `n=3,4,5,6`.

> **This is the exact non-abelian analog of the HN/LRC/UD spectral unification (S699g):** there the
> chromatic/independence bounds came from the Fourier transform of the forbidden-distance measure
> (the **Bessel `J_0`** on `ℝ²`, the **Dirichlet kernel** on the circle, the **structure factor**
> on the lattice). For `A_n` the "Fourier transform" is the **representation theory**, and the
> "Bessel function" is the **character ratio `χ_ρ(\text{3-cycle})/\dim ρ`**. The smallest one (the
> top irreps) sets the chromatic floor. *Same theorem, abelian → non-abelian: Fourier becomes
> characters.*

## 3. Two closures: geometry (icosahedral) and parity (Rédei)

- **Geometric closure.** `A_5 ≅` the **icosahedral rotation group**. Its `60`-vertex 3-cycle
  Cayley graph lives on the icosahedron/dodecahedron — back to the **sphere**, the natural
  curved-space sibling of the Hadwiger–Nelson **plane**. The unification's "forbidden-distance
  Cayley graph" on `A_5` *is* a spherical distance graph. (And `A_4 ≅` tetrahedral, `S_4 ≅`
  octahedral — the Platonic ladder.)
- **Parity closure (the repo's home).** `A_n` = the **even permutations** = the *alternating* =
  the **parity** structure of Rédei/OCF. A tournament's Hamiltonian paths are permutations; the
  even ones form `A_n`. The **signed** Ham-path count (over `A_n` minus the odd coset) is the
  **determinant** side; the unsigned is the **permanent** (S599e, the `z=−1` slice). So the
  alternating group is literally the home of the repo's determinant/parity face — and `AG_n`'s
  triangle (3-cycle) richness is the OCF's odd-3-cycle generators.

## 4. The unification, extended

> **HN / UD / LRC / Rédei are χ / |E| / α / parity of forbidden-distance Cayley graphs; the
> spectrum is the (abelian or non-abelian) Fourier transform — Bessel `J_0` (plane), Dirichlet
> kernel (circle), lattice structure factor (Eisenstein), and *character ratios* (`A_n`). The
> triangle/Eisenstein signature is the constant eigenvalue ratio `−1/2` (Hoffman `χ≥3`), shared by
> the triangular lattice (geometry) and `AG_n` (group). `A_5` = icosahedral closes the loop to the
> sphere; `A_n` = parity closes it to Rédei/the determinant.**

How it changes with `n`: `AG_n` keeps the `−1/2` ratio (triangle signature) but grows in degree
`2(n−2)` and `χ` (toward `2n−3`); the full 3-cycle graph's character-ratio spectrum spreads (more
irreps), raising the Hoffman `χ`-floor (`3,3,6,6,…`) — the chromatic content migrates into the
**top irreps** (large-dimension, negative character ratio), the non-abelian "high-frequency" modes.

## 5. Honest status

- **Verified:** `AG_n` spectra `λmax=2(n−2)`, `λmin=−(n−2)`, ratio `−1/2`, Hoffman `χ≥3`, `α≤N/3`
  (`n=3..6`); the full 3-cycle `A_5` spectrum `{20,0,0,5,−4}` = the character ratios (irreps
  `1,3,3,4,5`).
- **Established (standard, here mapped):** Cayley-graph spectrum = character ratios (Babai/Diaconis);
  Hoffman/ratio bounds; `A_5≅` icosahedral; `A_n` = even permutations = parity.
- **New (the synthesis):** `AG_n` shares the triangular-lattice `−1/2` ratio (the triangle/Eisenstein
  spectral fingerprint) — extending the HN/UD unification to a non-abelian group; the character
  ratio as the non-abelian Bessel; the icosahedral (sphere) and Rédei (parity) closures.
- **Open / directional:** the true `χ(AG_n)` (between `3` and `2n−3`); does the icosahedral `A_5`
  3-cycle graph's chromatic number connect to spherical-code / spherical-HN bounds? Is the
  character-ratio spectrum the right tool for the LRC tournament's dichromatic number (the
  permutation/parity Cayley graph)?

**Artifacts:** `04-computation/alternating_group_graph_s699h.py` (+`.out`). Builds on S699g
(HN/UD/LRC spectral unification), S699a (kissing/Eisenstein), THM-402 (dichromatic), Rédei/OCF,
S599e (`z=−1`/determinant), Babai (Cayley spectra). New: **HYP-2266**.
