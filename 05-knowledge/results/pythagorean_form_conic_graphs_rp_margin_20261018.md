# One quadratic form over Z, R and F_p: the Pythagorean tree, minimal models, Petersen / K5 / K6 / K3,3 / Coxeter, and the margin of the RP trick

**Session:** opus-2026-10-18-S44 (second part).

**Status by section.**
- §2 Berggren image: PROVED (elementary).
- §1, §2 checks: FINITE-EXACT / VERIFIED. The identifications are partly classical, and no novelty is claimed for the
  finite geometry.
- §3 margin lemma: PROVED (elementary) and sharp.
- Cross-domain readings: OBSERVATION / ANALOGY, as marked.
- Self-checked. The independent S44 audit predates this note.

**Owner prompt.**
> keep going as you were and also merge in related to https://github.com/ethanwsussman/Yang--Mills-LLM-experiment and
> the ternary tree of primitive pythagorean triples and the petersen graph, K_5, K_(3,3)

**What was read, and what is inherited.**
- **The Yang–Mills repository** (Ethan Sussman and Levent Alpoge: an LLM-assisted attempt to complete Balaban's
  renormalization-group construction; QED3/YM3 first, YM4 the goal). I read:
  - `README.md`;
  - `README.pdf` (Sussman's proposal for a working group, dated 2026-10-06);
  - `QED3_slop/YM_Rotation_Invariance_simplified.pdf` ("Rotation invariance from lattice symmetry and reflection
    positivity", about 20 pages).

  The two PDFs were fetched with the web tool, text-extracted in an isolated scratch folder, and nothing was executed.
  I did not read the 26 MB of `theorem_0_*.tex` sources, the YM4 zeroth draft, or the other QED3 PDFs. Nothing here
  assesses "Theorem 0".
- **Concurrent agents on the same prompt.**
  - [Codex graph-Pythagorean synthesis](graph_pythagorean_20261018_synthesis.md): Seidel radii (3,4,5) of
    (Petersen, K5, K3,3), the H_t family on the unit-gap spine, the snark isosceles geometry.
  - [clock_duality ternary ports](ternary_ports_golden_defect_20261018j.md): it independently replayed this session's
    conic dictionary and corrected its group-scope docstring.
  - [mac-mini unit-path / Kuratowski / Yang–Mills note](unit_path_ring_kuratowski_pythagorean_yang_mills_20261018.md):
    THM-4701 (no K5 or K3,3 in `G_q`; Clebsch/Petersen), Pell branches of the tree, and a README-only Yang–Mills
    reading. Its next obligation was to read the RP rotation PDF, which §3 does.
- Repo Berggren canon: THM-3333, THM-3334, THM-3357.
- S44 VOA results: THM-4695, THM-4710, HYP-9315 ([note](li_spectral_sequence_quantum_transfer_20261018.md)).

The organizing fact is that one quadratic form, `Q = x^2 + y^2 - z^2`, carries all three themes.
- **Over Z.** Its primitive null vectors are the Pythagorean triples, and Berggren's ternary tree is their orbit
  structure under the integral Lorentz group.
- **Over R.** Its null cone is the light cone. Its Lie algebra `so(2,1) = sl(2,R)` is the global conformal algebra
  `<L_(-1), L_0, L_1>` that acts on the singular support (THM-4710).
- **Over F_p.** Its conic organizes the Petersen graph, K5, K6, K3,3 and the Coxeter graph.

## 1. Over Z: the Berggren tree is a family tree of minimal models (FINITE-EXACT; elementary)

A primitive triple is `(q^2 - p^2, 2pq, q^2 + p^2) = (p + iq)^2`, with `q > p >= 1` coprime and `q - p` odd. Berggren's
matrices act on `(m, n) = (q, p)` by:
```
A: (m,n) -> (2m-n, m),     B: (m,n) -> (2m+n, m),     C: (m,n) -> (m+2n, n).
```
Attach to each node the Virasoro `(p,q)` model: a minimal model for `p >= 2`, and the logarithmic `(1,q)` triplet
series for `p = 1`. `berggren_minimal_models_20261018.py` checks the following exactly.
- **Bijection.** Every coprime pair `p < q <= 60` with `q - p` odd occurs exactly once.
- **Central charge from the triangle.** With `b = 2pq` the even leg and `h` the hypotenuse:
  ```
  c_(p,q) = 13 - 12 h / b,          c_eff = 1 - 6/(pq) = 1 - 12/b.
  ```
  Examples: Ising `(7, 24, 25)` gives `13 - 25/2 = 1/2`; Yang–Lee `(21, 20, 29)` gives `-22/5`.
- **Möbius moves.** On `t = q/p` the moves are `t -> 2 - 1/t`, `2 + 1/t`, `t + 2`, and `c = 13 - 6(t + 1/t)`.
- **Branches.**

| branch | models | `t -> ` | `c -> ` | known as |
|---|---|---|---|---|
| A, A, A, ... from the root | unitary `M(n, n+1)` | 1 | 1 | Codex's unit-gap spine (`h - b = 1`) |
| B, B, B, ... | `(1,2), (2,5), (5,12), (12,29), ...` (Pell) | `1 + sqrt 2` | `13 - 12 sqrt 2` | Codex's isosceles branch; mac-mini's `sqrt 2` Pell branch |
| C, A, C, A, ... | `(1,2), (1,4), (4,7), (4,15), (15,26), ...` | `sqrt 3` | `13 - 8 sqrt 3` | mac-mini's near-60-degree branch |
| C, C, ... from `(3,2)` and `(5,2)` | `(2, 2N+1)` | `∞` | `-∞` | the classically free series |

- **C fixes `p`.** So C-chains are lines of fixed `p`, the direction of the q-stability of transfer orders
  (HYP-9315(f)). Stability in fact holds across all `q` coprime to `p`, not only along one C-chain.
- **Index coincidence (NUMEROLOGY).** Codex's graphs `H_t` (Seidel radius `t`, order `t^2 + 1`) sit on the unit-gap
  spine, so `t = 2n + 1` matches `M(n, n+1)`. Petersen `= H_3` matches the root `(1,2)` (`c = -2`, symplectic
  fermions), and `H_7` matches Ising. No structural map is claimed.
- **Frames (OBSERVATION).** In the triple's angle `arctan(b/a)`, the A-chain tends to `π/2`, the C-chain to `0`, and the
  B-chain to `π/4`. These are the coordinate and diagonal reflection lines of `Z^2`, i.e. the two frames of the RP
  trick (§3). The unitary chain converges to `c = 1`, the free boson, along a coordinate frame.

## 2. Over F_p: the Pythagorean conic and the small symmetric graphs (FINITE-EXACT; partly classical)

In `PG(2,p)` take the conic `C: x^2 + y^2 = z^2`, the reduction of the Pythagorean light cone. A non-conic point is
*exterior* if its polar meets `C` in 2 points and *interior* if in 0. `G_int(p)` and `G_ext(p)` are the conjugacy graphs
(`B(P,Q) = 0`). The script is `pythagorean_conic_dictionary_20261018.py`.

**p = 5** (6 conic, 15 exterior, 10 interior points).

| conic object | graph object |
|---|---|
| `G_int(5)` | the **Petersen graph** |
| `G_ext(5)` | `5 K3`: the 5 self-polar exterior triangles, which are the vertices of **K5** |
| the 6 conic points | the vertices of **K6**; also the 6 perfect matchings of Petersen and the 6 pentagons of the hemi-dodecahedron |
| exterior point | edge of K6: its tangency duad; also an edge of Petersen, the third vertex of the self-polar triangle on it |
| interior point | a syntheme of K6 (10 of 15, the other 5 from the exterior triangles); equivariantly, a 3+3 split of the conic points |
| self-polar triangles | 20 = 5 (three exterior) + 15 (one exterior, two interior) |

**Kuratowski views of Petersen.**
- *Contracting* the perfect matching `M_x` gives **K5**, whose vertices are the 5 conic points other than `x`. Here
  `M_x` is the set of 5 edges labelled by the exterior points on the tangent at `x`. Verified for all 6 matchings.
- *Deleting* an interior point `P` gives a subdivided **K3,3**. But it is *not* equivariantly the conic split
  `T × T'` attached to `P`: `Stab(P) = S_3` acts regularly on the 6 far vertices and with two fixed orbits on the conic
  points. This is a negative result, recorded as such.

**p = 7** (8 conic, 28 exterior, 21 interior points).
- `G_ext(7)` is the **Coxeter graph** (girth 7). It was checked against an independent construction, `KG(7,3)` minus
  the 7 lines of a Fano plane.
- `G_int(7)` is the line graph of the **Heawood graph**.

**Survey through p = 13.**
- `G_int(p)` has degree `(p+1)/2` and `G_ext(p)` has degree `(p-1)/2`. So the only cubic cases are `G_int(5)` =
  Petersen and `G_ext(7)` = Coxeter.
- These are exactly the two famous connected vertex-transitive graphs without a Hamiltonian cycle, apart from
  `K_1`, `K_2` and the truncations of these two.
- Every other connected conic graph for `p <= 13` was found Hamiltonian, by explicit cycle: `p = 7` interior, and
  `p = 11, 13` both.
- Petersen's non-Hamiltonicity was re-checked exhaustively. Coxeter's is classical.
- Triangles: for `p ≡ 1 (mod 4)`, `G_int` is triangle-free and `G_ext` has triangles; for `p ≡ 3 (mod 4)` it is the
  reverse.

**Literature.**
- The PG(2,5) internal-point structure is classical. Brouwer's table describes the Desargues graph, the bipartite
  double of Petersen, as the configuration of 10 internal points and 10 elliptic lines of a conic in PG(2,5).
- Brouwer's Coxeter page uses Fano antiflags (Coxeter 1983; BCN section 12.3). The conic-exterior form was not located
  in the sources checked, but it is surely known in some form.
- No novelty is claimed for the finite geometry. The point here is that the Pythagorean form is their common home.

**Proposition (the tree mod p; PROVED, elementary).** Let `p` be an odd prime. The Berggren matrices generate in
`PGL(2,p) ≅ SO(3, F_p)` the subgroup `PSL(2,p)` if `p ≡ 1 (mod 4)`, and all of `PGL(2,p)` if `p ≡ 3 (mod 4)`.

*Proof.*
- The Veronese map `(m,n) -> (m^2 - n^2, 2mn, m^2 + n^2)` intertwines `A, B, C` with
  `A' = [[2,-1],[1,0]]`, `B' = [[2,1],[1,0]]`, `C' = [[1,2],[0,1]]`. This was checked directly.
- `C'^(-1) A' = S = [[0,-1],[1,0]]`, so `S C' S^(-1) = [[1,0],[-2,1]]`. Since 2 is invertible mod `p`, the two
  elementary matrices with entry `±2` generate `SL(2, F_p)`.
- `det B' = -1` is a square mod `p` exactly when `p ≡ 1 (mod 4)`. ∎

Data: image orders 24 (p=3), 60 (p=5), 336 (p=7), 1320 (p=11), 1092 (p=13).
- At `p = 5` the image is `A_5`. This is the rotation group of the dodecahedron modulo the antipodal map, i.e. the
  automorphism group of the hemi-dodecahedral map whose skeleton is Petersen. It has index 2 in `Aut(Petersen) = S_5`.
  The full `S_5` would need orientation reversal.
- clock_duality's matrix-closure audit confirms 60.

**Cross-links (CITED / ANALOGY).**
- **Modular images.** The modular representation of the Yang–Lee model `M(2,5)` has projective image
  `PSL(2,5) = A_5`, the icosahedral group: `T` has eigenvalue ratio of order 5, so the image is the (2,3,5) triangle
  group. For `M(2,7)` it is `PSL(2,7)` (Klein). These are the groups through which the Pythagorean tree acts on the
  Petersen and Coxeter conic graphs. The (2,5) and (2,7) models are also the first nodes of the classically free C-chains
  of §1. This is a coincidence of groups; no map is claimed.
- **Conway–Gordon–Sachs** (classical).
  - Label the vertices of K6 by the 6 conic points. The canonical equivariant bijection, interior point to 3+3 split,
    makes the 10 pairs of disjoint triangles of K6 the vertices of the Petersen graph.
  - CGS: in every spatial embedding of K6, the sum over these 10 pairs of the linking numbers is odd.
  - In abelian Chern–Simons language, with charges such that `<W_T W_T'> = (-1)^lk`, the product over the Petersen
    vertices of `<W(T_v) W(T'_v)>` is `-1` for every embedding (ANALOGY).
  - The Petersen family, the Y–Δ class of K6, is the set of forbidden minors for linkless embedding
    (Robertson–Seymour–Thomas). There K5 and K3,3 enter only as Kuratowski's planarity obstructions.

## 3. Over R: the light cone and the margin of the RP trick (PROVED, elementary; reading of the repository's argument)

**The repository's Theorem A**, from `YM_Rotation_Invariance_simplified.pdf`, paraphrased. Consider a law `μ` of closed
two-forms on `R^d`, `d >= 3`, with:
- Gaussian domination;
- invariance under translations and signed coordinate permutations;
- reflection positivity across a coordinate hyperplane `{t = 0}` and a diagonal hyperplane `{t + x = 0}`.

Then `μ` is invariant under all Euclidean isometries.

The proof has four steps.
1. The rotation defect `D[φ](θ) = S_N(R(θ)φ) - S_N(φ)` in the tx-plane has Fourier modes `ℓ ∈ 4Z`, because quarter
   turns are symmetries.
2. The two reflections give the spectral inequality `H >= |P|`. Joint spectrum in the forward light cone
   `{E >= |p|}` is the real Pythagorean cone. Together with Gaussian domination it gives an analytic continuation to
   imaginary angles with growth `e^(2N|η|)`, so the maximal surviving frequency satisfies `M_N <= 2N`.
3. Cauchy–Schwarz on the RP inner products gives `2M_(a+b) <= M_(2a) + M_(2b)`. Positivity of the leading coefficient
   in both frames forces `M_(2a) ∈ 8Z`, because the `π/4` frame change multiplies mode `ℓ` by `e^(iℓπ/4) = (-1)^(ℓ/4)`.
4. Convexity, integrality and slope then kill all frequencies.

**Lemma (the final step, abstracted; PROVED).** Let `M_n >= 0` be integers with
- (i) `M_0 = M_2 = 0`;
- (ii) `2 M_(a+b) <= M_(2a) + M_(2b)`;
- (iii) `M_(2a) ∈ μZ`;
- (iv) `M_n <= σ n`.

If `μ > 2σ`, then `M_n = 0` for all `n`.

*Proof.*
- By (ii) with `(a,b) = (j-1, j+1)`, the increments `Δ_j = M_(2j+2) - M_(2j)` are nondecreasing. They lie in `μZ` by
  (iii), and `Δ_0 = 0`.
- Suppose some `Δ_(j0) > 0`. Then `Δ_j >= μ` for all `j >= j0`, so `M_(2j) >= μ(j - j0)`. This contradicts
  `M_(2j) <= 2σ j` once `μ > 2σ`.
- So `M_(2j) = 0`, and (ii) with `a = 0` gives `2M_N <= M_0 + M_(2N) = 0`. ∎

**Sharpness.** At `μ = 2σ` the lemma fails: `μ = 8, σ = 4` with `M_n = max(0, 4n - 8)` satisfies (i)-(iv) and is
nonzero. Only the abstracted step is sharp; the full argument may have more structure.

**Margin.**
- In the theorem, `μ = 2m`, where `m = 4` is the order of the lattice rotation group and the second frame is at angle
  `π/m`.
- `σ = 1 + w`: one unit is the planar position cost, and `w` is the largest SO(2)-weight of the field's index
  representation. For two-forms `w = 1`, matching the note's `e^(2N|η|)`.
- So the step closes iff **`m > 1 + w`**. On square or hypercubic lattices (`m = 4`) it covers `w <= 2`: two-forms with
  margin, and symmetric 2-tensors such as a stress tensor.
- It does not by itself cover `w >= 3`, for instance spin-3 currents. A lattice with `m = 6` would cover `w <= 4`.
- This is a robustness remark about the argument, not a critique: the note's case `w = 1` has a factor-2 margin.

**The Pythagorean tree here.**
- The rational rotations of `Z^2` that map the lattice onto a coincidence sublattice are exactly the Pythagorean
  rotations. This is classical coincidence-site-lattice theory: the coincidence index is the hypotenuse.
- The tree enumerates them modulo the point group.
- A lattice law has no reflection positivity across such lines. If it had positivity in a frame at an angle
  `θ = arctan(4/3)`, which is irrational relative to `π` by Niven, then `e^(iℓθ) ≠ ±1` for every `ℓ ≠ 0` and the sign
  comparison alone would kill every frequency.
- The convexity–integrality step is the price of having only lattice frames.

**The so(2,1) thread (OBSERVATION).** The same Lorentz algebra appears three times:
- as Berggren's integral Lorentz group;
- in the light-cone spectrum condition `H >= |P|`;
- as the global conformal `sl_2` whose lowest-weight vectors are the singular-support generators (THM-4710).

In Theorem A, analytic continuation to imaginary angles turns rotations `R(iη)` into boosts. That is the passage from
`SO(2)` to `SO(1,1)` inside this form.

## 4. Ledger and next obligations

| claim | type |
|---|---|
| Berggren tree = minimal models with `pq` even; `c = 13 - 12h/b`; branch limits | FINITE-EXACT (elementary) |
| `G_int(5)` = Petersen, `G_ext(5)` = 5K3, K6/synthemes/K5-contraction views; `G_ext(7)` = Coxeter, `G_int(7)` = L(Heawood) | FINITE-EXACT (partly classical) |
| deletion K3,3 is not the conic split | FINITE-EXACT (negative) |
| Berggren image `PSL(2,p)` / `PGL(2,p)` by `p mod 4` | PROVED |
| cubic conic graphs = the two non-Hamiltonian vertex-transitive graphs; others Hamiltonian for `p <= 13` | OBSERVATION / FINITE |
| RP-trick final step: lemma and sharp margin `m > 1 + w` | PROVED (elementary) |
| frames, modular images, CGS/Wilson loops, so(2,1) thread | ANALOGY / OBSERVATION |

**Next obligations.**
1. With the owner's permission, read the Theorem 0 sources and the QED3 manuscripts. Check whether the field
   representations used ever have `w >= 3`. If so, the margin matters.
2. Is the conic-exterior construction of the Coxeter graph in the literature?
3. The minimal models along the √2-Pell and √3-Pell branches: do their transfer-order tables (HYP-9315) show any
   regularity in `p` beyond universality?
4. `G_int(p)` and `G_ext(p)` for larger `p`: Hamiltonicity (Lovász's conjecture), and the spectra as Kloosterman-type
   sums. Compare THM-4628 and THM-4701.
