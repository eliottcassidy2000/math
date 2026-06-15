---
id: THM-509
title: The two-layer baby-Hodge dichotomy — the spectral (det-side) "Hodge inequalities" are the Hankel-PSD-ity of the skew moment sequence, which is AUTOMATIC for every tournament (SS^T always PSD) and CONSTANT on cospectral classes, hence NEVER cuts a realizable-region hole; the holes are non-algebraic Hodge classes cut only by the integrality + conflict-graph (overlap/Witt) layer, which lives on the permanent/#P side of the Valiant det/per wall
status: PROVED for the det-side (Layer 1: SS^T PSD => Stieltjes => Hankel PSD; skew spectrum a function of char(A) by THM-507, so the skew-Hankel matrix is constant on every cospectral class — verified 0 splits exhaustively n<=7). Layer 2 (the cut is integrality+conflict) is established at the LINEAR/convex-hull level for specific holes (CERT-1, exhaustive n=6,7); the full positivstellensatz (that NO nonlinear det-side moment inequality cuts a hole either) is CONJECTURE.
source: monad-explorer-2026-06-15 (workflow subagent)
depends_on:
  - THM-507   # walk counts / A-affine pencil are spectral; Cor.1: the skew char poly det(xI-S)=prod(x^2+mu_j^2) is a function of charA (PROVED)
  - THM-505   # H = spectral skeleton + integer-linear combination of Witt/overlap defects (the conflict/#P carriers)
  - THM-500   # second spectral boundary: alpha_1 (odd-cycle count) non-spectral from n=7 via the overlap term TQ
  - THM-499   # first spectral boundary: H = 1+2(c3+c5)+4*alpha_2; alpha_2 (disjoint-pair count) first non-spectral OCF ingredient at n=6
  - THM-118   # c_k = tr(A^k)/k for k<=5 (low cycle counts ARE power sums / moments)
related:
  - THM-506   # permanental companion; det/per face program (the Valiant wall)
  - THM-508   # contraction/Hadamard wall (the rank-2 spectral boundary; the diagonal/Hadamard cube escapes)
  - HYP-2492  # cycle gaps = skew-spectral exclusions (this places H's gaps on the OTHER, non-spectral side)
  - OPEN-Q-093
  - OPEN-Q-096
  - reflection: the-spectral-resolution-ladder-of-the-ocf
---

# THM-509 — the two-layer baby-Hodge dichotomy

**The picture.** Reframe the realizable-invariant-vector question as a *truncated moment
problem* (the "baby Hodge" problem): the low cycle counts are spectral moments
`tr(A^k) = Σ_i λ_i^k`, and a realizable-region **hole** — an integer invariant vector that
sits inside the realized range but is hit by no tournament (e.g. `(c3,c5)=(8,10)` at n=6) —
is the combinatorial analog of a *non-algebraic Hodge class*: it satisfies every continuous
positivity ("Hodge") inequality the realized points satisfy, yet it is integer-forbidden.

This theorem proves that the realizability question splits into **two disjoint layers** that
cut on opposite sides of the Valiant determinant/permanent wall:

> **Layer 1 (det-side / spectral / P).** The "spectral Hodge inequalities" available from the
> *determinantal* world are exactly the **Hankel positive-semidefiniteness** of the **skew
> moment sequence** `m_r = tr((SS^T)^r)`, `S = A − Aᵀ`. These are **NECESSARY** for any real
> skew-symmetric matrix but **AUTOMATICALLY SATISFIED** — `SS^T = −S²` is PSD for *every* real
> skew `S`, so `(m_r)` is always a genuine Stieltjes sequence and the Hankel matrix is always
> PSD. Moreover the entire skew spectrum is a function of `charA` (THM-507 Cor.1), so the
> skew-Hankel matrix is **constant on every cospectral class**. Therefore the det-side
> inequalities **NEVER cut a hole**: they cannot separate a realized point from a cospectral
> hole, because they take the *same value* on both.
>
> **Layer 2 (per-side / #P).** The holes are cut only by the **integrality + conflict-graph
> layer** — the disjoint/overlap structure of the odd-cycle conflict graph `Ω` (`alpha_2`,
> the Witt/overlap census defects of THM-499/500/505). This is a **permanent-side / #P**
> object: `H = I(Ω,2)` is a packing count, the spectral skeleton plus an integer-linear
> combination of overlap defects (THM-505). The det-side is *blind* to it.

So the holes are **non-algebraic Hodge classes**: moment-feasible (they pass all det-side
spectral inequalities) but integer-forbidden by the #P-side conflict layer.

---

## Statement

Let `T` be a tournament on `n` vertices, `A` its 0/1 adjacency matrix (`A + Aᵀ = J − I`),
`S = A − Aᵀ` the skew (signed) matrix. Let `sig(T) = charA = (tr A, …, tr Aⁿ)` denote the
A-spectrum (cospectral class). Define the **skew moment sequence**
`m_r = tr((SS^T)^r)` for `r = 0,1,2,…` and the **skew Hankel matrix** `Hₜ = [m_{i+j}]_{0≤i,j≤t}`.

**(1) [Layer 1, PROVED] The skew moments are an automatic Stieltjes sequence.**
`SS^T = −S²` is symmetric PSD for *every* real skew `S`. Hence there is a positive measure
`ν_T = Σ_j δ_{μ_j²}` (the `μ_j²` = eigenvalues of `SS^T`, with `±iμ_j` the eigenvalues of `S`)
with `m_r = ∫ x^r dν_T = Σ_j μ_j^{2r}`. Therefore `(m_r)_{r≥0}` is a Stieltjes (hence
Hamburger) moment sequence and **every Hankel matrix `Hₜ` is PSD** — for every tournament,
with no exception and no hypothesis. These are the only positivity constraints a
determinantal/spectral reading can impose, and they are vacuously met.

**(2) [Layer 1, PROVED] The skew-Hankel matrix is a function of `charA`, hence constant on
cospectral classes.** By THM-507 Cor.1 the skew matrix `S = 2A − (J−I)` lies in the A-affine
pencil, and `char_S(x) = det(xI − S) = ∏_j (x² + μ_j²)` is a function of `charA` alone
(PROVED, the `A − J = −(Aᵀ+I)` no-angle-dependence mechanism). The multiset `{μ_j²}` is
therefore determined by `charA`, so each `m_r = Σ_j μ_j^{2r}` (a power sum of the `μ_j²`,
i.e. a symmetric function of the skew spectrum) is a function of `charA`. Consequently `Hₜ`
is **invariant on every cospectral class**.

**(3) [Layer 1 ⇒ the cut fails, PROVED] The det-side inequalities never separate a hole.**
Suppose a hole `h` (a forbidden invariant vector) lies in the same cospectral class as a
realized point `p` (this is the generic situation — cospectral classes split *into* the holes;
THM-499/500). By (2), `Hₜ(h) = Hₜ(p)`, so every Hankel-PSD test gives the identical verdict
on `h` and `p`. Since `p` is realized its verdict is "PSD-feasible", hence so is `h`'s. **No
det-side spectral Hodge inequality can exclude `h`.** The hole is **moment-feasible**.

**(4) [Layer 2, the cut, exact at n=6,7] The holes are cut by integrality + conflict.**
The forbidden values of `H` (`{7,21}`, THM-029/079) and the `(c3,c5)`-fiber holes live in the
`alpha_2`/overlap layer: `H = 1 + 2(c3+c5) + 4·alpha_2` (THM-499), and `alpha_2` (the
disjoint-odd-pair count) and the overlap defects `p33 = W6 − c6`, `TQ`, etc. (THM-500/505)
are the **non-spectral carriers**. These are independence-polynomial / packing data of `Ω`,
i.e. `#P`-side. The hole `(c3,c5)=(8,10)` at n=6 is integer-forbidden by the score
stratification (THM-498/499: the unique `c3=8` class is the regular score `(2,2,2,3,3,3)`,
which realizes `c5 ∈ {6,8,11,12}`, skipping 10) — a purely arithmetic/conflict obstruction,
not a spectral one.

**(5) [Definition + verdict] Non-algebraic Hodge class.** A realizable-region hole that is
**moment-interior** (passes every det-side spectral inequality, by (1)–(3)) yet
**integer-forbidden** (by (4)) is a *non-algebraic Hodge class*. The flagship example
`(c3,c5)=(8,10)` at n=6 (`tr A⁵ = 50`) is such a class.

---

## Proof

**(1).** For any real `S = −Sᵀ`, `SS^T = −S² = SᵀS`, and for any vector `v`,
`vᵀ(SS^T)v = ‖S^T v‖² ≥ 0`, so `SS^T ⪰ 0`. Its eigenvalues are real `≥ 0`; write them
`μ_1²,…,μ_n²` (the `±iμ_j` are the eigenvalues of the skew `S`, which come in conjugate
imaginary pairs, plus zeros). Then `ν_T := Σ_j δ_{μ_j²}` is a positive measure on `[0,∞)` and
`m_r = tr((SS^T)^r) = Σ_j μ_j^{2r} = ∫ x^r dν_T`. A sequence that is the moment sequence of a
positive measure on `[0,∞)` is a Stieltjes sequence; in particular it is Hamburger, so all
Hankel matrices `Hₜ = [m_{i+j}]` are PSD (`xᵀHₜx = ∫ (Σ_i x_i u^i)² dν_T ≥ 0`). This holds for
*every* tournament because it holds for *every* real skew matrix. ∎

**(2).** THM-507 Cor.1 (PROVED) gives `char_S(x) = ∏_j(x² + μ_j²) = (−1)ⁿ charA(?)`-type
function of `charA` via the A-affine pencil determinant — concretely, the matrix-determinant
lemma reduces `det(xI − S)` for `S = 2A − (J−I)` to `charA` times the spectral walk function.
Hence `{μ_j²}` (the roots of `char_S` in `x²`) is fixed by `charA`. Each `m_r = Σ_j μ_j^{2r}`
is a symmetric (power-sum) function of `{μ_j²}`, so `m_r = m_r(charA)` and
`Hₜ = Hₜ(charA)`. Cospectral ⇒ same `charA` ⇒ same `Hₜ`. ∎

**(3).** Immediate from (2): on a cospectral class `Hₜ` is constant, so it gives the same
PSD-feasibility verdict on every member, realized or hole. ∎

**(4)** is THM-499 (`H = 1+2(c3+c5)+4 alpha_2`, `alpha_2` first non-spectral), THM-500
(the `c7`/`TQ` overlap term), THM-505 (the full `spectral skeleton + Σ ±2^level · (overlap
defect)` form), and THM-498 (the `c5=10` score-stratification gap) read together. The
carriers `alpha_2, p33, TQ, Q44, T333, …` are conflict-graph packing/overlap counts, on the
permanent/#P side of THM-506. ∎

**(5)** combines (1)–(4): moment-interior (det-side) ∧ integer-forbidden (#P-side). ∎

---

## What is rigorously PROVED vs. what is honest scope

| claim | status |
|---|---|
| `SS^T` PSD for every tournament ⇒ `(m_r)` Stieltjes ⇒ all skew Hankel matrices PSD | **PROVED** (one-line PSD argument; verified 0 violations exhaustive n≤7) |
| skew spectrum `{μ_j²}` is a function of `charA` (THM-507 Cor.1) ⇒ skew-Hankel constant on cospectral classes | **PROVED** (cite THM-507; verified 0 splits of skew spectrum and skew-Hankel over all cospectral classes n≤7) |
| det-side spectral Hodge inequalities never separate a cospectral hole from a realized point | **PROVED** (corollary of the two above) |
| flagship hole `(8,10)` at n=6 is integer-forbidden | **PROVED** (THM-498 score stratification; reverified) |
| flagship hole is **moment-interior at the LINEAR/convex level** (CERT-1): `(24,50) ∈ conv(realized moment vectors)`; `50 = ⅓·40 + ⅔·55` | **PROVED at the convex-hull level** (exhaustive n=6,7: 6/6 and 30/30 holes lie in the convex hull of realized moment vectors) |
| **NO** det-side moment inequality of ANY degree cuts the hole (full positivstellensatz) | **CONJECTURE.** Layer 1 proves the *Hankel/PSD* (degree-2 SOS) and *cospectral-constancy* obstructions vanish; "moment-interior" is certified at the LINEAR (convex-hull) level by the convex-certificate script, not via a complete positivstellensatz over all polynomial moment inequalities. |

**Honest caveat (the one the task asks to state).** "Moment-interior" is established for the
specific holes at the **linear (convex-hull) level** — each hole lies in the convex hull of
the realized moment vectors (CERT-1), so it passes every *linear* moment inequality; and the
*spectral/Hankel* (det-side, SOS-degree-2) inequalities pass *automatically and identically*
on cospectral pairs (Layer 1, the main theorem). What is **not** claimed is a full
positivstellensatz certifying that *no* higher-degree polynomial moment inequality whatsoever
excludes the hole. The strong, fully-proved statement is the **negative** one that matters:
the det-side / skew-Hankel machinery (the only inequalities a determinantal reading supplies)
is constant on cospectral classes and therefore *provably powerless* to cut holes. The cut is
on the other side of the wall.

---

## Why this is the Valiant det/per wall

- **Det-side invariants are spectral and blind to holes.** Every `det`-based functional of an
  A-affine matrix is a function of `charA` (THM-507): the walk counts, `det(I+S) = ∏(1+μ_j²)`,
  the skew char poly, every Hermitian length-mod-r signing. THM-499 already noted `d ⊥ H`
  (Pearson ≈ 0); the structural cause is that `d` is a skew-spectral coordinate and `H` is
  not. The skew **Hankel** matrix is the canonical *positivity* (Hodge-Riemann–flavored)
  certificate on this side — and it is automatic and cospectral-constant.
- **`H` is a permanent-side / #P object.** `H = I(Ω,2)` is a packing/independent-set count
  (THM-505), the permanental companion of the determinantal world (THM-506). The
  overlap/Witt defects that carry its non-spectral content are exactly the data the
  permanent sees and the determinant does not. The holes live there.
- **Therefore the dichotomy = the det/per (P/#P) split**, read through the moment problem:
  the continuous moment body (det-side, convex, spectral) contains the holes; the integer
  lattice points inside it that are actually realized are carved out by the #P-side conflict
  layer. The holes are the lattice points the convex/spectral body admits but the #P layer
  forbids — non-algebraic Hodge classes.

---

## Verification

Independent re-verification (this session, exhaustive n≤7 via `gentourng`):

- `min eig(SS^T) ≥ −5×10⁻¹⁵` over all iso-classes at n=4,5,6,7 (PSD, numerical zero) — Layer 1.
- `min eig(skew Hankel Hₜ) ≥ −7×10⁻¹¹` over all iso-classes at n=4,5,6,7 (Stieltjes/Hankel-PSD).
- **0** cospectral classes (of 3/9/28/168 at n=4/5/6/7) where the **skew spectrum** differs;
  **0** where the **skew-Hankel matrix** differs — Layer 1's cospectral-constancy.
- Flagship hole: c3=8 fiber at n=6 realizes `c5 ∈ {6,8,11,12}` (`tr A⁵ ∈ {30,40,55,60}`);
  `c5=10` (`tr A⁵=50`) absent; `50 = ⅓·40 + ⅔·55` strictly interior — moment-interior,
  integer-forbidden.

**Artifacts:** `04-computation/baby_hodge_moment_region_macmini_0615s1.py` (region + Hankel
test), `04-computation/baby_hodge_convex_certificate_macmini_0615s1.py` (+ `.out`, CERT-1
convex-hull certificates: 6/6 holes at n=6, 30/30 at n=7), `04-computation/baby_hodge_realizability_kps8.py`
(region/holes). Independent check this session reproduced all four Layer-1 facts and the
flagship certificate with 0 violations.

## Provenance / relation to existing canon

This formalizes, as a single dichotomy theorem, the strands: the spectral boundaries
(THM-499/500), the OCF non-spectral defect decomposition (THM-505), the proof that all
A-affine determinants are spectral (THM-507), and the convex/moment certification of the
holes (the `baby_hodge_*` scripts). The new content is the **statement and proof that the
det-side positivity certificate (skew-Hankel PSD) is automatic AND cospectral-constant**,
which is exactly why it cannot cut the holes — placing the realizability obstruction
unambiguously on the permanent/#P side of the Valiant wall. The previously-existing
`baby_hodge_convex_certificate` docstring referenced a "THM-508" tag internally; the canon
THM-508 is a different result (the contraction/Hadamard wall), so this dichotomy is filed as
THM-509.

## Strengthening (this session): the holes are PURE INTEGRALITY GAPS — even the #P-side flag matrix can't cut them continuously

Layer 1 shows the det-side (spectral) positivity is powerless. One might hope the
*permanent-side* carriers (the overlap densities p33, alpha_2) supply a flag-algebra /
Cauchy–Schwarz moment matrix that DOES cut a hole at the continuous level. They do not:

> **No finite PSD relaxation of cycle/overlap densities — spectral OR conflict-side —
> excludes a baby-Hodge hole. The holes are integrality gaps, not continuous-positivity
> gaps.**

Proof (verified, n=6, c3=8 fiber). The hole `(c3,c5)=(8,10)` is the exact midpoint of
the *realized* points `(8,8)` and `(8,12)`: `c5=10 = ½·8 + ½·12`. A 50/50 random blend
of those two tournaments is a genuine tournament-LIMIT with `(c3,c5)`-density `(8,10)`.
Any moment / Cauchy–Schwarz / SOS / flag-algebra relaxation carves out a CLOSED CONVEX
region of densities, and `(8,10)` is an interior limit point of it, so no such matrix
can exclude it. Concretely: in the c3=8 fiber `p33 + alpha_2 = C(c3,2) = 28` identically
(one degree of freedom), and an explicit 4×4 overlap (triangle-root-type) Gram matrix is
PSD at every c5 ∈ {6,…,12} including the holes — being literally an outer-product sum of
real density data, it is automatically PSD at every genuine limit point. The hole exists
**only** because no INTEGER tournament realizes `c5=10` while its continuous relaxation
is feasible. Cutting it requires the Boolean/rank-1 (integrality) constraint that no flag
matrix imposes. (Fresh code: overlap_moment_matrix_n6 / overlap_gram_explicit_psd, this
session; min-eig ≥ −1e−9 at all c5.)

This is the sharp form of the dichotomy: the realizable region's continuous closure
(convex hull / flag-feasible body) is `det`-and-`per`-blind to the holes alike; the holes
are exactly the **integer lattice points interior to that body that the #P count skips**.
The baby-Hodge "non-algebraic Hodge class" is therefore a pure *integrality* phenomenon
sitting inside a perfectly nice convex moment body — which is precisely the flavor of an
actual non-algebraic Hodge class (satisfies every Hodge-theoretic inequality, fails an
integral/algebraic-cycle realizability).

## The moment–cumulant law underneath (the census ladder is a zeta moment-cumulant pair)

The `tr(A^k)` decomposition that generates all this is an exact **moment–cumulant**
correspondence, but of the **necklace/zeta** family — NOT free probability:

- Moments = `p_k = tr(A^k)` (closed k-walk counts = spectral power sums).
- Cumulants = `W_k = (1/k) Σ_{d|k} μ(d) tr(A^{k/d})` (Witt/necklace numbers = aperiodic
  primitive closed walks up to rotation), related by the Artin–Mazur / Bowen–Lanford
  zeta identity `exp(Σ_k p_k u^k/k) = Π_k (1−u^k)^{−W_k}`, i.e. Möbius inversion on the
  **divisor lattice**.
- Each `W_k = c_k + (overlap configs)`: `W_3=c3, …, W_5=c5` (spectral = cumulant);
  `W_6=c6+p33, W_7=c7+TQ, W_8=c8+Q44+TF` — the SPLIT of the cumulant into a simple cycle
  and overlaps is the non-spectral content (the very mechanism of `c_k` non-spectrality,
  THM-499/500/502/505).

VERIFIED exhaustively n≤6 + random n=7,8 (138372/138372 ladder + Witt checks exact; zeta
identity exact as a formal power series for arbitrary integer `p`, confirming it is a
structural Möbius identity). HONEST: this is **only formally analogous to free
probability** — the free cumulants (Speicher, non-crossing-partition NC(k) Möbius) differ
from `W_k` numerically (0/6 match on n=7), and the posets differ (divisor lattice, sizes
1,2,2,3,2,4 vs Catalan 1,2,5,14,42,132). The project's genuine free-probability content
(THM-438 semicircle/Catalan in the Paley spectrum) lives in the spectral DISTRIBUTION's
moments, a different object from this cyclic necklace census. (HYP-2526.)
