# Audit J1: THM-4910, THM-4911, HYP-9500 (seven-point sphere, signature simplex, equicentric tournaments)

2026-11-07. Independent audit of the mac-mini-2026-11-07-sevensphere packet at commit `14d7c11577`
(branch `mac-mini/session-20261107-sevensphere`). No packet file was edited.

Audit scripts and outputs: `04-computation/experiments/sevensphere_anticomplementary_20261107/audit_J1/`
(every script carries a one-line header; total runtime about 5 minutes; peak memory 0.47 GB).

**Overall verdict: PASS WITH CORRECTIONS.**
- The core mathematics is correct and was re-derived independently: THM-4910 statements 1-4 (symbolically for all
  `n`, with `S2` defined only by the corner powers, and exactly in coordinates for `n = 3..7`), and THM-4911
  statements 1-6 (proofs checked; census of all isomorphism classes `n <= 10` reproduces every number in the packet,
  with zero violations).
- Three claims are false as stated (errors):
  - THM-4910(6), "no flag extension", is stated universally. It holds only generically: the flag system is consistent,
    with a real sphere, on a real hypersurface of nondegenerate 4-simplices.
  - THM-4911's title says the eight-point sphere "transfers to every non-equicentric tournament". Two classes at
    `n = 5` and six at `n = 7` are non-equicentric and have no vertex in the regime.
  - HYP-9500(a) is typed OPEN at the orders `n = 0 (mod 4)` without known skew-Hadamard matrices. It is provable in a few
    lines at every even `n >= 6`: substituting an equicentric tournament into a regular one preserves equicentricity.
- Several overclaims, one non sequitur (the GF(2) reading of "complementary to V1V2V3") and missing credits
  (THM-1475(1), the June kernel-parity lemma, THM-472, Nozaki-Suda 2012).

## Verdict table

| item | verdict |
|---|---|
| THM-4910 (1) corner powers, existence, uniqueness | PASS (cosmetic: uniqueness citation) |
| THM-4910 (2) eight points, Vieta on the sides | PASS |
| THM-4910 (3) seven-point sphere with `C*` | PASS WITH CORRECTIONS (distinctness is generic only) |
| THM-4910 (4) defect `V0A·V0C* = Σ` | PASS |
| THM-4910 (5) mod 2 | PASS WITH CORRECTIONS (congruence correct; GF(2) "names `C*` correctly" is a non sequitur) |
| THM-4910 (6) no flag extension | FAIL as a universal statement; PASS once restated for generic simplices |
| THM-4910 remarks, status, credits | PASS WITH CORRECTIONS |
| THM-4911 (1) cube, Hamming distance `1 + δ` | PASS |
| THM-4911 (2) `s = 1/dist²`, `r² = n − 1/s`, bounds and equality cases | PASS (typing: producer check is floating point) |
| THM-4911 (3a)-(3d) equicentric classification, `λ ≡ n−1 (mod 4)` | PASS (credits; (3c) omits `n = 2`) |
| THM-4911 (4) `c ∈ {−1, 0}` | PASS (credit THM-472) |
| THM-4911 (5) QR7 = Fano line code | PASS |
| THM-4911 (6) transfer of the eight-point sphere | PASS (statement); title FAIL ("every non-equicentric tournament") |
| THM-4911 remarks | PASS WITH CORRECTIONS (overclaims) |
| HYP-9500 | PASS WITH CORRECTIONS (proved parts and table correct; (a) is in fact PROVED at every even `n >= 6`) |
| Note §0 | PASS WITH CORRECTIONS (credits) |
| Note §1 | PASS WITH CORRECTIONS (flag statement; unscripted scan claim) |
| Note §2 | PASS WITH CORRECTIONS (GF(2) non sequitur; "precisely ... regularity") |
| Note §3 | PASS WITH CORRECTIONS (HYP-9500 status; "By Reid-Brown") |
| Note §4 ledger | PASS WITH CORRECTIONS ("new invariant", "cube", coth caveat, ambiguous inclusion chain) |
| Note §5 ledger | PASS WITH CORRECTIONS |
| INDEX line | PASS WITH CORRECTIONS |

## A. THM-4910

**Re-derivation (j1_geometry_symbolic.py).** In THM-4742's Gram basis `(v, s, w)`, with `n, R², v·s, v·w, |s|², k`
symbolic, I defined `S2` *only* by the four corner powers `c_P (R² − |OP|²)`, `c = n+1, n, n−1, n`. Then:
- the unique solution is THM-4742's `S2` (all four coefficients agree);
- `Pow_S2` vanishes identically at all eight labels `A, B, C*, D, E, F, H', X`;
- `Pow_S2(V0) = Σ = 2nR² − 2v·s` and equals the signed product `V0A·V0C*`;
- the top quadruple `{D, F, C*, H'}` has determinant `n(n−1)` in basis coordinates;
- the bottom quadruple `{E, B, A, X}` has a determinant that is not identically zero. It is generically, not provably
  always, non-coplanar.

By hand: the corner powers follow from THM-4742's coefficients in two lines each. `Q(C*) = R²[(n−1)² + (n−1)(2−n) − n + 1] = 0`,
and `Q(X) = n r'² − nρ γ² = 0` with `ρ = r'²/γ²`.

**Exact coordinates (j1_geometry_exact.py).** I used 200, 60, 30, 15 and 8 random rational simplices for
`n = 3, 4, 5, 6, 7` in `Q^n`, in the regime, with `S2` again from the corner powers. In every case:
- all labels lie in the 3-flat and on `S2`, and the radius is real;
- on each side line the power restricted to the line has exactly the claimed roots (Vieta);
- both roots lie strictly beyond the far corner of their side, and all four corners have positive power;
- the defect equals `Σ` and `V0A·V0C*`, with `Σ > 0`.

**Proof check.** The Vieta argument is correct, including signs. All four corners lie inside the circumsphere and have
positive power with respect to `S2`. So on each side both points of `S2` lie on the same side of each corner, beyond the
far corner:
- `(n+1)g, R²/g > g` on `OG`;
- `(n+1)c, R²/c > c` on `OO'`;
- `nγ, r'²/γ > γ` on `O'G'`;
- `nm, a > m` on the median.

`|OG'|² = c² + γ²` holds because `OO'` is normal to the facet. Statement 4's parallel-axis computation is correct.
Statement 1's appeal to THM-4742(1) needs `D ≠ E` for uniqueness, but the identification holds anyway by comparing
coefficients (P1). The regime "`O, V0, G, O'` not coplanar" is equivalent to "`O, G, G', O'` not coplanar" (both mean
`v, s, w` independent).

**The regime excludes session23b's equilateral-face counterexample.** Take `V0 = (0,0,−1)`, `V1 = e1`, `V2 = e2`,
`V3 = e3`. Then `G' = O'`, so `s = 3w` and `O, V0, G, O'` are coplanar. Every facet with `G' = O'` is excluded.

**Distinctness is generic, not universal.**
- Labels on different sides never coincide in the regime. Non-adjacent sides of a skew quadrilateral are skew lines,
  adjacent sides meet only at a corner, and no label is a corner.
- Within a side, the two points coincide exactly at a tangency:
  - `D = E` iff `|OG|² = R²/(n+1)`;
  - `F = B` iff `|OO'|² = R²/(n+1)`;
  - `H' = X` iff `n|O'G'|² = r'²`;
  - `C* = A` iff `(n−1)|V0G'|² = R² − |OG'|²`.
- All four coincidences occur inside the regime:
  - `D = E` and `F = B` on session23b's exceptional family. Their example `V0 = (0,1,0)`, `Vi = (√3/2,0,½), (0,±√3/2,½)`
    has `det[v s w] = −√3/4 ≠ 0`.
  - `C* = A` for the rational tetrahedron `V0 = (2, 2/3, 4/3)` over the facet `(0,0,0), (4,0,0), (0,2,0)`, with 7 distinct
    labels, all on `S2`.
  - `H' = X` for an exact algebraic example. No rational example with an axis-parallel facet exists, because `3r'²`
    is never a sum of two rational squares.
- So "pairwise distinct" is correct as a statement about the packet's checks (FINITE-EXACT) but must be qualified as
  generic.

**Statement 5.** `C* = −(n−1)V0 + ΣVi` is correct, and for odd `n` it reduces to `{1, …, n}` mod 2. But the conclusion
does not follow. "So 'the vertex complementary to V1V2V3' names `C*` correctly over GF(2)" (THM-4910 line 94; note
§2 line 101) is a non sequitur:
- In the subset model of PG(3,2), `C* ≡ {1,2,3}`, which is the complement of `V0`'s point `{0}`.
- The subset complementary to `{1,2,3}` is `{0} = V0`, the literal reading.
- The title's wording, "'complementary' to `V0`", is correct.
- The GF(2) reduction gives no support to the anticomplementary reading of the pasted phrase.

The Fano remark "in characteristic 2 … the anticomplementary vertex reduces to the complementary point" holds only for
odd `n`. For the triangle itself (`n = 2`, the Fano picture), `C* ≡ V0 + V1 + V2`, the centroid point.

**Statement 6 (j1_flag_exact.py, j1_flag_locus.py).**
- The exact rank test confirms inconsistency for 10/10 random rational simplices at each of `n = 4, 5, 6`, and
  consistency at `n = 3` (THM-4910(1)).
- The minimum float residual over 500 random simplices is 1e-11 to 1e-12. Those near-zeros are near-degenerate flags,
  where a tiny 2-face has its circumcentre on `S`.
- There are also genuine exceptions. At `n = 4` the system (6 equations, 5 unknowns) is consistent iff `det[X|y] = 0`,
  and this determinant changes sign on straight paths between integer 4-simplices. On those paths:
  - every face stays nondegenerate (minimum face Gram determinant 18.7 or more);
  - `X` keeps full rank (`cond(XᵀX) <= 1.5e4`);
  - bisection to width 2e-14 leaves residual `<= 1e-13 R²`;
  - the flag sphere is real (radius² = 87R², 149R², 227R²).
- So consistent simplices form a real hypersurface. "There is no sphere …" and "it is inconsistent" are true only
  generically. The packet's "worst residual 0.16, 0.19, 0.12" is the wrong statistic for a NEGATIVE claim.
- Status: error in a NEGATIVE claim. The intended statement, that the rule is not an identity, is true.

**Unscripted claim (note §1).** "S49's question 2 …: a numerical scan of ten candidate coincidences found none …" has
no script or output. The October audit removed a claim of exactly this kind (MISTAKE-650).

**Credits.**
- Of the "new" items, the power `(n+1)R²` at `O`, the `ι`-pairs `{D,E}` and `{F,B}` on `OG` and `OO'`, and the circle
  `{D,E,F,B}` in plane `OGO'` are THM-4742's.
- THM-4742 also puts `ι(A)` and `ι(H')` on `S2`, so eight is not a closed count. Every corner inversion generates more
  points.
- New are `C* ∈ S2`, `X`, the corner powers at `G, G', O'`, the pairs on `GG'` and `G'O'`, and the defect as
  `V0A·V0C*`.
- No prior occurrence of "anticomplement*", `4G − 3V0`, `s − 2v`, "inverse of G'", "eight-point sphere" or "skew
  quadrilateral" was found in the repo (rg over all tracked text, including `session*/` and
  `05-knowledge/results/session*_*`).
- A brief web search found no "seven-point sphere" in the literature. The packet's terminology (homothety of ratio `−n`
  about `G`) matches the standard anticomplementary/medial tetrahedron.

**Reproducibility.** The producer's four scripts re-run and reproduce their committed outputs exactly (section F).

## B. THM-4911

**(1).** The proof is correct, and the coordinates `i, j` contribute exactly 1. Identity:
`δ_ij = (n − 2 − (MᵀM)_ij)/2`. Checked combinatorially on all classes `n <= 9` (Hamming distance of the columns
against the separator count from `A`).

**(2).** The projection argument is correct. It needs `|p_i|² = n` for all `i`, so that the circumcentre is the
projection of 0. `s_THM-4911 = s_THM-1950`, because `|(I+K)⁻¹u|² = uᵀ((I+K)(I−K))⁻¹u = uᵀ(I+KᵀK)⁻¹u`.

The bounds and both equality cases were checked:
- exactly (Fractions) on all 530 classes `n <= 7`;
- in floating point on all classes `n = 8, 9`, with an independent projection of 0 onto the affine hull.

The producer's (2) check is floating point, so "(1), (2), (4) FINITE-EXACT for `n <= 8`" over-types (2).

**(3).**
- `λ = |x|²/n` is correct.
- (3a) is correct: `λ·uᵀz = −uᵀMx = (Mz)ᵀx = 0` and `uᵀz` is odd.
  - The identity `adj M = zzᵀ` for odd skew `M` is classical but uncited. The argument only needs `Mz = 0`, or `rank M =
    n−1` plus the all-odd kernel vector.
  - Checked on all classes `n = 3, 5, 7`: `Mz = 0`, `adj M = zzᵀ`, all `z_i` odd, `Σz` odd. Also on 400 classes at
    `n = 9` and 60 random tournaments at `n = 11`.
  - Independently, `u ∉ range(M)` for every class of odd order `<= 9`.
- (3b) is correct, including `R_1` as a regular part.
- (3c) is correct. Two remarks:
  - `|S| = n/2` is automatic for almost-regular `T`.
  - The single arc at `n = 2` is almost regular and equicentric with `λ = 1`, so "exist exactly when `n ≡ 2 (mod 4)`,
    `n ≥ 6`" omits `n = 2`. Strong ones exist for `n >= 6`.
  - The circulant-cross construction was checked for `k = 3..13`: almost regular, equicentric with `λ = 1`, strong.
- (3d) is correct; the mod-4 computation was redone.

**Census (j1_tournament_census.py, vectorised; j1_labeled_crosscheck.py).**
- All 9,733,056 classes at `n = 10` were processed in 28 s.
- The census reproduces exactly:
  - equicentric counts 1, 2, 1, 5, 3, 14, 15, 165;
  - non-strong counts 0, 2, 0, 3, 0, 8, 0, 37;
  - almost-regular counts 0, 0, 0, 2, 0, 0, 0, 88;
  - the `λ`-sets.
- Strong equicentric counts are 2, 6, 128 at `n = 6, 8, 10`. They agree with `gentourng -c`. The strong `λ`-sets are
  `{1}`, `{3,7,11}` and `{1,5,9,13,17}`.
- There are zero violations of (1)-(4) (for (4), all classes `n <= 8` directly and all equicentric classes `n <= 10`).
- The labeled cross-check for `n <= 7`, independent of gentourng, gives labeled equicentric counts equal to
  `Σ n!/|Aut|` over the classes.

**(4).** The proof is correct (`MᵀM` commutes with `M`; `x ⊥ u`). The cases `c = −1` (`SSᵀ = nI − J`, the DRT engine)
and `c = 0` (`SSᵀ = (n−1)I`, skew conference, equivalently `I + S` skew-Hadamard) are THM-472's equality cases and
should be credited. The dichotomy itself (no other `c`) was not found in the literature in a brief search; it is likely
folklore. The first equivalent form of the skew-Hadamard conjecture is the definition, not Reid-Brown.

**(5).** Correct. The design parameters `2-(n, (n−1)/2, (n−3)/4)` are classical.

**(6), j1_transfer_exact.py.**
- 160/160 in-regime pairs `(T, v)` were checked exactly. The check used `S2` from the corner powers, `s` from its
  definition `uᵀ(I+MᵀM)⁻¹u`, and `Δ_v` by brute-force enumeration of directed 2-paths:
  - all eight labels lie on `S2`;
  - `Pow(O_T) = n r_T² = n(n − 1/s)`;
  - `Pow(p_v) = 4(n − 1 + Δ_v) = Σ_j |p_v − p_j|²`.
- The statement is correctly hedged ("such that … in THM-4742's regime at v").
- **The title is false.** j1_regime_failures.out lists non-equicentric classes with no vertex in the regime:
  - at `n = 5`, source ⇒ C3 ⇒ sink (scores 4,2,2,2,0; `G' = O'` at the source and sink, `O' = O` at the C3 vertices);
  - at `n = 5`, a second class with scores 3,2,2,2,1 (`v, s, w` dependent at every vertex);
  - six classes at `n = 7`.

**Remarks.**
- "The Fano tournament is the one where all of the sphere's incidences collapse" overclaims. Every equicentric
  tournament, and some non-equicentric ones, lack the sphere.
- The "hierarchy" sentence is true only for odd `n`. For even `n`, equicentric is a large class (165 at `n = 10`).

**Suggested upgrade (verified, j1_misc_checks.out).** `(M²u)_v = 2Δ_v − (n−1)² = 4SB(v) − 2(n−2)d⁺(v) − (n−1)`, where
`SB(v) = Σ_{v→k} d⁺(k)` is the Sonneborn-Berger (second out-) score. So:
- `T` is equicentric iff every vertex is an endpoint of the same number `Δ` of directed 2-paths, iff `2SB(v) − (n−2)d⁺(v)`
  is constant;
- then `λ = (n−1)² − 2Δ`.

This is the classical simplex fact (centroid = circumcentre iff every vertex has the same sum of squared edge lengths)
read through (1). It also shows that `Δ_v` is not a new invariant (note §4):
`Δ_v = 2SB(v) + (n−1)d⁻(v) + d⁺(v) − C(n,2)`.

## C. HYP-9500

**The proved parts are correct (j1_hyp9500.out).**
- Bordering a doubly regular `R` with a source gives a skew-Hadamard `I + M` (block computation:
  `I − M² = diag((m+1)I, m+1)`).
- Switching at `w` is `M → DMD` and preserves `(I+M)ᵀ(I+M) = nI`.
- The strongness proof is correct. A vertex `y` with `y → w` shares `(m−3)/4 >= 1` common out-neighbours `z` with
  `w`, and `y → z → w` survives the switch.
- Checked for every `w`: Paley `q = 7, 11, 19, 23, 31, 43, 47`, and the non-Paley order-15 core of the skew-Walsh
  matrix of order 16.
- The evidence table and the strong `λ`-sets are reproduced exactly. The strong `n = 8`, `λ = 3` class does have
  `x = (−3, −1, −1, −1, 1, 1, 1, 3)`.

**(a) is PROVED at every even `n >= 6`, so the OPEN typing is wrong.**
- *Substitution lemma.* Let `R` be regular of odd order `r >= 3`, and `T2` equicentric of order `m` with
  `−M2²u = λ2 u`. Then `R[T2]` (each vertex of `R` replaced by a copy of `T2`) is strong and equicentric with
  `λ = λ2`.
  - Proof: `M = M_R ⊗ J + I ⊗ M2` and `x_R = 0` give `Mu = u ⊗ x2`.
  - So `M²u = x_R ⊗ (uᵀx2)u + u ⊗ M2x2 = −λ2 u`.
  - `R[T2]` is strong because `R` is.
- *Construction at every even order.*
  - If `n = 2^e·r` with `r` odd, `r >= 3`, `e >= 1`, take `R_r[R_1 ⇒ R_(2^e − 1)]`, with `λ = 2^e − 1`. For `e = 1` this
    is `R_r[single arc]`, with `λ = 1`.
  - If `n = 2^e` with `e >= 3`, THM-447(2)'s skew-Walsh tower gives a skew-Hadamard matrix of order `2^e`. HYP-9500's own
    switching construction then applies.
  - Verified for every even `n` from 6 to 96.
- The orders "where skew-Hadamard matrices are not known" all have an odd factor (276 = 4·69, …). So (a) carries no
  information about the skew-Hadamard conjecture, and the "Why this matters" paragraph overclaims.

**(b) stays OPEN (CONJECTURE).** The data `n <= 10` support it. Two composition rules help: `λ(R[T2]) = λ(T2)` and
`λ(T1[R]) = r²λ(T1)` for regular `R`.

## D. Credits and priority

**In the repo** (rg over all tracked text, including `session*/` and `05-knowledge/results/session*_*`). Not found
anywhere outside the packet:
- "equicentric", "signature simplex", "anticomplement*";
- the `s = 1/dist²` reading (no hit linking "inverse response" to distance, affine hull or projection);
- an `M²𝟙 ∥ 𝟙` condition.

Prior in-repo ingredients that must be credited:
- **THM-1475(1)** (klein-2026-07-20-S339): `Pf(S(T))` is odd for every even-order tournament, because
  `S ≡ J − I (mod 2)` and `(n−1)!!` is odd. This is exactly (3a)'s "sub-Pfaffians are odd". THM-3729 (rooted Pfaffians,
  line 141) uses the same count, but THM-1475 is the direct statement.
- **07-reflections/border-remembers-doubling-forgets-the-gleason-ring-from-two-paley-tournaments-cbx0611s3.md**
  (2026-06-11): "the kernel vector w of the skew matrix of ANY odd-order tournament has ALL-ODD integer entries". This is
  the core of (3a).
- **THM-472**: odd-`n` equality `SSᵀ = nI − J` (DRT engine) and even-`n` equality `SSᵀ = (n−1)I` (skew conference,
  equivalently `I+S` skew-Hadamard). These are (4)'s two cases.
- THM-447, THM-1950, THM-211, THM-4742 and session23b are credited correctly. Nozaki-Suda is already in-repo
  (`03-artifacts/drafts/web-synthesis-opus-S5.md`, TANGENTS, INVESTIGATION-BACKLOG).

**Literature** (brief web search).
- Reid-Brown, *Doubly regular tournaments are equivalent to skew Hadamard matrices*, J. Combin. Theory Ser. A 12
  (1972) 332-338: credited correctly.
- **Nozaki-Suda**, *A characterization of skew Hadamard matrices and doubly regular tournaments*, arXiv:1202.5374
  (2012). They work with the main angles of the Hermitian skew Seidel matrix `√−1(A − Aᵀ)` of a tournament, so
  equicentric is exactly "at most two main eigenvalues `±√λ`". In their paper:
  - Lemma 2.1(3) is regular iff `S𝟙 = 0`;
  - Lemma 2.2: for even `n`, almost regular iff `Σθ²β² = 1`, i.e. `|x|²/n = 1`;
  - Theorem 2.7 concerns order `4k+2`, main eigenvalues `±1` (equicentric with `λ = 1`) plus two non-main eigenvalues
    `±√(4k+3)`. Such tournaments are equivalent to DRT(`4k+3`).
  - This is the closest prior work to (3c) and to HYP-9500's "skew-Hadamard shadow" motivation, and must be cited.
- No source was found classifying tournaments with at most two skew main eigenvalues ((3a)-(3d)), nor the
  `c ∈ {−1, 0}` dichotomy. Both are "new as far as a brief search shows". The dichotomy is likely folklore.

## E. Typing and overclaims

- **THM-4910 status.** "PROVED: statements 1-5" covers statement 5's interpretive sentence, which is STRUCTURAL and,
  as written, a non sequitur. "NEGATIVE (FINITE-NUMERICAL)" for (6) should be "NEGATIVE, generic (exact rank test);
  exceptions on a real hypersurface".
- **Fano paragraphs** (THM-4910 remarks; note §2). Correctly typed ANALOGY, except:
  - the char-2 reduction holds only for odd `n`;
  - "The Fano tournament is therefore exactly where the seven-point sphere does not exist. More generally the sphere
    sees precisely a tournament's deviation from regularity" is false. The sphere fails for every equicentric
    tournament, including non-regular even-order ones (every `R_a ⇒ R_b`, every skew-Hadamard tournament), and for some
    non-equicentric ones.
- **Note §4 ledger.**
  - Row 1: "the obstruction was a misreading" presumes the intended reading.
  - Row 2: `h_(G,−n)` fixes `G`, so the "cube of centres" is two skew quadrilaterals sharing the corner `G`. The power
    `Pow_(S2^W)(D) = n²(n+1)R²` is correct.
  - Row 3: the `d²` composition law is a one-line restatement of THM-1950 (KNOWN). It was checked on 40 exact pairs.
    The "coth-addition" form has real rapidity only when `s < 1`; 18 of the 40 pairs had some `s > 1`. The same caveat
    applies to THM-1950's "tanh of an additive rapidity".
  - Row 4: "Equicentric ⊂ equilateral ⊂ skew-Hadamard (as conditions)" is coherent only if `⊂` means "is a sub-condition
    of". Read as classes it is reversed, and `equilateral ⊂ skew-Hadamard` is false (DRTs). Write implications.
  - Row 6: "a new tournament invariant" overclaims (see the upgrade in B).
- **Note §3.** "By Reid-Brown, the skew-Hadamard conjecture … pairwise orthogonal signature vectors" is the definition.
  The HYP-9500 status there is wrong (C).

## Corrections (old → new; severity)

THM-4910 (`01-canon/theorems/THM-4910-…md`)
1. **error.** Line 11.
   - Old: `NEGATIVE (FINITE-NUMERICAL): statement 6 (no flag extension) at n = 4, 5, 6.`
   - New: `NEGATIVE, generic (exact rank test on random rational simplices, audit J1): statement 6 at n = 4, 5, 6; the
     flag system is consistent, with a real sphere, on a real hypersurface of nondegenerate 4-simplices (audit J1,
     j1_flag_locus.py).`
2. **error.** Lines 98-99.
   - Old: `- There is no sphere with … for all k.` / `- The worst affine-fit residuals are 0.16, 0.19, 0.12 at n = 4, 5, 6.`
   - New: `- For a generic simplex there is no sphere with … for all k: the 2(n−1) conditions on n+1 unknowns are
     inconsistent (exact rank test, n = 4, 5, 6). - The exceptions form a proper algebraic subset; for n = 4 it is a real
     hypersurface of nondegenerate simplices carrying a real flag sphere (audit J1).`
   - Also line 146: `it is inconsistent.` → `it is inconsistent for generic simplices.`
3. **overclaim.** Line 81.
   - Old: `In the exact checks all eight labels are pairwise distinct.`
   - New: append `This is generic: labels on different sides never coincide in the regime, but each side pair can merge
     into a tangency point inside the regime (D = E and F = B on session23b's exceptional family; C* = A for
     V0 = (2, 2/3, 4/3) over the facet (0,0,0), (4,0,0), (0,2,0); H' = X for an algebraic example; audit J1).`
4. **overclaim (non sequitur).** Lines 94-95.
   - Old: `So "the vertex complementary to V1V2V3" names C* correctly over GF(2). Over R it must be read as the
     anticomplementary vertex.`
   - New: `So over GF(2), C* is the face's own point {1,2,3}, the complement of V0's point {0}. In the same subset model
     the point complementary to {1,2,3} is {0} = V0, so the GF(2) reduction does not select C* as 'the vertex
     complementary to V1V2V3'; the anticomplementary reading remains an interpretation.`
5. **overclaim (ANALOGY scope).** Line 162.
   - Old: `anticomplementary vertex reduces to the complementary point (statement 5).`
   - New: `anticomplementary vertex reduces to the complementary point for odd n (statement 5); for the triangle itself
     (n = 2) it reduces to the centroid point V0 + V1 + V2.`
6. **cosmetic.** Line 113.
   - Old: `By THM-4742(1), S2 is the unique such sphere, so the two coincide.`
   - New: `By THM-4742(1) (when D ≠ E), or in general by comparing coefficients in the Gram basis (P1), the two coincide.`
7. **cosmetic.** Lines 168-169.
   - Old: `The 'top' and 'bottom' quadruples … are not coplanar. In THM-4742's basis the determinant is n(n−1) ≠ 0`
   - New: `The 'top' quadruple {D,F,C*,H'} is never coplanar (determinant n(n−1) ≠ 0 in THM-4742's basis), so there is no
     Miquel cube; the 'bottom' quadruple {E,B,A,X} is generically non-coplanar.`
8. **cosmetic.** Line 152.
   - Old: `valid over any field in THM-4742's regime`
   - New: `valid over any field in THM-4742's regime (with |G' − O'|² ≠ 0 for X)`
9. **credit.** Line 21.
   - Old: `the corner-power description and the side-pair structure;`
   - New: `the corner powers at G, G', O' and the side pairs on GG' and G'O' (the power (n+1)R² at O, the ι-pairs {D,E},
     {F,B} and the circle {D,E,F,B} are THM-4742's);`

THM-4911 (`01-canon/theorems/THM-4911-…md`)

10. **error.** Title.
    - Old: `THM-4910's eight-point sphere transfers to every non-equicentric tournament, with`
    - New: `THM-4910's eight-point sphere transfers to every tournament and vertex v in THM-4742's regime (never for an
      equicentric T, and not for some non-equicentric T such as source => C3 => sink), with`
11. **typing.** Line 7.
    - Old: `(1), (2), (4) for n <= 8;`
    - New: `(1), (4) for n <= 8; (2) in floating point for n <= 8 (audit J1: exact for n <= 7, floating point for n <= 9);`
12. **credit.** Line 19.
    - Old: `sub-Pfaffians of +-1 skew matrices are odd (the rooted-Pfaffian technique of THM-3729).`
    - New: `even-order tournament Pfaffians are odd (THM-1475(1), klein-2026-07-20-S339; the same count in THM-3729), and
      the kernel vector of an odd-order tournament's skew matrix has all-odd entries (07-reflections …-cbx0611s3.md,
      2026-06-11);`
    - Add: `the cases c = −1 (SSᵀ = nI − J) and c = 0 (SSᵀ = (n−1)I, skew conference ⟺ skew-Hadamard) are THM-472's
      equality cases; the main-angle framework is Nozaki–Suda (arXiv:1202.5374, 2012): equicentric = at most two main
      eigenvalues ±√λ of √−1·M; their Theorem 2.7 (order 4k+2, main eigenvalues ±1, i.e. λ = 1) is the closest prior
      result to (3c).`
13. **cosmetic.** Line 70.
    - Old: `Such tournaments exist exactly when n ≡ 2 (mod 4), n ≥ 6, and they have λ = 1.`
    - New: `Such tournaments exist exactly when n ≡ 2 (mod 4) (n = 2: the single arc; strong ones for every n ≥ 6), and
      they have λ = 1.`
14. **cosmetic/credit.** Line 85.
    - Old: `Consequently (KNOWN in substance, Reid–Brown)`
    - New: `Consequently (the first form is the definition of a skew-Hadamard matrix; the second is Reid–Brown)`
15. **overclaim.** Lines 194-195.
    - Old: `the Fano tournament is the one where all of the sphere's incidences collapse into symmetry.`
    - New: `the Fano tournament is one of the many (every equicentric tournament, e.g. every regular one) for which the
      sphere's construction collapses (G_T = O_T).`
16. **overclaim.** Line 188.
    - Old: `So the classical tournament hierarchy (regular, doubly regular, skew-Hadamard) is the hierarchy of simplex
      symmetry (equicentric, equilateral).`
    - New: `For odd n, regular = equicentric; for all n, equilateral = doubly regular or skew-Hadamard; for even n,
      equicentric is a much larger class with no classical name.`
17. **upgrade (optional).** Add to (3): `T is equicentric iff every vertex is an endpoint of the same number Δ of directed
    2-paths, iff 2SB(v) − (n−2)d⁺(v) is constant (SB = Sonneborn–Berger score); then λ = (n−1)² − 2Δ.`

HYP-9500 (`05-knowledge/hypotheses/HYP-9500-…md`)

18. **error (mis-typed status).** Line 12.
    - Old: `OPEN: (a) at the orders n = 0 (mod 4) where skew-Hadamard matrices are not known; (b) in general.`
    - New: `PROVED (audit J1): (a) at every even n >= 6. R_r[T2] is strong and equicentric with λ(T2) for any regular R_r
      (r odd >= 3) and equicentric T2; take T2 = R_1 => R_(2^e − 1) when n = 2^e·r, and at n = 2^e use THM-447's
      skew-Walsh tower and the switching construction. OPEN: (b).`
    - Update line 5 (`OPEN (CONJECTURE), with large proved parts.`) and the title sentences `This is proved for n = 2
      (mod 4) and … should not need the skew-Hadamard conjecture` accordingly.
19. **overclaim.** Lines 81-83.
    - Old: `(a) at all n ≡ 0 (mod 4) is therefore implied by the skew-Hadamard conjecture, and is far weaker. A direct
      construction would be a small, honest step in the direction of that open problem…`
    - New: `(a) is elementary (substitution into a regular tournament preserves equicentricity), so it carries no
      information about the skew-Hadamard conjecture; (b) is the open content.`
20. **credit.** Add Nozaki–Suda 2012 (Theorem 2.7: `λ = 1` tournaments of order `4k+2` with an extra spectral condition
    are equivalent to DRT(`4k+3`)).

Note (`05-knowledge/results/sevensphere_anticomplementary_signature_simplex_20261107.md`)

21. **error.** Lines 85-86.
    - Old: `does not continue down a flag of faces for n ≥ 4 (residuals 0.12–0.19).`
    - New: `does not continue down a flag of faces for a generic simplex, n ≥ 4 (exact rank test; consistent on a real
      hypersurface at n = 4, audit J1).`
22. **typing.** Lines 87-88.
    - Old: `a numerical scan of ten candidate coincidences found none … It stays open.`
    - New: `It stays open (no scan is recorded).` Alternatively, commit the scan script.
23. **overclaim (non sequitur).** Lines 101-102.
    - Old: `So the pasted phrase "complementary to V1V2V3" names the right point over GF(2), the field where Fano's
      seventh line is straight.`
    - New: as in correction 4.
24. **overclaim.** Lines 110-111.
    - Old: `The Fano tournament is therefore exactly where the seven-point sphere does not exist. More generally the
      sphere sees precisely a tournament's deviation from regularity.`
    - New: `So the seven-point sphere does not exist for QR7, as for every equicentric tournament (all regular ones, every
      R_a ⇒ R_b, every skew-Hadamard one) and for some non-equicentric ones (source ⇒ C3 ⇒ sink). It needs
      G_T ≠ O_T plus THM-4742's non-coplanarity, not irregularity as such.`
    - Also `all its centres coincide` → `O_T = G_T and O'_v = G'_v`.
25. **cosmetic/credit.** Lines 133-134.
    - Old: `By Reid–Brown, the skew-Hadamard conjecture (OPEN) says`
    - New: `By definition, the skew-Hadamard conjecture (OPEN) says … (Reid–Brown gives the doubly regular form at n − 1)`
26. **error (status).** Lines 136-140.
    - Old: `Proved so far: n ≡ 2 (mod 4); every n ≡ 0 (mod 4) with a skew-Hadamard matrix … Open at the remaining orders.`
    - New: `(a) is proved at every even n ≥ 6 (audit J1); (b) is open.` Line 177's obligation is then moot.
27. **overclaim.** §4 row 1.
    - Old: `In two of the three cases the obstruction was a misreading or a characteristic.`
    - New: `In two of the three cases the obstruction disappears under another reading (C*) or another characteristic.`
28. **cosmetic.** §4 row 2.
    - Old: `gives a cube of centres.`
    - New: `gives a second skew quadrilateral D, G, V0, O'_W sharing the corner G (h fixes G).`
29. **typing.** §4 row 3. Status `PROVED` → `KNOWN (THM-1950), restated`; add `(real rapidity only when s < 1)`.
30. **cosmetic.** §4 row 4.
    - Old: `Equicentric ⊂ equilateral ⊂ skew-Hadamard (as conditions on I + MᵀM)`
    - New: `skew-Hadamard ⟹ equilateral ⟹ equicentric (equilateral = doubly regular or skew-Hadamard)`
    - Column 3: `HYP-9500(b) is the open part.`
31. **overclaim.** §4 row 6.
    - Old: `A new tournament invariant attached to each vertex by the sphere.`
    - New: `A classical quantity: Δ_v = 2SB(v) + (n−1)d⁻(v) + d⁺(v) − C(n,2) (SB = Sonneborn–Berger score); equicentric ⟺
      Δ_v constant.`
32. **§5 ledger.** Update the rows to match corrections 1, 12, 18 and 10:
    - `No flag extension … | NEGATIVE (numerical)` → `NEGATIVE, generic (exceptions on a real hypersurface at n = 4)`;
    - `Equicentric … | PROVED (new)` → `PROVED (new as stated; ingredients THM-1475(1), cbx0611s3 kernel lemma;
      framework Nozaki–Suda)`;
    - `Equilateral … | … (Reid–Brown)` → add `, THM-472`;
    - `HYP-9500 | OPEN (…)` → `(a) PROVED at every even n ≥ 6 (audit J1); (b) OPEN`;
    - `Eight-point sphere of the signature simplex … | PROVED` → add `(at vertices in the regime)`.
33. **credit.** §0: add THM-1475, THM-472, the cbx0611s3 reflection and Nozaki–Suda 2012.

INDEX line (`05-knowledge/hypotheses/INDEX.md`, line 2)

34. Make these replacements:
    - `audit PENDING` → `audit J1 PASS WITH CORRECTIONS`;
    - `no flag extension (n = 4..6, NEGATIVE)` → `no flag extension for generic simplices (n = 4..6; exceptions on a real
      hypersurface at n = 4)`;
    - `(the Fano tournament is where the seven-point sphere degenerates)` → `(QR7, like every equicentric tournament, has
      no seven-point sphere)`;
    - `the eight-point sphere transfers with` → `the eight-point sphere transfers at vertices in the regime, with`;
    - `(proved for n = 2 mod 4 and skew-Hadamard orders; open otherwise)` → `((a) PROVED at every even n >= 6 by
      substitution R[T2] and skew-Walsh towers, audit J1; (b) lambda-progressions OPEN)`;
    - `Builds on THM-4742 (S49), session23b, THM-1950, Reid-Brown.` → add `, THM-472, THM-1475, Nozaki-Suda 2012.`

## Genuine mathematical errors (MISTAKE entry)

1. **THM-4910(6)** (also note §1 and INDEX) asserts that no sphere has the flag corner powers for `n >= 4`. This is
   false on a real hypersurface of nondegenerate 4-simplices, where a real flag sphere exists. The claim is true only
   generically. A "worst residual" over random samples cannot support a universal NEGATIVE; an exact rank test plus
   a search for sign changes can.
2. **THM-4911 title.** "The eight-point sphere transfers to every non-equicentric tournament" is false. Source ⇒ C3
   ⇒ sink has no vertex in the regime; there are 2 such classes at `n = 5` and 6 at `n = 7`.
3. **HYP-9500(a)** is typed OPEN, as a "shadow" of the skew-Hadamard conjecture, but it is a short theorem. Substitution
   of an equicentric tournament into a regular one preserves equicentricity and makes it strong; powers of 2 are covered
   by the skew-Walsh tower already in THM-447. Lesson: before calling an existence question open, try the standard
   compositions (substitution, ordered sum, doubling).
4. **Logic slip** (STRUCTURAL remark, THM-4910(5), note §2). "C* ≡ {1,2,3} mod 2, so 'complementary to V1V2V3' names
   C*" inverts the complement: `{1,2,3}` is complementary to `V0`.

The overclaims (the Fano tournament as "the one"/"exactly where"; "precisely deviation from regularity"; "a new
invariant") and the missing credits are not mathematical errors.

## F. Not checked, and limits

- Literature priority was checked only briefly, through web search and abstract/HTML pages; no PDFs were read. Reid-Brown
  was not consulted in the original. The novelty of (3a)-(3d), of the `c`-dichotomy and of the seven-point sphere is
  therefore "not found", not "established".
- Whether the flag system has exceptions at `n = 5, 6` (codimension 2 and 3) was not determined.
- The bottom quadruple was shown only to be not identically coplanar.
- HYP-9500(b) was not tested beyond `n = 10` (no local search at `n = 12+`).
- The pasted text itself is not in the repo. The reading of "(C)" was judged only against the wording quoted in the
  packet.
- Reproducibility of the producer's scripts: `eight_point_sphere_exact.py`, `flag_negative.py`, `switching_paley.py`
  and `signature_simplex.py` were re-run, with outputs to scratch only. All four exit 0 and reproduce their committed
  `.out` files exactly (timing stamps aside).

## Audit files

`04-computation/experiments/sevensphere_anticomplementary_20261107/audit_J1/`:
- `j1_geometry_symbolic.py`: THM-4910 (1)-(4) symbolically for all `n`, with `S2` from the corner powers (2 min).
- `j1_geometry_exact.py`: exact coordinates for `n = 3..7`; signs; coincidences inside the regime.
- `j1_flag_exact.py`, `j1_flag_locus.py`: statement 6, the exact rank test and the real consistency hypersurface at
  `n = 4`.
- `j1_tournament_census.py`: all classes `n <= 10` (28 s, 0.47 GB); `j1_labeled_crosscheck.py`, labeled check
  `n <= 7`.
- `j1_pfaffian.py`: (3a) kernel and parity.
- `j1_transfer_exact.py`, `j1_regime_failures.py`: (6) and the title counterexamples.
- `j1_hyp9500.py`: the switching construction, the substitution lemma, the construction at every even `n <= 96`, and
  the (3c) construction.
- `j1_misc_checks.py`: the velocity/coth line, `s(DT) = s(T)`, `n = 2`, and the `Δ_v` and Sonneborn-Berger identities.
