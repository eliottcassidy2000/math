# A second reading of openai/math: crossing numbers, lattices, Littlewood, plane colouring, Ramsey, Hindman, Barnette, cancellation (2026-10-06)

**Session:** mac-mini-2026-10-06-oaimath2.
**Owner prompt:** "spend another similar session picking out another handful or two of the papers from that repo and looking for connections and extensions and solutions to problems".

**What was read.** github.com/openai/math, read-only through the GitHub API. No OpenAI code was run or built. Twelve manuscripts, in six parallel lanes:

| lane | manuscripts | our fronts |
|---|---|---|
| crossing | #165 crossing numbers of `K_(m,n)` and `K_n` (Lean) | THM-913, THM-922 book drawings |
| lattice | #090 triangular universal optimality (Lean); #022 weak inhomogeneous Duffin–Schaeffer | THM-412/431 unit distances; LRC |
| Littlewood | #076 ultraflat Littlewood polynomials and two companions; #155 aperiodic monotile in `Z³` | the skew-Hadamard tower (HYP-9162 / THM-4557); THM-4551 |
| plane | #158 the plane is not 5-colourable (Lean); #172 classification of Euclidean Ramsey sets (Lean) | Hadwiger–Nelson field tower (HYP-2276–2278); THM-440 |
| Ramsey | #164 Hindman's FS∪FP conjecture; #189 cycle–clique Ramsey numbers; the sharp log-exponent Ramsey pair | Erdős 592 (THM-453/469/470/521); book Ramsey; 3-smooth numbers |
| Barnette | #180 Barnette's conjecture; #047 a counterexample to affine cancellation in dimension 4 | knight tours (THM-4550/4552, HYP-9211); Jacobian/Dixmier (THM-1300, THM-3605) |

**Typing.**
* PROVED (a proof is given or cited);
* FINITE-EXACT (exact computation);
* NUMERICAL;
* CONDITIONAL (rests on an unrefereed openai/math manuscript; "Lean" means its release ships a formalization whose comparator statement we read but did not build);
* ANALOGY / NUMEROLOGY;
* KNOWN (the result is in the literature; credit given).

## Outcomes at a glance

| item | type | where |
|---|---|---|
| THM-922 (III): the class-colouring minimum for `K_n` is `Z(n)` for every `n` | PROVED (DDS count + Ábrego et al. 2012) | §1, THM-922 |
| The parity-bipartite drawing of `K_(m,m)` has exactly `Z(m,m)` crossings | PROVED (new telescoping proof) | §1, THM-922 |
| The 2-page Zarankiewicz conjecture `ν_2(K_(m,n)) = Z(m,n)` | CONDITIONAL on #165 (Lean) | §1 |
| THM-913's drawing is the classical DDS construction | correction (MISTAKE-574) | §1 |
| The triangular lattice attains `u(21) = 57` (THM-431 said 47) | FINITE-EXACT correction (MISTAKE-575) | §2 |
| The skew tower's rows are Thue–Morse-like; Rudin–Shapiro switching flattens all rows to `≤ 8.24√N` | PROVED (THM-4561) | §3 |
| Explicit cyclotomic replacement for #155's random stacking partition | PROVED (its use in #155 CONDITIONAL) | §3 |
| Residue colourings of field planes (Lemma R); exact `χ` for the HN field tower | PROVED (THM-4558); technique KNOWN (Madore 2015) | §4 |
| Our Heegner roadmap for `χ(R²)` is false; "measure cannot reach 5" is false | correction (MISTAKE-576) | §4 |
| Every algebraic spherical set is Euclidean Ramsey | CONDITIONAL on #172 (Lean); deduction PROVED (THM-4559) | §4 |
| Every gap-determined triangle-free graph on `N^n` has an independent binary subgrid; THM-521 D unconditional | PROVED (THM-4560) | §5 |
| THM-470's Finv is not THM-453 F's game (cutoffs 4 vs 5 at `n = 2`) | correction (MISTAKE-577) | §5 |
| Davenport law for FS-sets in one `v_p`-level; 3-smooth Schur thresholds 5 and 13 | PROVED / FINITE-EXACT | §5 |
| Prescribed-path blocking law on knight tori | HYP-9215 (FINITE-EXACT evidence) | §6 |
| Slice lemma; THM-1300's counterexample is fully non-slice | PROVED (THM-4562) | §6 |
| `D(X) ≅ A_4` for #047's exotic 4-space? | HYP-9216 (open question) | §6 |
| The rational plane reduces mod 7 onto the 7 × 7 knight torus | PROVED (one line) | §7 |

## 1. Crossing numbers: THM-913 and THM-922 meet the Harary–Hill and Zarankiewicz theorems (#165)

**The papers.** openai/math #165, *The crossing number of complete bipartite graphs* and *The crossing number of complete graphs* (2026-09-23), claim:

    cr(K_(m,n)) = d_m d_n,   d_r = floor(r/2) floor((r−1)/2)        (Zarankiewicz),
    cr(K_n) = Z(n) = (1/4) floor(n/2) floor((n−1)/2) floor((n−2)/2) floor((n−3)/2)   (Harary–Hill).

* **The lower bound is linear algebra.** If subspaces `U_i, V_i` of `R^N` satisfy `U_i ∩ V_i = 0`, then `Σ_(i,j) dim(U_i ∩ V_j) <= floor(m^2/4) N`. At `N = 1` this is `|A||B| <= floor(m^2/4)` for disjoint `A, B ⊆ [m]`, the source of the factor `d_r`.
* **The detector.** A linear map sends the comparison-sign matrix of every linear order to an invertible matrix, and kills row-plus-column terms.
* **The bridge to drawings.** The signed crossing matrix of a drawing, plus half the incident-order signs at shared endpoints, is block-additive because closed plane curves have zero algebraic intersection. Its blocks then have ranks bounded by crossing counts.
* **Upper bound for `K_n`.** A two-page "endpoint-sum" drawing.

Both papers are unrefereed. Everything below that uses their lower bounds is typed CONDITIONAL.

**Their status is stronger than most of the collection.**
* Both results are claimed **Lean-formalized**: `OAI.Zarankiewicz.mainTarget_proof` in `lean/OAI/Combinatorics/Crossing` (339 files) and `OAI.Paper170.complete_graph_crossing_number` in `lean/OAI/Combinatorics/CompleteCrossing` (46 files).
* We read the comparator statement `ComparatorChallenges/BipartiteCrossing.lean`. It is faithful: admissible drawings by injective continuous paths, finitely many proper crossings, no triple points; crossings counted as points; both the attaining drawing and the universal lower bound.
* A textual scan of all 385 files finds no `sorry`, `axiom` or `admit`. We did not rebuild them.
* So "CONDITIONAL" here means conditional on a formally verified theorem whose statement we checked, not on an unchecked manuscript.

**Prior art for THM-913 (found this session).** The upper-bound drawing of the `K_n` paper is our THM-913's "parallel-class book drawing" exactly. Edge `ij` lies in the matching `M_s`, `s = i + j mod n`, and the pages are contiguous blocks of matchings.
* This is the classical **DDS construction**: Damiani–D'Antona–Salemi 1994, after Blažek–Koman 1964, in the geometric form of Shahrokhi–Sýkora–Székely–Vrt'o.
* It is described in de Klerk–Pasechnik–Salazar, arXiv:1207.5701, Section 5.1, which also computes its `k`-page crossing count.
* THM-913's crossing-profile proof is ours, but the construction is not. This is recorded in THM-913 and in MISTAKE-574.

**What closes (PROVED unless typed).**
1. **The DDS drawing has exactly `Z(n)` crossings for every `n`, even included.** This is classical (DPS 2012), re-proved in the `K_n` paper, and checked here directly for `3 <= n <= 120`.
   * With Ábrego–Aichholzer–Fernández-Merchant–Ramos–Salazar (2012: `ν_2(K_n) = Z(n)`), THM-922 (III)'s class-coloring minimum is `Z(n)` for every `n`, **unconditionally**. THM-922 (III) had only `n <= 14`.
2. **The parity-bipartite drawing has exactly `Z(m,m)` crossings for every `m` (new proof; the construction was not found in DPS 2012/2014).**
   * Put `Z_(2m)` on the spine with parts = parities. `K_(m,m)` is then the union of the `m` odd sum classes. Split them contiguously into two pages, `floor(m/2)` classes and the rest.
   * A crossing is a 4-set whose two alternating diagonals are both bipartite and lie on one page.
   * With cyclic gaps `g_1..g_4`, bipartiteness forces `k = g_1 + g_3 = 2κ`. There are `G(κ) = 2κ(m−κ) − m` gap choices, and the number of same-page residues is `m − 2 min(κ, m−κ)`.
   * So the count is a telescoping sum:

         c = Σ_(κ=1)^K (2κ(m−κ) − m)(m − 2κ) = K^2 (m − K − 1)^2,    K = floor((m−1)/2),

     which is `(floor((m−1)/2) floor(m/2))^2 = d_m^2 = Z(m,m)` exactly.
   * The drawing is the DDS drawing of `K_(2m)` restricted to its even–odd edges. So in that drawing exactly `d_m^2` of the `Z(2m)` crossings are bipartite–bipartite; for example 4 of 18 at `n = 8`, and 81 of 315 at `n = 14`.
3. **CONDITIONAL on #165 (`cr(K_(m,n)) = Z(m,n)`).**
   * The **2-page Zarankiewicz conjecture** `ν_2(K_(m,n)) = Z(m,n)` holds for all `m, n`. De Klerk–Pasechnik–Salazar (2014) state it as open, verified for `min(m,n) <= 6` and a few `(7..8, 7..10)` cases. The upper bound is classical, and `ν_2 >= cr`.
   * THM-922 (I)'s class-coloring minimum is `Z(m,m)` for every `m`, attained by the contiguous split (item 2). Enumeration had reached `m <= 7`.
4. **CONDITIONAL on the `K_n` paper.** THM-913's drawing is optimal among all plane drawings, not only 2-page ones.

**Typing of the method's relevance to our other problems.**
* The rank inequality is a linear-algebra lift of `|A||B| <= m^2/4`. It is the kind of statement our tournament machinery uses (bipartite splits of spectra).
* We found no immediate second application. A possible one, not attempted: lower bounds for `k`-page crossing numbers via block ranks on `k` pages.

Script: [oai2_20261006_crossing_books.py](../../04-computation/experiments/oai2_20261006_crossing_books.py). Its sections:
* exact counts;
* the symbolic telescoping identity;
* the DDS restriction;
* exhaustive class-coloring minima for `n <= 11`;
* a numerical sanity test of the rank inequality.

## 2. The triangular lattice: universal optimality (#090), and a correction to our canon at N = 21

**The paper.** #090 (Lean-formalized in the collection; statement read, not rebuilt) claims that the density-one triangular lattice minimizes the lower energy per particle for every completely monotone potential of squared distance. It also settles Sandier–Serfaty and the 2-D Brauchart–Hardin–Saff linear term.

The method is a linear-programming bound with cardinal interpolation and computer-certified sign checks. The node sets are periodic supersets of the shell values: residues mod 12, and mod 36 in the "atomic certificate" version, verified in exact rational arithmetic. There are no modular forms. The representation numbers `r_Q(D)` enter only the value of the minimum.

**Correction to THM-431 (FINITE-EXACT; MISTAKE-575).** Reading #090 against THM-412 and THM-431, the session's reader found that the triangular lattice attains `u(21) = 57`.
* The 21 Eisenstein integers `(2+ω)p + (3−ω)q`, with `p ∈ {0, 1, ω}` and `q ∈ {0} ∪ units`, have exactly 57 pairs at distance `√7`.
* This is the Erdős product triangle × `W_6` at the resonant angle `arccos(11/14)`.
* THM-431's claim (2) ("NOT optimal; max section 47; gap 10") came from a disk-patch-only search. Corrected in both THM-431 files and HYP-2267.
* The same product gives `W_6 × W_6`: 49 lattice points with `168 > 3N` unit distances.
* The reader's searches find the lattice attaining `u(N)` for every `N <= 12` and at `N = 21`, and one edge short for `13 <= N <= 20` (NUMERICAL; OPEN).

**A torus corollary (CONDITIONAL on #090; reader's derivation, not audited).** For `N` points on a torus of area `N`, the periodic energy for every completely monotone potential is at least the triangular lattice's. Equality holds when the period lattice lies inside the triangular lattice.

At `N = 7` the optimum is the `K_7` triangulation of the torus (FINITE-EXACT, reader's check). Its two triangle classes are the two Fano planes. Orienting edges along the three 120° directions gives the Paley tournament `P_7`, with symmetry `F_21`.

**Weak inhomogeneous Duffin–Schaeffer (#022) and the Lonely Runner (typed negative).**
* With finitely many speeds, the lonely set at the threshold is a finite set of points, invisible to almost-everywhere theorems.
* "LRC for almost all speed vectors" is classical (Kronecker; Czerwiński 2012).
* For the deep well `{1..12, 182}` at `14/183`, the bad sets cover everything while their pairwise overlap is 1.33 times the independent value. So the quasi-independence that drives Duffin–Schaeffer proofs points the wrong way (FINITE-EXACT).
* Wall A needs additive (near-AP) structure, where #022 finds multiplicative structure (ANALOGY only).

## 3. Littlewood polynomials (#076 and companions) and the aperiodic monotile (#155), against the skew-Hadamard tower

**The papers (CONDITIONAL; unrefereed).**
* The three Littlewood manuscripts give successively stronger flatness for `±1` polynomials.
  * The first: `min max_(|z|=1) |P|/√N → 1`. Since `1/F ≤ m_N² − 1`, the best merit factor tends to infinity, which disproves Turyn's conjecture. A Morse-shift corollary follows via Downarowicz–Lacroix.
  * The second: two-sided flatness `√N/16 ≤ |P| ≤ (1+η)√N`.
  * The third, #076: ultraflat polynomials, `(1−ε)√N ≤ |P| ≤ (1+ε)√N` for every `N ≥ N₀(ε)`.
* **The method of #076** is not Rudin–Shapiro doubling, Kahane's random signs, Legendre or Golay. It transplants the Kahane / Bombieri–Bourgain chirp construction to real coefficients:
  * an auxiliary real trigonometric polynomial with `‖F‖_∞ ≤ 1+δ`;
  * chirps on disjoint arcs from a Pippenger–Spencer packing;
  * stationary phase, so each index sees one value of `F`;
  * piecewise-quadratic phase gaps with exterior tails `≤ K_δ/N`;
  * Spencer / Lovett–Meka partial colouring to round the coefficients to signs.
  * No rate; `ε ≍ δ^(1/4)√log(1/δ)`.
* **#155** gives a finite `T ⊂ Z³` that tiles but admits no periodic tiling (also `T + [0,1]³` in `R³`). The method is Greenfeld–Tao:
  * a p-adic Sudoku with infinite descent under `m ↦ pm`, for `p > 200`;
  * CRT encodings in `Z² × V`;
  * stacking by a random partition of `Z/q` into classes with `E − E = Z/q`;
  * rigidification to `Z³`.
  * The size of `T` is not stated. The reader's rough estimate is `|T| > 10^(7·10⁷)`.

**The tower meets Littlewood (THM-4561, PROVED; script `oai2_20261006_littlewood_tower.py`, 12 checks).** The doubling `H ↦ [[H, H], [−H^T, H^T]]` of HYP-9162 / THM-4557:
* takes a **Morse step** on rows, `x ↦ x(1 ± z^m)`, up to one flipped entry;
* takes a Rudin–Shapiro-shaped step on columns, `c_j ∓ w r_j`, applied to the anti-aligned pair `(c_j, r_j)`, so it does not flatten.

So the tower is Thue–Morse in effect:
* **Closed form.** The rows are Walsh rows plus sparse dyadic spike trains, `r_i = w_i − 2Σ_l [y ≡ i mod 2^l](−1)^(Σ_(m≥l) i_m y_m)`.
* **Merit factor.** Every row has merit factor `< 2`, by the transfer `R′ = (3/2 + 2εs)R`, `s′ = (1+4εs)/(6+8εs)` and Cauchy–Schwarz.
* **No Golay pairs.** No two rows are Golay complementary.
* **Rate.**
  * Walsh rows have `F ≤ 1/((3/2)^⌊k/2⌋ − 1)`.
  * Tower rows have `max F → 0` (0.236, …, 0.055 for `k = 8..13`) and `sup/√N → ∞`.
  * The best Walsh row decays like `0.905·(4/(1+√17))^k` (NUMERICAL).
* **The positive statement.** Conjugating by the Rudin–Shapiro diagonal `D = diag((−1)^(Σ y_m y_(m+1)))` gives the same doubly regular tournament `T_k`, since `D H D` renormalizes to `H`. In that switching every row is a sum of Davis–Jedwab Golay pieces on dyadic progressions:

      sup_(|z|=1) |row of D H_(2^k) D| ≤ (3 + 2√2)·√(2N) ≈ 8.24 √N,   F ≥ 1/(33 + 24√2).

  * Numerically, the worst row is `≤ 3.0√N` and `F ∈ [0.75, 3.0]` for `k ≤ 12`.
  * So the same tournament is Thue–Morse in its natural signing and √2-Golay in the Rudin–Shapiro signing.
* **FINITE-EXACT.** The tower's `H_16` has no bent switching: no `±1` vector `d` with `|H_16 d| ≡ 4`, while Sylvester's has 448.
* **Open.** Does every `T_k` have a switching with all rows two-sided flat, or ultraflat? Getting one ultraflat row is trivial given #076.

**Paley rows and Barker (PROVED; Turyn–Storer 1961 plus a finite check).**
* A rotation of a Legendre sequence is a Barker sequence exactly for `p ∈ {3, 7, 11}`.
* So `P_7 = T_3` holds a Barker-7 row in cyclic vertex order (`F = 49/6`). In tower order its best row has `F = 0.889`: the merit factor of a tournament's rows depends on the vertex order.
* Bordered Paley rows reach max `F ≈ 6` at `p = 8191` (FINITE-EXACT, consistent with Høholdt–Jensen).

**#155: a deterministic stacking partition (PROVED; its use inside #155 CONDITIONAL).**
* For a prime `q ≡ 1 (mod s)`, the order-`s` cyclotomic classes, with 0 added to `C_0`, all satisfy `E − E = Z/q` once `q − 3s + 1 > (s − 1)(s − 2)√q`.
  * Proof: the Jacobi-sum bound `s²·(t, t) ≥ q − 2 − 3(s − 1) − (s − 1)(s − 2)√q` on the cyclotomic numbers.
* For `s = 2` this is the Paley split `{0} ∪ QR | NQR`, valid for every prime `q ≥ 7`. When `q ≡ 3 (mod 4)` it is the out/in split of a vertex of the Paley tournament.
* Every failure for `s ≤ 8`, `q < 6000` lies inside the bound (FINITE-EXACT). The largest failures are 5, 13, 41, 101, 277, 491, 337.
* This replaces the random partition of #155's stacking step by an explicit one.

**What does not transfer.**
* #076's sequences are one per length, non-algebraic and without a rate. The tower needs a family closed under its dyadic twists, and Golay is the only such family we have; it gives only `√2`-flatness.
* #155's aperiodicity is p-adic (`p > 200`; `p = 2` is excluded). THM-4551's monotile boundary words are Sturmian. That link is ANALOGY only.
* In the plane, translation tiles admit periodic tilings (Bhattacharya; Beauquier–Nivat), so a translation-only Fibonacci monotile cannot exist there.
* Nothing reaches the H-spectrum (THM-1370); the shared "7" (Barker 7 = `P_7`) is NUMEROLOGY there.

## 4. The plane is not 5-colourable (#158) and the Euclidean Ramsey classification (#172), against our Hadwiger–Nelson field tower

### 4.1 The papers (CONDITIONAL; both ship Lean developments, not built here)

* **#158: no proper 5-colouring of `R²` in ZFC, so `6 ≤ χ ≤ 7`.**
  * *Transfer theorem.* A proper `k`-colouring exists iff a weak measurable one does (same-colour unit pairs null). It is proved by averaging over the real-algebraic plane and its algebraic rotations, then a rigidity theorem: an invariant measure giving the continuous characters zero mass is Haar.
  * *No weak 5-colouring.* The steps:
    * circle palettes and a sixty-degree exclusion;
    * at most two labels per palette;
    * local finiteness from `J₀` decay and the irrationality of `arccos(5/6)`;
    * a covering-dimension lemma (Brouwer), which forces a label cycle;
    * angular analysis, ending with a Moser spindle placed by `z ↦ (35+12i)/37·(z + (−290+149i)/250)`.
  * The finite 6-chromatic graph exists only by compactness (de Bruijn–Erdős): no size, no coordinates.
  * *Checked here (session reader, FINITE-EXACT):* the seven-point placement certificate (11 exact unit edges, not 3-colourable, all points strictly inside the region) and the combinatorial core of the five-cycle step.
* **#172: a finite `A ⊂ R^d` spanning `R^d` affinely is Ramsey iff a tensor certificate exists.**
  * The certificate is `P ∈ Mat_(d+1)(F ⊗_Q F)` with `(p_i ⊗ 1)^T P (1 ⊗ p_i) = 0` and `m(P_spatial) = I`, where `F` is the coordinate field.
  * Consequences in the paper:
    * every subset of a transitive set is Ramsey;
    * every `≤ 5` concyclic points are Ramsey;
    * the transcendental Leader–Russell–Walters kites are Ramsey but not subtransitive;
    * nine generic concyclic points are not Ramsey;
    * a twelve-point Liouville configuration is not Ramsey.

### 4.2 Residue colourings of field planes (THM-4558, PROVED; technique KNOWN)

**Lemma R.** If `L ⊂ C` is a number field with `c(L) = L`, and `𝔓` is a prime of `L` fixed by complex conjugation, then:
* every unit vector (`u c(u) = 1`) is a `𝔓`-unit;
* reduction mod `𝔓`, done coset by coset, maps the unit-distance graph on `L` into `Cay(k_𝔓, U_𝔓)`.

| `𝔓` over `L⁺` | bound |
|---|---|
| ramified | 2 (residue characteristic 2) or 3 |
| inert, residue field `F_q` in `L⁺` | `κ(q) = χ(Cay(F_(q²), μ_(q+1)))` |

**`κ(q)` (FINITE-EXACT).** For `q = 2, 3, 4, 5, 7, 8, 9, 11, 16` it is `4, 3, 4, 4, 4, 4, 3, 5, 4`. Two independent SAT computations agree (`oai2_20261006_fields_ramsey_checks.py`, part 1).

**Credit.**
* For `K²` with the form `x² + y²` and an anisotropic residue form, this is **Madore's reduction** (arXiv:1509.07023, 2015). Madore gets `χ(Q(√2)²) = 2`, `χ(Q(√3)²) = χ(Q(√7)²) = 3`, and `4 ≤ χ(Q(√3, √11)²) ≤ 5` (the upper bound via `F_11`), and lists the exact value as open.
* **That value is 4: K. G. Fischer, 1994** (Congressus Numerantium 104). The decalion89 repository (v1.1.0, 2026; Lean) gives a new short proof by reduction at a place over 2. That repository also has `χ(Q(√2, √3)²) = 4` and `χ(Q(√−3, √−11, √−247)) = 5`.
* Gibbs (Polymath16) 4-coloured the *ring* `Z[ω₁, ω₃]` by the same `F_4` homomorphism.
* Our formulation works with an arbitrary conjugation-stable `L` instead of `K(i)`. That is what lets primes over 2 (where `x² + y²` degenerates) and fields not containing `i` be handled uniformly.

**Exact values for our fields (PROVED).**

| plane | χ |
|---|---|
| `Q(√−d)` | `≤ 3`, `≤ 2` if `d ≡ 1, 2 (mod 4)` (field version of THM-418) |
| Moser field `Q(√−3, √−11)`, and `Q(√3, √11)²` | 4 (KNOWN, Fischer 1994) |
| `Q(√−3, √−(4N−1))`, `N ≡ 2 (mod 3)`: the rungs `√−7, √−19, √−43, √−67, √−163` | **3** |
| the same with `N = 3n`, `n` odd and Loeschian (`N = 3` is Moser) | 4 |
| the Heegner compositum `Q(√−3, √−11, √−19, √−43, √−67, √−163)` | **4** (2 is inert in every factor and `Frob_2 = c`) |
| `Q(ζ_6, ω_t : t odd)`, finite subgraphs | `≤ 4` |
| the Polymath field `Q(√−3, √−11, √−15)` | **5** (upper: the prime over 11, `κ(11) = 5`; lower: Heule's graphs, cited) |

* `G_7`, the 7 × 7 knight torus of THM-4552, is `Cay(F_49, μ_8)`, so `χ(G_7) = κ(7) = 4` (§7).
* A Hilbert-90 test (2000 unit vectors `z/c(z)` in `Q(√−3, √−D)`, `D = 7, 19, 43, 67, 163`) finds every residue in `μ_4 ⊂ F_9`, as the lemma predicts (FINITE-EXACT).
* The reader also checked each colouring on exact-arithmetic unit-distance patches (865–1500 points), with no monochromatic edge.

**Corollary.** If a field plane contains a graph with `χ ≥ 4`, then `L/L⁺` is unramified at every finite prime, so `h⁺(L⁺)` is even. For example, the Moser field is the genus field of `Q(√33)`.

### 4.3 Corrections to our canon (MISTAKE-576)

* **The Heegner roadmap is false.**
  * HYP-2277 conjectured "each chromatic step adjoins a class-number-one rotation field; χ = 5 ↦ √−19". But `Q(√−3, √−19)` is 3-colourable.
  * HYP-2278 (4) conjectured "χ = 2 + #Heegner rotations". But all six together give 4.
  * The 5-chromatic step needs the class-number-two field `Q(√−15)` (Heule's `ω₄`).
  * What governs `χ` is how 2, 3, 5, 7 and 11 decompose in `L/L⁺`, not class numbers.
* **"Measure cannot reach 5" (HYP-2278 (1)) is false.**
  * Falconer (1981) proved `χ_m ≥ 5` by measure theory.
  * `m₁(R²) ≤ 0.247 < 1/4` (Ambrus et al., arXiv:2207.14179) gives it by density alone.
  * #158 claims 6 by ergodic and topological means.
  * What survives: **density alone cannot give 6** (PROVED; Croft's `m₁ ≥ 0.2293 > 1/5`). So #158's sixth colour is necessarily topological, and `J₀` enters #158 only through its decay.
* Our spectral floor `3.48` (HYP-2278) would, under #158's transfer theorem, apply to all colourings; it still rounds to 4.

### 4.4 Where #158's graph can live (CONDITIONAL on #158)

* The compactness graph `H` has real-algebraic coordinates (Tarski). By Lemma R its field avoids the Moser, Heegner and Polymath fields and the whole odd-`t` compositum.
* The smallest Polymath-ladder fields that the residue test (`p ≤ 400`) does not exclude are `Frac Z[ω₁, ω₂, ω₃]`, `Frac Z[ω₁, ω₃, ω₄, ω₅]` and `Frac Z[ω₁, ω₂, ω₃, ω₄]`. Their conjugation-stable residue fields give only `κ(17) ∈ [5, 7]` and `κ(29), κ(41) ≥ 5` (Hoffman).
* These are the natural places to search for an explicit 6-chromatic graph.

### 4.5 Algebraic spherical sets are Ramsey (THM-4559, CONDITIONAL on #172; the deduction PROVED)

* For a number field `F`, the separability idempotent `e ∈ F ⊗_Q F` has `m(e) = 1` and `e·(x ⊗ y) = e·(xy ⊗ 1)`.
* So `P = e·(H ⊗ 1)`, with `H` the sphere matrix, is a certificate for every spherical `A` with algebraic coordinates.
* **Graham's spherical conjecture therefore holds for algebraic configurations; every non-Ramsey spherical set is transcendental.** This matches Pálvölgyi's heptagon (arXiv:2609.23327, non-Ramsey for every transcendental radius) and #172's own examples, both driven by derivations, which vanish on number fields.
* We did not find it stated in #172 or in Pálvölgyi's abstract.
* *Our configurations.*
  * Every spherical configuration in the repo's unit-distance work (Eisenstein, Moser, Heegner fields, `P_7` heptagon subsets) is Ramsey.
  * THM-431/440's extremal sets are not spherical, hence not Ramsey (unconditional, EGMRSS).
  * No bearing on `u(22)`.
* *Open, and decidable by the criterion:* 6 generic concyclic points.

## 5. Hindman's FS∪FP conjecture (#164), cycle–clique (#189), sharp log exponents: Erdős 592, 3-smooth numbers, books

### 5.1 The papers (CONDITIONAL)

* **#164.** For every finite colouring of `N` and every `m`, some `a_1 < … < a_m` have all subset sums and products in one colour. A super-separated, finite version holds inside `qN`.
  * The proof uses no idempotent ultrafilters; a nonprincipal ultrafilter appears only as a limit device in an auxiliary parameter.
  * Its finite core (after Alweiss's proof over `Q`) combines Ramsey on index subsets with Folkman–Rado–Sanders on block sizes, giving "product chains".
  * Its analytic core:
    * harmonic variables at separated scales;
    * divisor weights and Cauchy–Schwarz, reducing to additive cube tests;
    * Green–Tao–Ziegler and Tao–Ziegler concatenation, giving nilsequence models;
    * a quantitative Leibman theorem and nilpotent recurrence, giving uniform alignment.
* **#189.** `R(C_m, K_n) = (m−1)(n−1) + 1` for all `m ≥ n ≥ 3` except `(3, 3)`.
  * The method: a minimal counterexample with expansion `|N[I]| ≥ k|I| + 1`, a large clique from BFS layers, path systems on a maximum clique, and 3,099 finite pattern instances (`5 ≤ k ≤ 17`) killed by two implementations.
* **Sharp log exponents (skimmed).** `r(s, t) = t^(s−1)/(log t)^(s−2+o(1))` for fixed `s ≥ 6` (and `s = 5` in a companion). The lower bound uses random ordered incident flags in `PG(d, q)`, the upper bound Ajtai–Komlós–Szemerédi.
  * The only link to our triangle-case front is an ANALOGY: the height-1 row of our tree-grid numbers is `r(3, b)`.

### 5.2 Erdős 592: invariant witnesses always die (THM-4560, PROVED)

* *The idea comes from Hindman's own original method, not from #164's.* Take a chain of minimal idempotent ultrafilters `q_1 ≤ … ≤ q_n` on the level semigroups `S_ℓ` (lex-positive vectors of level `≥ ℓ`), with `P^(ℓ) ∈ q_ℓ`.
* A sum-free gap set `E` lies in no idempotent (Galvin–Glazer). So gaps of any prescribed level word can be chosen with **every contiguous sum outside `E`**.
* With the ruler word `n − v_2(i)`, the prefix sums are the leaves of a binary subgrid whose pairwise gaps are exactly those contiguous sums. **Every gap-determined triangle-free graph on `N^n` has an independent binary subgrid.**
* Consequences:
  * `t_dead(Finv) < ∞` for every `n` (with THM-470 A2);
  * every gap-determined feature algebra dies at a finite `t` (A3);
  * **THM-521 D's strong-Specker barrier holds unconditionally**, without HYP-2396 or "invariant = valuation gradings".
  * A strong witness, if one exists, must be value-dependent (HYP-2558, open).
  * The 2-adic seam is exactly what keeps `E` out of every idempotent: it can delay the death of an invariant witness but not prevent it.
* *Finite data (FINITE-EXACT).*
  * Fully gap-determined at `n = 2`: SAT at `t = 3`, UNSAT at `t = 4` (CaDiCaL here, and the reader).
  * At `n = 3`: SAT for `t = 4, 5, 6`, with `|E| = 34, 70, 109`; witnesses were brute-verified, including all `1.7·10^8` binary subgrids at `t = 6`.
  * `(3, 7)` is undecided after three CaDiCaL runs of 15–32 minutes, plus the repo's 2-hour run. So `t_dead(Finv) ∈ [7, ∞)` at `n = 3`, finite by the theorem.
* *Correction (MISTAKE-577).* THM-470 identified Finv with "the translation-invariant game of THM-453 F", which is only row-invariant. At `n = 2` the row-invariant cutoff is 5 (THM-453 G, recomputed) and the fully gap-determined cutoff is 4.

### 5.3 The Davenport law for FS-sets in one valuation level (PROVED, elementary)

* The largest `m` such that all nonempty subset sums of `a_1, …, a_m` have the same `p`-adic valuation is **`p − 1 = D(Z/p) − 1`**, the Davenport constant of `Z/p` minus one. Taking all `a_i ≡ 1 (mod p)` attains it.
* If FS∪FP lies in one level and `m ≥ 2`, that level must be `v = 0`.
* The graph `v_p(|x − y|) = v` has clique number exactly `p`.
* **THM-469 A1 ("`L_v` is sum-free for every `v` iff `p = 2`") is exactly the case `D(Z/2) − 1 = 1`.** #164's Cor. 2.4 runs on the same mechanism (`v = 2v mod q` forces `v = 0`): any homomorphic grading pushes FS∪FP sets into its identity class.

### 5.4 3-smooth numbers: no Hindman theorem (PROVED / FINITE-EXACT)

* **Primitive solutions.** The primitive solutions of `x + y = z` in the 3-smooth numbers `S` are `1+1, 1+2, 1+3, 1+8` (Gersonides; rechecked below `10^12`).
  * Hence two distinct elements of `S` whose `Ω` has the same parity never sum into `S`.
  * So colouring `S` by `Ω mod 2` leaves no monochromatic Schur triple inside `S`.
* **Three colours.** `Ω mod 2` on `S`, plus a third colour elsewhere, avoids every `{a, b, a+b}` with `a, b ∈ S`.
* **Two colours, quartets.** Two colours avoid every quartet `{a, b, a+b, ab}` with `a, b ∈ S`.
* **Two colours, Schur forced.** Two colours do force a monochromatic `{a, b, a+b}` with `a, b ∈ S`. The least forcing `N` is **5** if `a = b` is allowed and **13** if not (FINITE-EXACT; reader and this session agree).
* **Folkman.** With `m = 3` and 3-smooth generators, Folkman is 2-colourable up to `2^20` (SAT), and an explicit rule works to `2^40` (FINITE-EXACT, reader). Extending it to all `N` needs every solution of `x + y = z + w` in `S`.
* **#164 avoids 3-smooth numbers.** Each `a_d` in #164's construction has a prime factor `> w`.
* **Collatz.** THM-4555's root collision is the relation `1 + 2 = 3`, i.e. `−1/2` is the fixed point of `x ↦ 3x + 1`.
  * Ramsey theory cannot force such `{2,3}`-unit relations; their sparsity is S-unit theory.
  * Parity-vector colourings are residue colourings mod `2^K`, on which Hindman only returns the zero class.

### 5.5 Book Ramsey (PROVED; FINITE-EXACT for `n ≤ 10`)

* **Turán-type colourings fall short.** Take red disjoint cliques and blue complete multipartite graphs, for `(B_(n−1), B_n)`. They reach only `N ≤ max(2n, 3n − 3)`. At `n = 100` that is 297, against the 398 we need. So our (94, 104) search is quasirandom territory that #189's method does not see.
* **The stored witness.** The `n = 12` witness (`N = 46`) verifies. #189's consequences for it (red `C_7`, `C_8`; blue `C_6`–`C_10`) were all already known.

## 6. Barnette's conjecture (#180) and affine cancellation (#047): knight tours and the Jacobian ecosystem

### 6.1 The papers (CONDITIONAL)

* **#180: every cubic, bipartite, planar, 3-connected graph is Hamiltonian.**
  * The work happens in the dual triangulation. Tutte states (black faces choose corners bijectively) come in pairs (Prop 3.1: planar density plus integral max-flow).
  * A disk identity makes the winding `J` an integer. The exponential sum `Z(x) = Σ_pairs i^J e^(xω(s))` cancels the cyclic selections and is nonzero by a double-dimer regrouping.
  * A forest selection then dualizes to a Hamiltonian cycle.
  * No Kasteleyn or Pfaffian orientation appears.
  * Corollaries: every edge is avoided by some Hamiltonian cycle; in the Pfaffian case, every 3-edge path lies on one.
  * *Checked here (reader, weak):* on duals of 7 Barnette graphs, `J` is integral and forest pairs exist. Lemma 5.2's cancellation was tested only vacuously.
* **#047: `X = Spec A`, `A = C[p, s, u, F, J]/(H)`, has `X × A¹ ≅ A⁵` but `X ≇ A⁴`.**
  * Here `x = s² + u³ + p²F` and `H = x²F − (1+2sx)J − p²J² − pu`.
  * `H` is a stable coordinate but not a coordinate, and its fibrations are not locally trivial.
  * The cylinder comes from a characteristic-0 residual-coordinate route: `Δ = −p²∂_u`, `ΔH = p³`. Every identity was checked in sympy by the reader.
  * The non-isomorphism uses a Derksen-type invariant on the `p`-adic associated graded ring. Lifting to `SL_2 × A¹` and applying Mason–Stothers to the (3,2)-cusp `d² + a²u³` gives Makar-Limanov's idea for the Russell cubic.

### 6.2 Knight tours: one blocking law from `k = 3` to `k = 8` (HYP-9215)

* **The dictionary.**
  * `β_P` is the fewest edges, disjoint from a prescribed path `P`, whose deletion kills every Hamiltonian cycle through `P`.
  * `λ_P` is the cheapest local obstruction: starve one square, or force an early closing.
  * #180's corollaries are `β_∅ = λ_∅ = k − 1 = 2` and `β_(P_4) = λ_(P_4) = k − 2 = 1` in the cubic case. HYP-9211 is `β_∅ = 7 = k − 1` at `k = 8`.
* **Data on the knight torus `G_n`** (FINITE-EXACT, reader):
  * `β_P = λ_P` for every class of 1-move paths at `n = 6`, and every class of 2- and 3-move paths at `n = 5..8`;
  * minimum blocking sets are exactly the local ones in three sampled classes at `n = 6` (43, 45 and 84 sets);
  * every path of `≤ 5` moves lies in a closed tour, `n = 5..8`.
  * This is filed as **HYP-9215**.
* **"Stars only" is not Barnette-like.** In a 14-vertex Barnette graph with a nontrivial 3-edge cut, 15 of 57 minimum blocking pairs are not stars. On the knight torus the locality comes from its restricted edge connectivity 14 (THM-4552 v).
  * #180's disk identity cannot transfer either: `genus(G_n) ≥ 1 + n²/2`.
* **On the real 8 × 8 board,** every non-extendable path of 1 to 5 moves (0, 11, 101, 809 and 5123 classes) is caught by iterated local forcing (FINITE-EXACT).
  * The first *global* obstruction is our weave law (THM-4550 W): 9 moves inside rings 1–2 force `k ≥ 9 > 8`.
  * Example: the path `(6,2),(5,4),(3,5),(1,6),(2,4),(3,2),(1,1),(2,3),(4,2),(6,3)` passes every local test and still extends to no tour (PROVED by the weave law).
  * So the weave law is a genuinely global, Pósa-type obstruction.

### 6.3 The Jacobian ecosystem

* **Slice lemma (THM-4562, PROVED).** A Keller map's dual field `V_j` is locally nilpotent iff `C[x] = R[F_j]` with every other `F_i ∈ R`.
  * Then `C^n ≅ X_j × A¹`, and `F` is injective iff the Keller map `X_j → A^(n−1)` is.
  * So JC_n with a slice is **cancellation in dimension `n − 1` plus JC on `X_j`**.
  * For `n = 3`, surface cancellation makes such a counterexample a JC₂ counterexample.
  * For `n = 5`, #047's `X` admits no injective Keller map to `A⁴` (CONDITIONAL).
* **THM-1300 is fully non-slice (THM-4562, PROVED).** Fibre counts over `F_q`:
  * zero fibres `2q² − 2q + 1`, `q² − q + 1`, `2q² − q`;
  * nonzero fibres `q² − q + 1`, `q² + 1`, `q² − q`.
  * Checked by brute force here and by the reader at `q = 5..17`.
  * Hence no dual field is locally nilpotent and no component is a (stable) coordinate. The counterexample's non-injectivity cannot be pushed down a dimension, which fits JC₂ being open.
* **The same threefold (PROVED, a change of variables).** #047's base `{xy = z(z+1)} × A¹` is THM-3605's Russell cylinder (`B = 4z`, `C = x`, `Y = 16y`). So the stabilized THM-3561 map is an étale, non-injective map `A³ → Spec R`, over the threefold where #047 principalizes its ideal.
* **Differential operators (HYP-9216).**
  * From #047, projective `O(X)`-modules are free. So `T*X ≅ A⁸`, `D(X) ⊗ A_1 ≅ A_5` and `gr D(X) ≅ C^[8]`. K-theory and Hochschild homology agree with `A_4`.
  * **Is `D(X) ≅ A_4`?** Either answer is new: "no" would be a noncommutative cancellation failure; "yes" would give non-isomorphic varieties with the same differential operators.
* **Analogies only.**
  * The (3,2)-cusp sits beside THM-1340's cuspidal cubic and THM-2570's cusp cylinder, with no map between them.
  * Cancellation does not pad, whereas JC failure does.

## 7. Cross-links

**7.1 The rational plane reduces mod 7 onto the 7 × 7 knight torus (PROVED, one line).**
* 7 is inert in `Q(i)`. So by Lemma R every rational unit vector `(a/c, b/c)` (`a² + b² = c²`) has `7 ∤ c`, and reduces to a point of the circle `x² + y² = 1` over `F_7`: `(±1, 0), (0, ±1), (±2, ±2)`.
* Multiplication by `1 + 2i` in `F_49 = F_7(i)` (norm 5) maps this circle onto the knight's moves `(±1, ±2), (±2, ±1)`, the circle `x² + y² = 5` of THM-4552.
* So **`z ↦ (1 + 2i) z mod 7` maps every rational unit-distance pair to a knight move on the 7 × 7 torus**: a graph homomorphism from the unit-distance graph of `Q²` (each coset of the unit-vector group) onto `G_7`.
* The same reduction gives:
  * mod 3: the `3 × 3` rook torus `K_3 □ K_3`;
  * mod 11: a 5-chromatic residue graph.
* The exceptional knight torus of THM-4552 is the residue graph of the plane at the prime 7, which is why `χ(G_7) = κ(7) = 4` appears in both places.

**7.2 Idempotents as certificates.** Three of the session's results turn on an idempotent:
* the separability idempotent `e ∈ F ⊗ F`, which turns sphericity into #172's tensor certificate (THM-4559);
* minimal idempotent ultrafilters on the level semigroups, which choose gaps avoiding a sum-free set (THM-4560);
* reduction modulo a conjugation-stable prime, a ring projection onto a finite field, which turns an infinite field plane into a finite Cayley graph (THM-4558).

#164 proves Hindman's FS∪FP conjecture *without* idempotent ultrafilters, while #172's sufficiency direction uses Ellis's idempotent lemma. In each case the infinite statement is certified by one finite or algebraic object: a residue graph, a tensor, or a level word. That is the "information compression" reading of these results.

**7.3 Sevens (ANALOGY, with PROVED or classical ingredients).** The prime 7 recurs across the session:
* `G_7 = Cay(F_49, μ_8)` (§7.1, PROVED);
* `#090`'s `N = 7` torus is `K_7`, whose two triangle classes are the two Fano planes; oriented at 120°, it is `P_7` with symmetry `F_21` (FINITE-EXACT, reader);
* `u(21) = 57` is attained by the Minkowski product of the norm-7 directions `2 + ω, 3 − ω` at angle `arccos(11/14)` (§2);
* the Legendre/Barker 7 is a row of `P_7 = T_3` (§3), and `Aut(T_k) = F_21` (THM-4557).

Isbell's 7-colouring uses the *split* prime `2 − ω` of norm 7, so it is a lattice quotient, not a Lemma R reduction (split primes are not conjugation-stable). Reading these as one phenomenon is NUMEROLOGY.

**7.4 The 2-adic seam (ANALOGY).** Several 2-adic objects appear in the session:
* Erdős 592's 2-adic seam (THM-469; the Davenport law at `p = 2`);
* the prime over 2 with residue field `F_4` that caps the whole odd-`t` Hadwiger–Nelson compositum at 4;
* the Rudin–Shapiro path form that flattens the 2-adic tower.

In each, a quadratic or valuation structure at 2 is either the obstruction or the certificate. No common theorem is claimed.

## 8. Verdicts by manuscript

| manuscript | verdict |
|---|---|
| #165 crossing numbers (Lean) | **Closes repo items:** THM-922 (III) for every `n`; THM-922 (I) upper bound for every `m` (new proof). The 2-page Zarankiewicz conjecture and THM-913's plane optimality are CONDITIONAL. Prior art found for THM-913 (MISTAKE-574). |
| #090 triangular universal optimality (Lean) | **Correction:** `u(21)` is attained by the lattice (MISTAKE-575). Torus corollary CONDITIONAL. Silent on unit distances (non-CM). |
| #022 weak inhomogeneous Duffin–Schaeffer | **Typed negative** for LRC: finite speed sets are invisible to almost-everywhere theorems; the deep well's overlaps run 1.33 times independent. |
| #076 + companions (Littlewood) | **New theorem** about the tower (THM-4561); the method does not transfer. |
| #155 aperiodic monotile in `Z³` | Explicit cyclotomic stacking partition (PROVED); the monotile link is ANALOGY. |
| #158 plane not 5-colourable (Lean) | **Corrections** (MISTAKE-576) and **exact field-plane values** (THM-4558; technique Madore's). Where a 6-chromatic graph can live is CONDITIONAL. |
| #172 Euclidean Ramsey classification (Lean) | **Conditional extension:** algebraic spherical sets are Ramsey (THM-4559). |
| #164 Hindman FS∪FP | **Method transfer** (the classical idempotent route): THM-4560 makes THM-521 D unconditional and gives `t_dead(Finv) < ∞`. Davenport law. No 3-smooth Hindman theorem. |
| #189 cycle–clique | No help for our book-Ramsey search (Turán-type colourings stop at `3n − 3`). |
| sharp log exponents | ANALOGY only. |
| #180 Barnette | A blocking-law dictionary from `k = 3` to `k = 8`: HYP-9215. The weave law is a global obstruction. |
| #047 affine cancellation | Slice lemma and the non-slice THM-1300 (THM-4562); Russell-cylinder identity; HYP-9216. |

## 9. Audit record

(pending: two independent audits, filled in below before canon)

## 10. Reproduction

Session scripts, each with a `.out` and ending in ALL CHECKS PASSED:
* [oai2_20261006_crossing_books.py](../../04-computation/experiments/oai2_20261006_crossing_books.py): §1.
* [oai2_20261006_lattice_u21.py](../../04-computation/experiments/oai2_20261006_lattice_u21.py): §2.
* [oai2_20261006_littlewood_tower.py](../../04-computation/experiments/oai2_20261006_littlewood_tower.py): §3, THM-4561.
* [oai2_20261006_fields_ramsey_checks.py](../../04-computation/experiments/oai2_20261006_fields_ramsey_checks.py): §4–5, THM-4558, THM-4560, MISTAKE-577.

The six readers' own scripts, logs and witnesses, as run, are in [oai2_20261006_readers/](../../04-computation/experiments/oai2_20261006_readers/), one folder per lane. Their absolute scratch paths may need adjusting. The OpenAI LaTeX sources the readers worked from are not committed: they are public at github.com/openai/math.

**Sources.**
* openai/math manuscripts #165, #090, #022, #076 (+2), #155, #158, #172, #164, #189, #180, #047, accessed 2026-10-06.
* D. A. Madore, arXiv:1509.07023.
* K. G. Fischer, Congr. Numer. 104 (1994).
* github.com/decalion89/chromatic-number-of-the-plane (v1.1.0).
* D. Pálvölgyi, arXiv:2609.23327.
* G. Ambrus, A. Csiszárik, M. Matolcsi, D. Varga, P. Zsámboki, arXiv:2207.14179.
* de Klerk–Pasechnik–Salazar, arXiv:1207.5701.
* Ábrego et al., arXiv:1210.2918.
* Alexeev–Mixon–Parshall, arXiv:2412.11914.
* Hindman–Strauss, *Algebra in the Stone–Čech Compactification*.
* Davis–Jedwab (1999); Turyn–Storer (1961); Falconer (1981); Fujita (1979); Miyanishi–Sugie (1980).
