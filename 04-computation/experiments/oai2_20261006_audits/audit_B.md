# Audit B — THM-4560, MISTAKE-577, crossing numbers, u(21), section 5, HYP-9215

Auditor: independent adversarial audit (read-only on the repo). Worktree `math-wt-chessboard-20261006` at `9258bf9d52` (every quoted sentence below was re-grepped in the committed files).
All checks use my own code in this folder (`auditB_*.py` + `.out`). The session's `oai2_20261006_fields_ramsey_checks.py` was read only to compare encodings. OpenAI material was read as text (LaTeX/PDF of #164), never run.

## Summary of verdicts

| # | item | verdict |
|---|---|---|
| 1 | THM-4560 (ultrafilter proof, corollaries, finite data) | **PASS WITH CORRECTIONS**: the proof is correct; the THM-521 D claim is too broad; some wording |
| 2 | MISTAKE-577, THM-470 correction, THM-521 update | **PASS WITH CORRECTIONS**: the cutoffs 4 and 5 are confirmed; one false THM-470 line is left uncorrected; THM-521 D's scope |
| 3 | crossing numbers (§1, THM-922/913 updates) | **PASS WITH CORRECTIONS**: all counts are confirmed; one wrong arXiv number; one missing reference |
| 4 | u(21) = 57 in the triangular lattice (§2, MISTAKE-575) | **PASS WITH CORRECTIONS**: 57 and 168 are confirmed; the convention for ω is not stated or is inconsistent |
| 5 | §5.3–5.5 (Davenport, 3-smooth Schur, books) | **PASS WITH CORRECTIONS**: everything is confirmed; one citation number is wrong (#164 Cor. 2.4 → 2.6) |
| 6 | HYP-9215 | **PASS WITH CORRECTIONS**: well-posed and consistent with HYP-9211; spot check passed; the analogy is overtyped |

---

## 1. THM-4560 — PASS WITH CORRECTIONS

**Proof, step by step (all correct).**
* **S_ℓ, P^(ℓ).** `S_ℓ` is closed under `+`, since the level of a sum is the minimum of the levels. `P^(ℓ) + S_ℓ ⊆ P^(ℓ)`, and `S_ℓ` is commutative, so `P^(ℓ)` is a two-sided ideal.
  * Its closure is a two-sided ideal of `βS_ℓ`. The right-ideal half follows from continuity of `ρ_q`; the left-ideal half from continuity of `λ_x` for `x ∈ S`, then of `ρ_q`. So the closure contains `K(βS_ℓ)`.
* **Chain.** `cl_{βS_ℓ} S_(ℓ+1) ≅ βS_(ℓ+1)` (same operation), so `q_(ℓ+1)` is an idempotent of `βS_ℓ`.
  * The "minimal idempotent below a given idempotent" theorem gives a minimal `q_ℓ ≤ q_(ℓ+1)`. It applies because `βS_ℓ` is compact right topological, so Ellis gives an idempotent in a minimal left ideal.
  * `q_ℓ ∈ K(βS_ℓ) ⊆ cl P^(ℓ)` gives `P^(ℓ) ∈ q_ℓ`.
  * `≤` is transitive, and all `q`'s live in `βS_1 = βP_n` with one operation. So `q_a + q_b = q_b + q_a = q_min(a,b)`.
* **Galvin–Glazer.** If `A ∈ p = p + p`, then `A* = {x ∈ A : −x + A ∈ p} ∈ p`. Pick `x ∈ A*` and `y ∈ A ∩ (−x + A)`; then `x, y, x+y ∈ A`, where `y = x` is possible.
  * So exactly the theorem's notion of sum-free (`d1 = d2` allowed) excludes `E ∩ S_ℓ` from `q_ℓ`, which gives `F ∈ q_ℓ`.
* **Induction.**
  * I checked `C_j = {y ∈ F : δ_i+…+δ_j+y ∈ F ∀ i ≤ j}`.
  * For each future level `b′ = min(w_(j+2..m))`, `min(a, b′) = min(w_(j+1..m))` is an invariant level. So `C_j ∈ q_min(a,b′) = q_a + q_(b′)`, i.e. `{x : −x + C_j ∈ q_(b′)} ∈ q_a`. Both cases in the file are this one identity.
  * There are finitely many `b′`, and `C_j, P^(a) ∈ q_a`, so the intersection is nonempty. Any `δ_(j+1)` in it lies in `C_j`, so every sum ending at `j+1` is in `F`, and it keeps the invariant.
* **Ruler corollary.**
  * `max_{i<m≤i′} v_2(m)` is the top differing bit `r` of `i, i′`, so `s_(i′) − s_i` has level `n − r`.
  * Points in one block of the top `d` bits share their first `d` coordinates. Each block splits into two sub-blocks that differ, and strictly increase, at coordinate `d+1`.
  * So the prefix tree is binary and the points are the leaves of a binary subgrid in the sense of THM-453 C. A large `s_0 = (M,…,M)` makes them nonnegative.
  * All pairwise gaps are contiguous sums, so the subgrid is independent.
* **König.**
  * For Finv, the classes realizable at size `t` are the lex-positive vectors in `(−t,t)^n`.
  * Realizability of a triple is per coordinate: `spread{0,a,a+b} = max(|a|,|b|,|a+b|)`. So the box version of sum-freeness is exactly triangle-freeness on `[t]^n`.
  * A2's argument is dimension-free. Hence `t_dead(Finv) < ∞` for every `n`.
* **Definitions match.** THM-4560's "gap-determined" (any `E ⊆ P_n`) is exactly the class of Finv-measurable rules of THM-470 (`kp0611` code: `Finv(x,y) = y − x`).

**Finite data, checked independently.**
* `n = 2` gap game: SAT at `t = 2, 3`, UNSAT at `t = 4, 5`.
  * Two solvers, CaDiCaL 1.5.3 and Glucose 4, on a point-triple encoding (`auditB_n2_games`).
  * A solver-free exhaustive check: 75 winning `E` among all `2^12` at `t = 3`, and 0 of `2^24` at `t = 4` (`auditB_gap_bruteforce`).
* `n = 3` witnesses: `finv_witness_n3_t{4,5,6}.txt` have `|E| = 34, 70, 109`. All are triangle-free and hit all `6^7`, `10^7` and `15^7 = 170,859,375` binary subgrids (`auditB_verify_n3_witnesses`, hierarchical checker written independently).

**Errors and overclaims.**
1. **Scope of the THM-521 D claim (title, Cor 2, and everything that copies it).** "the strong-Specker barrier of THM-521 D holds unconditionally" and "**THM-521 D … holds without HYP-2396 and without "invariant witnesses = valuation gradings".**" are proved only for the *fully gap-determined* reading.
   * THM-521 D's premise "the `p=2` grading runs out at the linear wall `t=2n+1`" is true at `n = 2` only for THM-453 G's **row-invariant** dyadic family (`B_g` through `v_2(g)`, arbitrary column relations). Its cutoff is 5 = 2n+1, which I re-verified: SAT at 4, UNSAT at 5 (`auditB_n2_more`).
   * The gap-determined bare 2-adic algebra (sign, `v_2`) per coordinate has cutoff **4 = 2n** at `n = 2`. It is SAT at 3 and UNSAT at 4, by SAT and by exhaustive enumeration of all `2^12` class sets (`auditB_F2_bruteforce`; FINITE-EXACT, new).
   * So THM-521 D's "translation-invariant (valuation-grading)" is at least partly the row-invariant notion, which THM-4560 does not cover.
   * Corrected text: "…holds unconditionally for fully gap-determined witnesses (in particular every valuation grading of the gap vector). For row-invariant witnesses (THM-453 F/G, including the dyadic `B_(v_2(g))` family) it remains conditional (HYP-2396) / open."
2. **Cor 3.** "Sum-freeness is exactly what keeps `E` out of every idempotent ultrafilter."
   * Sum-freeness is sufficient, not exact. A set lies in no idempotent iff it contains no IP set `FS(⟨x_k⟩)`.
   * Corrected: "Sum-freeness (already the absence of one Schur triple, `d1 = d2` allowed) keeps `E` out of every idempotent ultrafilter."
3. **Notation clash, the very confusion MISTAKE-577 is about.** THM-4560 uses `Q_inv(n,t)` for the Finv game without defining it. THM-453 G uses `invQ(n,t)` for the row-invariant game, and these are different games.
   * Define `Q_inv` in the Setting, or better rename it (e.g. `Q_gap(n,t)`), and add the warning.
4. **Status line.** "Hindman-Strauss Ch. 1-4 facts cited", but Thm 5.8 (Galvin–Glazer) is in Ch. 5. Write "Ch. 1–5".
   * I could not check the numbers Thm 1.60 and Cor 4.18 against the book. Their content is standard (the closure of a left/right ideal is a left/right ideal; a minimal idempotent lies below any idempotent).
5. **Missing cross-reference in the finite data.** "`Q_inv(3,t)` is SAT for `t = 4, 5, 6`" was already implied in the repo.
   * THM-470 B has F2J, a gap-determined algebra, SAT at (3,4), (3,5), (3,6). By A3, Finv is SAT there; the `kp0611` docstring says "F2 SAT there already implies it".
   * Say the new runs are confirmations.
6. **Optional prior-art remark.** The device of a chain of minimal idempotents with `q_a + q_b = q_min(a,b)` is standard in Milliken–Taylor-type and variable-word proofs (e.g. Bergelson–Blass–Hindman 1994). The theorem's statement is repo-specific; no novelty problem.

**Note §5.2 (same material).**
* "*The idea comes from Hindman's own original method, not from #164's.*" is a misattribution. Hindman's 1974 proof is combinatorial; the idempotent method is Galvin–Glazer's, which THM-4560 itself says correctly. Corrected: "The method is the Galvin–Glazer idempotent-ultrafilter proof of Hindman's theorem (not Hindman's original combinatorial proof, and not #164's)."
* "The 2-adic seam is exactly what keeps `E` out of every idempotent" is garbled. Corrected: "Sum-freeness keeps `E` out of every idempotent; the 2-adic seam (one source of sum-free gradings, THM-469) can delay the death of a gap-determined witness but not prevent it."
* "`(3, 7)` is undecided after three CaDiCaL runs of 15–32 minutes" does not match the committed logs:
  * `finv_3_7.log` is a 3600 s TIMEOUT (60 min);
  * `finv_3_7_b1024.log` and `finv_3_7_sym.log` stop at 1138 s and 842 s with no final line.
  * Corrected: "…after the reader's runs (one 60-minute timeout, two runs stopped at about 14 and 19 minutes)…".
* "**THM-521 D's strong-Specker barrier holds unconditionally**", the Outcomes row "THM-521 D unconditional", and §8's "THM-4560 makes THM-521 D unconditional" all need the qualifier from error 1.

## 2. MISTAKE-577 and the THM-470 / THM-521 blocks — PASS WITH CORRECTIONS

**Definitions, checked against the original code.**
* THM-453 F/G is **row-invariant**.
  * `erdos592_invariant_quotient_macmini_s1.py` and `erdos592_satverifier_frontier_macmini_s2.py` both use `key = (a'−a, x[1:], y[1:])`, and `(0, sorted suffixes)` within a row.
  * So this is translation in the top coordinate only, with column dependence arbitrary.
* THM-470's Finv is **fully gap-determined**: `erdos592_invariant_wall_kp0611.py`, `Finv(x,y) = y − x`.
* Gap-determined ⊂ row-invariant (`R(c,c′) = [(0,c′−c) ∈ E]`, `B_g(c,c′) = [(g,c′−c) ∈ E]`), so the row-invariant cutoff is at least the gap-determined one.

**Cutoffs (own encodings, two solvers; `auditB_n2_games.out`).**

| game on `[t]^2` | t=3 | t=4 | t=5 | cutoff |
|---|---|---|---|---|
| free (THM-453 E) | SAT | SAT | UNSAT | 5 (= R(2,2), matches) |
| row-invariant (THM-453 F/G) | SAT | SAT | UNSAT | **5** |
| fully gap-determined (Finv) | SAT | **UNSAT** | UNSAT | **4** |
| gap-determined (sign, `v_2`) | SAT | UNSAT | — | 4 (new) |
| row-invariant dyadic `B_(v_2 g)` (THM-453 G) | SAT | SAT | UNSAT | 5 (matches THM-453 G) |

The session's part 3 encodes the same constraints (I checked its `solve_game`; its unused inner `rec` in `binary_subgrids` is dead code).

**Errors.**
1. **An uncorrected false sentence in THM-470.** The Honesty section says "The invariant = free cutoff equality is proved only at n = 2 (THM-453 G)". Read with THM-470's own "invariant" (= Finv), this is false: Finv has cutoff 4 and the free game 5. The correction block does not mention it.
   * Add: "This equality is for THM-453 G's row-invariant game only; for Finv the `n = 2` cutoff is 4 < 5."
2. **MISTAKE-577 says "Both are SAT, so no false conclusion followed"; THM-470's block says "so nothing false followed".** Too strong. The identification also fed:
   * THM-470's Honesty line above;
   * the `kp0611` docstring's "R(2,2)=5 (where invariant = free, THM-453 G)", cited as support for a `2n+1` gap-determined wall (HYP-2396);
   * THM-521 D / HYP-2558's "the `p=2` grading runs out at the linear wall `t=2n+1`", which at `n = 2` holds for the row-invariant family only (gap-determined walls: 4 = 2n).
   * Corrected: "The INV/invQ comparison itself drew no false conclusion. But the identification also produced THM-470's 'invariant = free at `n = 2`' line, and the use of a `2n+1` gap-determined wall as HYP-2396 evidence. For gap-determined rules the `n = 2` wall is 4 = 2n, strictly below the free cutoff 5."
3. **THM-521 update.**
   * "**Corollary D now holds unconditionally.**" needs the scope qualifier (item 1, error 1).
   * "What remains open is HYP-2558: a non-invariant, value-dependent strong witness" should read "a strong witness that is not fully gap-determined (this includes row-invariant ones)".
   * Two sentences of D are left standing:
     * "explains the t=7 wall (the invariant algebra is exactly the sum-free `p=2` grading…)" is false for Finv, and `(3,7)` is undecided;
     * "can use **only the `p=2` grading**" is true only for bare valuation gradings: THM-469 A2/E's leading-digit gradings are sum-free for every `p` and alive at (3,5).
   * Mark both withdrawn or qualified.
4. **Side note, pre-existing.** The historical index has two different HYP-2558 entries: the strong-Specker barrier, and "residue 1 mod 8 … H-spectrum". THM-4560/521 cite "HYP-2558" without disambiguation.

## 3. Crossing numbers — PASS WITH CORRECTIONS

**Checked (own brute force over 4-sets; `auditB_crossings.out`, `auditB_dds_offsets.out`).**
* The DDS contiguous half-split of `K_n` has exactly `Z(n)` crossings for `4 ≤ n ≤ 60`.
  * This holds for both block sizes `⌊n/2⌋, ⌈n/2⌉`, and for every block offset (`n ≤ 20`).
* The parity-bipartite drawing (odd sum classes on `Z_(2m)`, contiguous `⌊m/2⌋ | ⌈m/2⌉` split) has exactly `d_m^2 = Z(m,m)` crossings for `2 ≤ m ≤ 40`. This equals the telescoping sum.
* **Derivation re-done by hand.**
  * Cyclic gaps: `g1+g3 = s_2 − s_1 = 2κ`. With `g1 + g2` odd, the count is `#{(g1,g2)} = κ(m−κ−1) + (κ−1)(m−κ) = 2κ(m−κ) − m = G(κ)`.
  * Bookkeeping: there are `2m` starting points `a`, each 4-set is counted 4 times, and `κ ↔ m−κ` under rotation.
  * Result: **every unordered pair of classes at cyclic distance κ carries exactly G(κ) crossings** (verified for all pairs, `m ≤ 25`). Pairs inside one class carry 0.
  * Same-page pairs at distance κ number `m − 2κ` (verified).
  * Telescoping: with `f(κ) = κ²(m−κ−1)²`, `f(κ) − f(κ−1) = (m−2κ)(2κ(m−κ)−m)`, so the sum is `K²(m−K−1)²` for **every** K. Also `m−K−1 = ⌊m/2⌋`.
* DDS of `K_8` and `K_14`: 18 and 315 crossings, of which 4 and 81 are bipartite–bipartite.
* The class-colouring minimum of (I) equals `Z(m,m)` for `m = 3..8` (1, 4, 16, 36, 81, **144**), extending the repo's `m ≤ 7`.
* (III) for every `n`: class colourings are 2-page drawings, so `≥ ν_2(K_n) = Z(n)` (Ábrego et al.), and DDS attains `Z(n)`. PROVED, correct.

**Attributions (checked online).**
* de Klerk–Pasechnik–Salazar, arXiv:1207.5701 (SIAM J. Discrete Math. 27 (2013) 619–633):
  * Section **5.1** is titled "The DDS construction…". DDS = Damiani–D'Antona–Salemi 1994 (adjacency matrices); Blažek–Koman 1964 gave 2-page `Z(n)` drawings.
  * They follow "the more geometrical viewpoint of Shahrokhi et al."; their Prop. 5 is the `k`-page count.
  * THM-913's prior-art section and MISTAKE-574 are accurate on these points.
* The 2014 paper is de Klerk–Pasechnik–Salazar, *Book drawings of complete bipartite graphs*, **arXiv:1210.2918**, Discrete Appl. Math. 167 (2014) 80–93. It says:
  * Zarankiewicz's drawings "can be easily adapted to 2-page drawings";
  * "Zarankiewicz's Conjecture implies the (in principle, weaker) conjecture `ν₂(K_(m,n)) = Z(m,n)`";
  * Zarankiewicz is verified for `min ≤ 6` and (7,7), (7,8), (7,9), (7,10), (8,8), (8,9), (8,10).
* That paper gives no explicit parity/sum-class drawing and no telescoping count. The parity-bipartite count is not there, as far as its text shows.

**Errors.**
1. **Wrong arXiv number.** The note's Sources say "Ábrego et al., arXiv:1210.2918."; 1210.2918 is the DPS 2014 bipartite paper. Corrected: "Ábrego, Aichholzer, Fernández-Merchant, Ramos, Salazar, *The 2-page crossing number of K_n*, arXiv:1206.5669 (Discrete Comput. Geom. 49 (2013) 747–777)". Add DPS 2014 = arXiv:1210.2918.
2. **Missing reference.** "De Klerk–Pasechnik–Salazar (2014)" (note §1.3, MISTAKE-574, THM-922 update) is never given a reference. Add the one above.
3. **Novelty wording.** "(new proof; the construction was not found in DPS 2012/2014)" should add: "the upper bound `ν_2(K_(m,n)) ≤ Z(m,n)` is classical (Zarankiewicz drawings adapt to 2 pages, DPS 2014); ours is this explicit sum-class form and its count".
4. **Typing.**
   * "bipartiteness forces `k = g_1 + g_3 = 2κ`": `k` is undefined. Write `s_2 − s_1 = g_1 + g_3 = 2κ`.
   * State the per-class-pair statement explicitly; the sketch jumps from gap counts to the sum.
   * "conditional on a formally verified theorem whose statement we checked" should read "on a *claimed* formal verification (not rebuilt here); only the bipartite comparator statement was read".
5. **Optional.** THM-922 (I)'s lower bound is unconditional for `m ≤ 8` (Kleitman; Woodall (7,7), (8,8)), not only `m ≤ 7`.
6. **Pre-existing.** Two files carry `id: THM-922` (`…cyclic-zarankiewicz…` and `…route-A-signoff`). The update was correctly placed in the first.

## 4. u(21) = 57 — PASS WITH CORRECTIONS

**Checked exactly in `Z[ω]`, ω = e^{iπ/3} (`auditB_lattice.out`).**
* The 21 points `(2+ω)p + (3−ω)q` (`p ∈ {0,1,ω}`, `q ∈ {0} ∪ units`) are distinct and have **exactly 57** pairs at squared distance 7. This equals the generic product count `3·7 + 12·3`, so there are no extra coincidences.
  * The next most frequent squared distances are 21 (18 pairs), then 25, 27, 16, 3 (12 each).
* `W6 × W6`: 49 distinct points and exactly **168** pairs (generic count `12·7 + 12·7`), against `3N = 147`.
* `cos∠(α, ᾱ) = Re(α·conj ᾱ)/7 = (11/2)/7 = 11/14`, and `3 − ω = conj(2 + ω)`.
* `u(21) = 57` matches OEIS A186705, `a(21) = 57`. Alexeev–Mixon–Parshall (arXiv:2412.11914) determine `u(n)` exactly for `n ≤ 21`.

**Errors (typing).**
1. Note §2 and MISTAKE-575 never fix ω; THM-431's correction does (`ω = e^(iπ/3)`). With ω = e^{2πi/3}, `{0,1,ω}` is not a unit triangle and `N(2+ω) = 3`. State ω = e^{iπ/3} in both places.
2. Note §7.3: "Isbell's 7-colouring uses the *split* prime `2 − ω` of norm 7". Under §2's convention `N(2 − ω) = 3`. The norm-7 primes are `2 + ω` and `3 − ω`, with `(2+ω)(3−ω) = 7`. Fix the element or state the other convention there. The note also overloads ω (`ω₁, ω₃, …` Polymath rotations in §4; `ω(s)` in §6).

## 5. §5.3–5.5 — PASS WITH CORRECTIONS

**Checked by brute force, no SAT (`auditB_item5.out`).**
* **Davenport law.** The maximum `m` is `p − 1` for `p = 2, 3, 5, 7, 11` at levels `v = 0, 1`.
  * Integer witnesses `(1,1)` (`p = 3`) and `(1,1,1,1)` (`p = 5`); none of size `p` with entries ≤ 30.
  * Proof: same valuation `v` for all subset sums ⟺ the residues `a_i/p^v mod p` form a zero-sum-free sequence. Length ≤ `D(Z/p) − 1 = p − 1`, attained by all-ones.
* **FS∪FP in one level with `m ≥ 2`** forces `2v = v`, i.e. `v = 0`. Correct.
* **Clique number** of `v_p(|x−y|) = v` is `p` (`p = 2, 3, 5`, `v = 0, 1`). Correct.
* **THM-469 A1 = the case `p − 1 ≤ 1`**: correct, with `x = y` allowed on both sides.
* **Primitive `x+y=z` in 3-smooth numbers** up to `10^15`: `(1,1), (1,2), (1,3), (1,8)`. An elementary proof: pairwise coprime forces a 1 among the three, and Gersonides gives the rest.
  * Ω-parity colouring: 0 monochromatic triples (`a = b` included).
  * The 2-colouring {S with even Ω} | {everything else} has 0 monochromatic quartets `{a, b, a+b, ab}` for `a, b ∈ S ≤ 10^7`, with a one-line proof.
* **Schur thresholds** by exhaustive `2^N` enumeration: **5** with `a = b` allowed, **13** with `a ≠ b`. Correct.
* **Book Turán bound.** Red disjoint cliques plus blue complete multipartite, avoiding red `B_(n−1)` and blue `B_n`: the maximum `N` is exactly `max(2n, 3n−3)` for `n = 3..12, 16, 20`.
  * Attained by `(n,n)` or by three parts `(n−1)³`. `n = 100` gives 297 < 398 = 4n − 2. Correct as stated.

**Errors.**
1. **Wrong citation number.** "#164's Cor. 2.4 runs on the same mechanism (`v = 2v mod q` forces `v = 0`)". In #164 (LaTeX and PDF), 2.4 is *Principle 2.4 (Alignment)*. The mechanism ("so `v = 2v` modulo `q` and `v = 0`") is in the proof of **Corollary 2.6** (*Finite interval form with prescribed divisibility*). Corrected: "Cor. 2.6".
2. **Minor.** State the quartet 2-colouring explicitly: colour A = `{x ∈ S : Ω(x) even}`, colour B = everything else.
3. **Minor.** §5.5 is orientation-specific. The swapped orientation (blue cliques, red multipartite) gives `max(2n+2, 3n−6)`, so the overall Turán-type maximum is `max(2n+2, 3n−3)`. This differs only for `n ≤ 4`; irrelevant at `n = 100`.
4. **Not checked.** The Folkman `2^20`/`2^40` claims and "each `a_d` in #164 has a prime factor `> w`".

## 6. HYP-9215 — PASS WITH CORRECTIONS

**Well-posedness.**
* β_P is defined whenever P lies on a Hamiltonian cycle; the reader reports every path of ≤ 5 moves extends (`n = 5..8`).
* λ_P's values 7, 6, 6 are right for **n ≥ 6**, because G_n is triangle-free there: bipartite for even n, and three knight moves cannot sum to 0 mod n for odd n ≥ 7.
  * 1 move: starving costs 7; early closing would need a triangle.
  * 2 and 3 moves: a non-path neighbour of an interior square has exactly one unusable move, cost 6. The endpoint and early-closing options also cost ≥ 6.
  * The status line "PROVED: β_P ≤ λ_P" is correct.
* `genus ≥ 1 + n²/2` is correct for n ≥ 6 (girth 4).
* **Consistent with HYP-9211**: the 0-move case is λ = β = 7 by stars (n ≥ 5). The evidence counts match C(142,5) = 448,072,338 and C(254,5) = 8,468,125,050.

**Spot check, n = 6, P = (0,0)–(1,2) (`auditB_knight_1move*.out`, independent of the reader's `kp_*` code).**
* P lies on a Hamiltonian cycle.
* Deleting the other 7 moves at (0,0) blocks every tour through P (λ_P ≤ 7).
* **CEGAR proof that no ≤ 6 moves avoiding P block**: the hitting-set instance with 412 Hamiltonian cycles through P and `|D| ≤ 6` is UNSAT (CaDiCaL, 531 s). So **β_P = λ_P = 7**, agreeing with the reader's exhaustive C(143,6) = 10,679,057,389-set run.
  * Control: with `|D| ≤ 7` the same code immediately returns the 7-star.
* **14-vertex Barnette remark confirmed.** The cube with one vertex replaced by Q3 − v (all 6 attachments) is planar, bipartite and 3-connected, with no single blocking edge.
  * It has 57 blocking pairs: 42 stars, **15 non-stars**.
  * Every one of its 84 three-edge paths lies on a Hamiltonian cycle.

**Errors and overclaims.**
1. "In a cubic graph (`k = 3`), #180's corollaries say `β_∅ = k − 1 = 2` and `β_(P_4) = k − 2 = 1`." Also in the note §6.2: "…in the cubic case".
   * These are #180's corollaries for **Barnette graphs** (cubic, bipartite, planar, 3-connected). They are CONDITIONAL on #180, and the `P_4` statement is "in the Pfaffian case" only (note §6.1). Arbitrary cubic graphs need not be Hamiltonian.
   * Corrected: "In a Barnette graph (`k = 3`), #180's corollaries (CONDITIONAL) give `β_∅ = 2`, and in the Pfaffian case `β_(P_4) = 1`."
2. "So the locality here is a consequence of the knight torus's restricted edge connectivity 14 (THM-4552 v)" is unproved. HYP-9211 only shows that no other 7-set kills every **2-factor** (even n ≤ 12). Corrected: "plausibly reflects (heuristic, not proved)".
3. **λ_P definition.** For non-spanning P, the move joining P's two ends is also unusable. Add it next to "moves into P's interior do not count". The values are unchanged.
4. **Evidence sizes.** "`4.48·10^8` sets per class at `n = 6` and `8.47·10^9` at `n = 8`" are the 2-move figures. 3-move classes are C(141,5) ≈ 4.32·10^8 and C(253,5) ≈ 8.30·10^9. The n = 6 1-move class needed C(143,6) ≈ 1.07·10^10 six-sets. Say so.
5. **Name the graph.** "a 14-vertex Barnette graph with a nontrivial 3-edge cut" should be named: the cube with a vertex replaced by Q3 − v.
6. **Optional.** The title says "at most 3 moves" but lists values for 1–3 only. Add "7 for 0 moves (HYP-9211)".

---

## Required corrections

1. **Scope of THM-521 D** (THM-4560 title and Cor 2, THM-521 update, note §5.2 bullet, Outcomes row, §8 #164 row). Replace "THM-521 D holds unconditionally" with "THM-521 D holds unconditionally for fully gap-determined witnesses (all valuation gradings of the gap vector); for row-invariant witnesses (THM-453 F/G, including the dyadic `B_(v_2 g)` family, `n = 2` cutoff 5 = 2n+1) it remains conditional on HYP-2396 / open".
2. **New FINITE-EXACT datum, THM-521 / THM-470.** The gap-determined (sign, `v_2`) algebra at `n = 2` is SAT at `t = 3` and UNSAT at `t = 4` (cutoff 4 = 2n; exhaustive). So D's premise "the `p=2` grading runs out at the linear wall `t=2n+1`" holds at `n = 2` only for the row-invariant dyadic family. Mark D's "explains the t=7 wall / invariant algebra is exactly the `p=2` grading" as withdrawn. Qualify "can use only the `p=2` grading" as "only the bare valuation `p = 2`; leading-digit gradings are sum-free for every `p` (THM-469 A2/E)".
3. **THM-470 correction block.** Add a fix of the Honesty line "The invariant = free cutoff equality is proved only at n = 2 (THM-453 G)": true for the row-invariant game only; for Finv the cutoff is 4 < 5.
4. **MISTAKE-577 and THM-470's block.** Replace "Both are SAT, so no false conclusion followed" / "so nothing false followed" by the sharper statement in §2, error 2 (the identification fed THM-470's Honesty line and the `2n+1` gap-determined-wall evidence for HYP-2396).
5. **THM-4560 Cor 3.** "Sum-freeness is exactly what keeps E out" becomes "Sum-freeness keeps E out (exact criterion: E contains no IP set)". In note §5.2, replace "The 2-adic seam is exactly what keeps `E` out of every idempotent" accordingly.
6. **Note §5.2.** "The idea comes from Hindman's own original method" becomes "the Galvin–Glazer idempotent-ultrafilter method (not Hindman's 1974 combinatorial proof, not #164's)".
7. **Notation.** Define `Q_inv(n,t)` (the Finv game) in THM-4560's Setting, or rename it (`Q_gap`), with a warning that it is not THM-453 G's `invQ(n,t)` (row-invariant).
8. **THM-4560 status line.** "Ch. 1-4" becomes "Ch. 1–5" (Thm 5.8). Optionally check Thm 1.60 and Cor 4.18 against the book.
9. **THM-4560 finite data.** Note that `Q_inv(3,4..6)` SAT was already implied by THM-470 B (F2J SAT, A3).
10. **Note §5.2.** Fix the "(3,7)" run description: one 3600 s timeout and two runs stopped at 842 s and 1138 s, per the committed logs. Not "15–32 minutes".
11. **Note sources.** "Ábrego et al., arXiv:1210.2918" becomes arXiv:1206.5669 (DCG 49 (2013) 747–777). Add de Klerk–Pasechnik–Salazar, *Book drawings of complete bipartite graphs*, arXiv:1210.2918, DAM 167 (2014) 80–93, wherever "DPS (2014)" is cited (note §1, MISTAKE-574, THM-922 update).
12. **Note §1 item 2.** Add that `ν_2(K_(m,n)) ≤ Z(m,n)` is classical (Zarankiewicz drawings adapt to 2 pages, DPS 2014). Define `k` in "`k = g_1 + g_3`" (it is `s_2 − s_1`), and state that each class pair at distance κ carries `G(κ)` crossings.
13. **Note §1.** "conditional on a formally verified theorem whose statement we checked" becomes "conditional on a claimed formal verification (not rebuilt; only the bipartite comparator read)".
14. **Note §2 and MISTAKE-575.** State ω = e^{iπ/3}. In note §7.3, "split prime `2 − ω` of norm 7" becomes `2 + ω` (or `3 − ω`), since `N(2 − ω) = 3` under §2's convention.
15. **Note §5.3.** "#164's Cor. 2.4" becomes "#164's Cor. 2.6"; 2.4 is Principle 2.4 (Alignment).
16. **HYP-9215 and note §6.2.**
    * "In a cubic graph … #180's corollaries" becomes "In a Barnette graph … (CONDITIONAL on #180; the `P_4` case only for Pfaffian ones)".
    * "is a consequence of the knight torus's restricted edge connectivity 14" becomes a heuristic remark.
    * Add the end-to-end move to the unusable moves in λ_P's definition.
    * Correct the per-class enumeration sizes (1-move class C(143,6) ≈ 1.07·10^10; 3-move classes C(141,5), C(253,5)).
    * Name the 14-vertex Barnette graph (cube with a vertex replaced by Q3 − v).
17. **Optional.** THM-922 (I) is unconditional for `m ≤ 8` (Kleitman; Woodall (8,8)), and the class-colouring minimum at `m = 8` is 144 (audit B enumeration). The note §5.5 orientation remark (overall Turán-type maximum `max(2n+2, 3n−3)`). Name the HYP-2558 number collision in the historical index.

**What survives unchanged:**
* THM-4560's theorem and proof, Corollary 1, and `t_dead(Finv) < ∞` for every `n`.
* The `n = 2` cutoffs 4 (gap-determined) and 5 (row-invariant) and the `n = 3` witnesses.
* Every crossing count and the telescoping identity, and THM-922 (III) for all `n`.
* `u(21) = 57` in the lattice, and the 49-point/168-pair set.
* The Davenport law, the clique number, the 3-smooth thresholds 5 and 13, the primitive solutions, and `max(2n, 3n−3)`.
* HYP-9215's λ values and the `n = 6` one-move case (`β_P = λ_P = 7`).
