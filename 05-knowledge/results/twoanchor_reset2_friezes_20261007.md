# Two anchors for the residual first-reset-2 branch; chirotopes, friezes and a pair-independence law

**Session.** mac-mini-2026-10-07-twoanchor.

**Owner prompt.** "work towards a variable-depth rule for the residual first-reset-2 branch, think cluster algebras, total positivity and how friezes and conway-coxeter connect with our tournament adjacent work, again take inspiration from the attached papers but go much further beyond them". The four attached OpenAI preprints are treated as trusted; section 9 says what was used from them.

**Starting point.** Section 4 of [terminal lifts](collatz_terminal_lifts_20261007.md) leaves the residual branch open. Concurrent work on it:
* [reset-two rules](collatz_reset2_rules_20261007.md): guarded K → K−2 and K → K−6;
* [reset-two synthesis](collatz_reset2_synthesis_20261007.md): ternary turns and a ledger;
* [frieze package](frieze_collatz_reset2_20261007.md): positive marked minors, and the ear-move obstruction;
* [child closure](collatz_child_closure_synthesis_20261007.md).

This note is complementary. The aim is uniformity in the second unbounded parameter of the branch, the length of its two-run, which no earlier rule controls.

COLLATZ IS STILL OPEN. Nothing here covers every positive integer.

**Audits.** Two independent audits ran on 2026-10-07: A (Collatz) and B (tournaments, friezes, partner law). The core mathematics was confirmed. Their corrections are applied below and recorded in MISTAKE-586. Reports: `04-computation/experiments/twoanchor_20261007/audit_{A,B}/REPORT.md`.

## 0. Results and types

| Item | Statement | Type |
|---|---|---|
| [THM-4600](../../01-canon/theorems/THM-4600-anchored-debts-are-run-transparent-the-borel-torus-picture-of-pair-chains.md) | Pair-chain states anchored at the letter-a fixed point `c_a = 1/(2^a−3)` pass unchanged through common runs of a. They are the dilations commuting with the run map: a torus in the Borel group. `a = 1` is THM-4555's collision at −1; `a = 2` (the trivial cycle +1) is new in the repo. | PROVED |
| [THM-4601](../../01-canon/theorems/THM-4601-two-anchor-reduction-of-the-residual-first-reset-two-branch.md) (i) | The residual branch is exactly the lag-one translation pair (x, x−1) with x ≡ 1 mod 8. | PROVED |
| THM-4601 (ii) | **+1 barrier.** For D ≤ 8000 no (generalized) deletion collision with a residual source completes by Terras time j after the run end, in particular not inside the two-run. | PROVED for D ≤ 8000 (exact limit chains); HYP-9240 for all D |
| THM-4601 (iii) | **Universal residual state.** If the two-run has J ≥ 3 letters, the K → K−3 and K → K−4 chains end the run exactly in `(3, 1−27)`, i.e. `X − 1 = 27(Y − 1)`, whatever K, t and J are. | PROVED |
| THM-4601 (iv) | **Ladder compiler at +1.** Head collisions with `|v| = |u|+3` and a 4z+1 ladder condition (child: `Σu = Σv+2i`, `F_v(1) = 4^i F_u(1) + (4^i−1)/3`; or source) give certificates uniform in both run lengths K and J. Explicit i = 1 pairs `(10)~(3,1,1,3)` and `(1,10)~(1,1,2,2,3)`. | PROVED (an identity of affine maps; the child's word is forced) |
| THM-4601 (v) | **The variable-depth rule.** For every D the end-of-run state depends only on (D, J) and is constant for J ≥ J_0(D). So certification is a fixed pair-chain problem in the post-run bits, uniform in K, t and J ≥ J_0(D), and its Haar coverage tends to 1 (THM-4581). | PROVED |
| THM-4601 (vi) | Mersenne line: in every shell `v2(K−1) = m ≥ 4` the density certified by K → K−3 within post-run depth s is one number `ρ_r(s)`. `2^1889−1 ⇝ 2^1886−1` and `2^5249−1 ⇝ 2^5246−1` use the same head pair. | PROVED (shell law), FINITE-EXACT (examples) |
| [THM-4602](../../01-canon/theorems/THM-4602-tournaments-as-chirotopes-pfaffian-criterion-and-legendre-frieze-strips.md) | A tournament is a Gr(2,n) sign pattern iff every 4-Pfaffian is ±1 iff it is locally transitive. The transitive tournament is the totally positive part; positive friezes are a slice of it and Conway–Coxeter friezes are its integer points. Over F_p (p ≡ 3 mod 4) the Legendre pattern of P^1(F_p) is Paley + sink up to switching, and Singer strips close as friezes. | PROVED; statements (1), (2), (3)(a)–(b) KNOWN (Babai–Cameron 2000; Gunderson–Semeraro 2017, 2022) |
| [HYP-9241](../hypotheses/HYP-9241-pair-independence-product-law-for-collatz-translation-partners.md) | Pair independence for translation partners: `q_R(T) ≈ ∏_(pairs) q^(lag)(T)`, hence exponent R(R+1)/4. The data reject Karlin–McGregor repulsion. | NUMERICAL |

## 1. The branch, exactly

**Setup.**
* `n = 2^K t − 1`, t odd, K ≥ 2.
* After `K − 1` letters 1 the orbit is at `x = 2·3^(K−1) t − 1`. Write `x = 1 + 2^j u` with u odd.
* First reset 2 means `v2(3^K t − 1) = 1`. That is equivalent to x ≡ 1 mod 8, i.e. j ≥ 3.
* The source then has `J = ⌊(j−1)/2⌋` letters 2: the two-run, ending at Terras time `2J = j − r`. It is followed by letter 1 if r = 1, or by a letter ≥ 3 if r = 2, where `r = j − 2J`.
* On the Mersenne line (t = 1), `j = 3 + v2(K−1)`.

**Children.**
* The deletion children are `h_D = (n+1)/2^D − 1`, with run ends `y_D = (x+1)/3^D − 1`. At the run ends the pair chain of THM-4581 is in state `(D, 3^D − 1)`, with Terras shift D: anchored at −1. The reset child is D = 1.
* Generalized deletions `3^a(n+1)/2^b − 1` (`3^a < 2^b`) give (b−a)-chains. These are exactly the K-uniform collision certificates on unrefined 2-adic classes.
* 3-adic refinements give others: for 3 | t, `(2n−1)/3 < n` maps to n in one step.

**Lag-one form (THM-4601 (i)).**
* `x − 1 = 3z + 1`, where `z = 2·3^(K−2) t − 1` is the last point of the reset child's ones-run.
* The reset rule is the translation pair (x, x−1) merging at once when x ≡ 5 mod 8.
* The residual branch is the same pair with x ≡ 1 mod 8.

## 2. Anchors (THM-4600)

**The lemma.**
* For a letter a put `c_a = 1/(2^a − 3)` (c_1 = −1, c_2 = 1, c_3 = 1/5, …). Then `3c_a + 1 = 2^a c_a`.
* If odd u, v satisfy `u − c_a = 3^k (v − c_a)`:
  * v has U-letter a ⟺ `v2(v − c_a) ≥ a + 1` ⟺ u has U-letter a;
  * in that case `U(u) − c_a = 3^k (U(v) − c_a)`.
* Hence common runs of a of any length are transparent, and both runs end at the same step.

**General form.**
* For a periodic word w read at its rational cycle point `c_w = B_w/(2^(Σw) − 3^|w|)`, anchored states are carried around the cycle as long as `v2(v − c_w) > Σw`.
* The condition is sharp. From another cycle point, read the rotated word.

**Borel–torus reading.**
* A pair-chain state is the affine map `A(x) = 3^k x + e` with u = A(v), and a common step τ acts by conjugation.
* The states fixed under conjugation by the run map `τ_a(x) = c_a + (3/2^a)(x − c_a)` are exactly its centralizer: the dilations about c_a, a torus.
* The translation coordinate `ε = e − c_a(1 − 3^k)` is multiplied by `3/2^a` at each common letter:
  * archimedean-expanding for a = 1, contracting for a ≥ 2;
  * 2-adically expanding by `2^a`, so non-anchored states leave a run after about `v2(ε)/a` letters.
* This is the affine-group analogue (ANALOGY) of THM-4553's "PSL(2,7) > Borel > torus" ladder. The finite-field ladder itself follows from THM-4602 (3)(b) and Babai–Cameron's Lemma 3.1.

## 3. The +1 barrier (THM-4601 (ii), HYP-9240)

**The limit chain.**
* As j → ∞ the run end x tends 2-adically to 1, and `y_D → y_D* = (2 − 3^D)/3^D`.
* For s < j the actual parities agree with those of 1 and y_D*. So the actual chain *is* the limit chain up to time j.
* The bound is exact: the states always differ at j + 1 (8947 chains, audit A).

**Exact computation of the limit chain** (twoanchor_core.py for D ≤ 3000; barrier_extend.py for 3001..8000; reproduced by audit A):
* D = 1, 2: the child's limit reaches 0, at Terras time 1 and 3 respectively. The debt then grows by one every two steps.
* D ≥ 3:
  * the child's limit makes exactly D denominator-clearing odd steps;
  * it lands on a positive integer `N_D` (≤ 880 for D ≤ 3000, ≤ 2527 for D ≤ 8000);
  * it enters the trivial cycle with debt `k_∞(D) ≥ 3`.
* In phase / out of phase: 2288 / 710 for D ≤ 3000, and 3610 / 1390 for 3001..8000.
* NUMERICAL: in-phase `k_∞/D ∈ [0.75, 1.47]`, and `[0.966, 1.045]` for D ≥ 2000.
* **No D ≤ 8000 absorbs.**

**Consequence.**
* No generalized deletion certificate (D ≤ 8000) completes by Terras time j after x, in particular not inside the two-run.
* Every K-uniform collision certificate on unrefined 2-adic classes therefore has depth beyond the two-run length. The explicit heads below complete at j + 11 and j + 9, the earliest possible.
* This is the positive-cycle twin of the sign barrier of THM-4594.
* 3-adic refinements are not bound by it. For example `(2n−1)/3` when 3 | t.

**HYP-9240 for all D.**
* Absorption would need `odd(N_D) = ⌈(s_0 + σ_T(N_D))/2⌉`, i.e. an odd-step excess of about D in the orbit of the small landing N_D.
* For N ≤ 880 the excess is at most 9.

## 4. Universal state and ladder compiler (THM-4601 (iii)–(iv))

**Clearing identity.**

    F_(2,2,2)(x) − 1 = 27 (F_(4,1,1)(y_3) − 1),   y_3 = (x+1)/27 − 1,

* The child's side equals `(x+63)/64`.
* The child's letters (4,1,1) are actual iff x ≡ 1 mod 128, i.e. J ≥ 3.
* The reset partner y_4 reaches the same value with letters (2,2,1,1).
* THM-4600 with a = 2 then carries `X − 1 = 27(Y − 1)` through the remaining J − 3 letters 2.
* Verified on 840 sources with K = 4..1500 and J = 3..150 (audit A).

**Ladder compiler.** Take heads u, v with `|v| = |u| + 3` and either
* a child ladder: `Σu = Σv + 2i` and `F_v(1) = 4^i F_u(1) + (4^i − 1)/3`, or
* a source ladder: `Σv = Σu + 2i` and `F_u(1) = 4^i F_v(1) + (4^i − 1)/3`.

Then for all K ≥ 4, J ≥ 3, c ≥ 1, as affine maps (child ladder),

    F_(1^(K−4) 4 1 1 2^(J−3) v (c+2i)) ((n+1)/8 − 1)  =  F_(1^(K−1) 2^J u c)(n).

* The final letters come from the 4z+1 ladder: `z ↦ 4z + 1` preserves U and adds 2 to the letter.
* **Why the child's word is actual.** The common endpoint is an odd integer. That forces every intermediate child value to be an odd integer with exactly the stated valuations.
* This is the compiler of [reset-two rules](collatz_reset2_rules_20261007.md) §2 (anchor −1, i = 1), moved to the anchor +1 and composed with the clearing identity.

**Explicit heads** (i = 1; `F_u(1)` → `F_v(1)`):
* r = 1: `(1,10) ~ (1,1,2,2,3)`, 7/1024 → 263/256. Post-run Terras depth 12 (time j + 11 after x).
* r = 2: `(10) ~ (3,1,1,3)`, 1/256 → 65/64. Post-run Terras depth 11 (time j + 9).
* These are the earliest possible absorptions.
* They certify every one of 560 random residual sources tested (K ≤ 1001, J ≤ 150) in their native classes: 8 (r = 1) and 32 (r = 2) of the 8192 odd classes of w′ mod 2^14.

**Head grammar** (heads_plus1.py, FINITE-EXACT).
* The child-ladder i = 1 pairs with `Σu ≤ 21` number 6237. The first are:
  * `(10)~(3,1,1,3)`;
  * `(1,10)~(1,1,2,2,3)`;
  * `(9,2)~(3,1,1,3,1)`;
  * `(1,9,2)~(1,1,2,2,3,1)`;
  * `(8,4)~(3,1,1,4,1)`;
  * `(9,1,2)~(3,1,1,3,1,1)`.
* No source head is a prefix of another. A pair's source class has measure `2^(r−Σu)` within its run type.
* Union: 0.0139 (r = 1) and 0.0295 (r = 2).
* The i = 1 grammar is not complete: child ladders i ≥ 2 (first `(1,2,9)~(1,1,1,1,1,1)`, i = 3) and source ladders (first `(12,1,1)~(3,1,1,3,2,6)`) also occur.
* Exact absorbed mass by post-run depth 22: 0.0169 (r = 1) and 0.0326 (r = 2). The i = 1 grammar carries 82% and 90% of it (audit A).
* This is the +1 counterpart of the complete −1 parser in [collatz_collision_dp_20261007.md](collatz_collision_dp_20261007.md).

## 5. The variable-depth rule (THM-4601 (v)) and its numbers

**The rule.**
* For every D, the D-chain state at the end of the two-run is the limit chain's state at Terras time 2J.
* It depends only on (D, J), is the same for r = 1 and 2, and is constant for `J ≥ J_0(D)`. J_0 = 3, 3, 6, 6, 8, 8, 10, 10, 15, 15 for D = 3..12.
* r enters only through the child's forced post-run prefix: (1, 1) if r = 1, (1, 0, 0) if r = 2.
* The certificate is the first absorption of one of these fixed chains, driven by the post-run bits. It is uniform in K, t and J ≥ J_0(D).
* Every such chain absorbs Haar-almost surely (THM-4581).

**End-of-run states** (short_types.py, audit A):

| D | J = 1 | J = 2 | J = 3 | J = 4 | J = 5 | from J_0 on |
|---|---|---|---|---|---|---|
| 1, 2 | (D, (3^D−1)/2) | (2, 1) | (3, 1) | (4, 1) | (5, 1) | never constant: (J, 1) |
| 3, 4 | (D, (3^D−1)/2) | (4, 10) | (3, −26) | (3, −26) | (3, −26) | **(3, −26)**, J_0 = 3 |
| 5, 6 | (D, (3^D−1)/2) | (6, 91) | (6, −296) | (6, −404) | (6, −485) | (6, −728), J_0 = 6 |

* For J = 1 the states are anchored at −1/2, the fixed point of the unhalved map 3x+1. For J = 2 they are anchored at −1/8.

**Haar coverage within post-run depth 400** (child_compare.py, 300 sources per row, 400-bit t):

| J | D = 1 (reset child) | D = 3 (= D = 4) | any D ≤ 8 |
|---|---|---|---|
| 1 | 0.45 / 0.55 | 0.35 / 0.32 | 0.67 / 0.72 |
| 3 | 0.35 / 0.38 | 0.34 / 0.38 | 0.56 / 0.53 |
| 10 | 0.16 / 0.20 | 0.37 / 0.42 | 0.59 / 0.58 |
| 16 | 0.08 / 0.10 | 0.37 / 0.41 | 0.57 / 0.56 |
| 32 | 0.003 / 0.007 | 0.35 / 0.40 | 0.52 / 0.56 |

(each cell: r = 1 / r = 2)

* The reset child collapses like a random walk started at height J. h_3 and h_4 are flat in J.
* Weighted by the Haar mass of the run types, the deletion children D ≤ 8 certify 0.646 ± 0.014 of the whole residual branch within 400 post-run steps (audit A).

**From the universal state with realisable post-run bits** (audit A):

| s | 30 | 100 | 200 | 300 | 400 | 1000 |
|---|---|---|---|---|---|---|
| ρ_1(s) | 0.038 | 0.172 | 0.269 | 0.328 | 0.375 | 0.519 |
| ρ_2(s) | 0.053 | 0.184 | 0.278 | 0.338 | 0.384 | 0.529 |

* The session's first numbers used unconstrained fair bits from (3, −26): 0.18, 0.53, 0.72 at s = 100, 1000, 4000. The shortest such pattern has length 9, but it starts with an even bit and cannot follow a two-run.

## 6. Mersenne line

**Shells.** For odd K with `m = v2(K−1) ≥ 4`:
* `J = 1 + ⌊m/2⌋ ≥ 3` and `r = 1 + (m mod 2)`.
* The post-run class `w′ = 3^(J−3)(3^(K−1) − 1)/2^(m+2) mod 2^L` runs bijectively over the odd classes as `(K−1)/2^m` does.
* So the density of K certified by K → K−3 within post-run depth s is the same `ρ_r(s)` in every shell.

**Observed** (K < 12000, depth 200): m = 4..8 give 0.283, 0.283, 0.255, 0.319, 0.304.

**Two exponents, one head pair** (verified by direct orbits; common odd endpoints of 2983 and 8310 bits):
* `2^1889 − 1 ⇝ 2^1886 − 1`: m = 5, J = 3, words `1^1888 2^3 (10) 3` and `1^1885 (4,1,1) (3,1,1,3) 5`.
* `2^5249 − 1 ⇝ 2^5246 − 1`: m = 7, J = 4, words `1^5248 2^4 (10) 1` and `1^5245 (4,1,1) 2 (3,1,1,3) 3`.
* The r = 1 pair occurs too: `2^6129 − 1 ⇝ 2^6126 − 1`, m = 4, J = 3, words `1^6128 2^3 (1,10) 1` and `1^6125 (4,1,1) (1,1,2,2,3) 3` (endpoint of 9706 bits).

**Short shells.** m = 1, 2, 3 have J = 1, 2, 2. 2^1459 − 1 has m = 1. The end states for these J are the fixed states of the table in §5.

## 7. Limits (what the rule does not do)

* **Heavy tail.** Under Haar, coverage tends to 1 but slowly.
* **Positive integers.** A finite pair whose orbits both reach 1 with unequal odd-step counts never merges at equal time. So a single D certifies only about 35–42% of large random residual sources within 400 post-run steps, and about 28% of Mersenne shells within 200.
* **Orphans.** THM-4556 (vi): 2% of the odd exponents (50 of the 2500 odd a in [10^3, 6·10^3]) have no equal-time Mersenne partner. For them the deletion rule fails, and other certificates are needed: 3-adic ones, non-Mersenne children, or descent of depth ∝ K.
* **Re-anchoring is possible in both directions.** The source's own limit chain at x* = 1 is positive and never reaches a negative anchor. After the two-run, though, a long post-run ones-run of the child (`Y ≡ −1 mod 2^L`, r = 1) takes (3, 1−27) to `(−5, 3^(−5) − 1)`, anchored at −1, in 12 Terras steps via the limit pair (−53, −1). The rest of that ones-run is then transparent.
  * Corrected after audit A: the session first claimed that no +1 → −1 transition exists.

## 8. Cluster algebras, total positivity, friezes and the tournament work

**8.1 Borel and torus** (section 2).
* Every state with k ≠ 0 is the dilation by 3^k about its fixed point `c = −e/(3^k − 1)`. So the Borel group minus translations is the union of the tori `T_c`, which meet only at the identity.
* On a common step τ the new state `τAτ^(−1)` has fixed point τ(c). **The anchor moves by the same affine step as the two orbits**: it is a virtual third orbit, the point where u and v would coincide.
* Run transparency is the case where the anchor sits at the fixed point of the run's letter. Debt-changing steps change k by ±1 and make the anchor jump.
* Cluster-algebra reading (ANALOGY only):
  * seeds correspond to anchors;
  * cluster tori correspond to the tori T_c;
  * mutations correspond to debt-changing steps.

  The clearing identity `(2,2,2) ~ (4,1,1)` is the analogue of a mutation path from the −1 torus to the +1 torus.

**8.2 Friezes and Collatz.**
* The tiling package showed that the prefix carries of a word give a positive configuration in Gr(2, r+1), whose Plücker coordinates are subword carries.
* Here, the head condition is the reset-two compiler moved from the anchor −1 to +1, with the 4z+1 ladder.
* An honest frieze identity for these heads was not found; the tiling note's ear-move obstruction stands. Type: ANALOGY at most.

**8.3 Tournaments as chirotopes (THM-4602; statements largely KNOWN).**
* A tournament is a skew sign matrix. It is realisable by vectors in R^2 iff all 4-Pfaffians `b_ij b_kl − b_ik b_jl + b_il b_jk` are ±1, iff it is locally transitive (Babai–Cameron 2000).
  * `|Pf| = 3` exactly on the 4-sets with one 3-cycle.
  * Labelled counts are (n−1)!·2^(n−1) for n = 4..9.
* The transitive tournament is the totally positive part. Positive real SL2 friezes are its slice `p_(i,i+1) = p_(1n) = 1`, and Conway–Coxeter friezes are the integer points of that slice.
* Over F_p, the standard chart of P^1(F_p) has Legendre pattern Paley + sink. Scaling a vector by a non-residue switches its vertex, so the switching class is an invariant (KNOWN: Paley 1933; Gunderson–Semeraro).
  * Any ordering of the p+1 points gives the class.
  * New here: the PGL2 Singer ordering closes as a frieze of width p−2 with 2-periodic quiddity.
* Constant quiddities: elliptic ones cover at most (p+1)/2 points, hyperbolic ones at most (p−1)/2. The parabolic quiddity −2 gives a closed frieze through p points with the QR_p class. This corrects the session's first attempt.
* ANALOGY: the real chart gives the transitive order, and the F_p chart gives Paley. Total positivity has a different counterpart: χ-positive configurations are transitive subtournaments of Paley + sink, with at most tt(QR_p) + 1 points (attained for p ≤ 43; logarithmic growth observed but unproved).

**8.4 Catalan coincidence.**
* THM-438's leading coefficient C_k is a signed Möbius / free-cumulant sum over even-series patterns, not a count of plane-tree tours (MISTAKE-060/061, THM-438 ADDENDUM-2).
* Conway–Coxeter friezes with k+2 points are also counted by C_k. No bijection is known.
* Type: NUMERICAL COINCIDENCE. The session first called this a dictionary from THM-438's superseded headline (MISTAKE-586).

**8.5 Pair independence (HYP-9241).**
* q_R(T) is the probability that y, y+1, …, y+R have pairwise distinct Terras orbits up to time T (Haar y).
* `q_2/q_1^3` = 1.11–1.19 and `q_3/q_1^6` = 1.7–1.9 for T = 16..512, while q_3 falls by a factor 35. The value 1.30 at T = 1024 is a fluctuation. Local-slope ratios are about 1 : 2.9 : 5.6.
* Audit B: `q_2 = q_01 q_12 q_02` within 3% for T ≥ 64, so `c_2 ≈ q_02/q_01 ≈ 1.2` is a lag effect.
* So the R(R+1)/4 exponents are the vicious-walker exponents, but the constant is that of **independent pairs**. Karlin–McGregor repulsion would give π/4 for Brownian walkers.
* No determinantal or total-positivity structure is exhibited. The session's first "Karlin–McGregor law" is withdrawn (MISTAKE-586).

**8.6 Rejected or not pursued.**
* Real-rootedness of the tournament F-polynomial is not universal (about 89% in the Worpitzky synthesis), and F depends on the labelling.
* Quiver mutation does not preserve tournaments in general. At a source or sink it is switching, which is the symmetry of (3)(b).
* The frieze ear move is not a Collatz rewrite (tiling note).

## 9. The four attached papers

* **Primitive sextic torus packets** (equidistribution in SL6(Z)\SL6(R), §§6–7). Type: ANALOGY; it shaped section 2.
  * Long diagonal segments correspond to long runs.
  * The cross-root group Z_J (normal in VA, h = zd) corresponds to the translation coordinate ε relative to the anchor.
  * An inactive cut corresponds to an anchored state, the only kind that survives arbitrarily long runs.
  * Their arithmetic tube bound `e^(−ηt)` corresponds to the 2-adic loss of precision `2^(−a)` per letter.
* **Exact orders and invariant anticanonical systems** ("why finite generation attains the limit"). Type: ANALOGY.
  * The certificate family is finitely generated over the free run parameters (K, J): finitely many heads per post-run depth, with transparent runs.
  * Exact attainment for every integer, the analogue of their theorem, is precisely what remains open.
* **Kakeya maximal (3D) and Kakeya (4D).** Multiscale tube incidence. No usable transfer was found, so they were not used.

## 10. Next obligations

1. HYP-9240 for all D: bound the landing `N_D`, or its odd-step excess against `s_0/2`.
2. Compile the short-J states (J = 1, 2) and the universal states into complete ladder grammars (child and source ladders, all i). Then give the exact Haar coverage of the whole residual branch per depth.
3. A positive-integer statement for the universal state: the probability of absorption before both orbits reach 1, as a function of size.
4. Prove pair independence (exponent 3/2) for R = 2 from THM-4581's weight.
5. Anchors at −5 and −17 (contracting inverse cycles) and transitions into them.
   * The concurrent [child-closure synthesis](collatz_child_closure_synthesis_20261007.md) uses the generators `G_1(n) = (2n−1)/3` and `G_5(n) = (8n−5)/9`.
   * These are the contracting inverse branches of the periodic anchors −1 (word 1) and −5 (word 12) in THM-4600 (2). Its mixed words are walks between negative anchors.

## 11. Reproduction

All scripts are in `04-computation/experiments/twoanchor_20261007/`:
* `twoanchor_core.py [Dmax]`: parts A–E; ALL CHECKS PASSED (119,914 assertions, about 9 s).
* `barrier_extend.py 8000 3001` (+ .out): the barrier for 3001 ≤ D ≤ 8000.
* `short_types.py`, `heads_plus1.py` (+ .out): end states and the i = 1 head grammar.
* `tournament_chirotope.py`: the 4-Pfaffian criterion and Legendre frieze strips; ALL CHECKS PASSED.
* `child_compare.py`, `mersenne_shells.py` (+ .out): the J-table and the Mersenne shells.
* `km_exponent.py` (+ `_n3000.out`, `_n20000.out`) and `km_analysis.py` (+ .out): the partner law.
* `audit_A/`, `audit_B/`: independent audit scripts, outputs and reports.
