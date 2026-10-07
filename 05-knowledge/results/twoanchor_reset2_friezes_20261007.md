# Two anchors for the residual first-reset-2 branch; chirotopes, friezes and a Karlin–McGregor law

**Session.** mac-mini-2026-10-07-twoanchor.

**Owner prompt.** "work towards a variable-depth rule for the residual first-reset-2 branch, think cluster algebras, total positivity and how friezes and conway-coxeter connect with our tournament adjacent work, again take inspiration from the attached papers but go much further beyond them". The four attached OpenAI preprints are treated as trusted; section 9 says what was used from them.

**Starting point.** Section 4 of [terminal lifts](collatz_terminal_lifts_20261007.md) leaves the residual branch open. Concurrent work on it:
* [reset-two rules](collatz_reset2_rules_20261007.md): guarded K → K−2 and K → K−6;
* [reset-two synthesis](collatz_reset2_synthesis_20261007.md): ternary turns and a ledger;
* [frieze package](frieze_collatz_reset2_20261007.md): positive marked minors, and the ear-move obstruction.

This note is complementary. The aim is uniformity in the second unbounded parameter of the branch, the length of its two-run, which no earlier rule controls.

COLLATZ IS STILL OPEN. Nothing here covers every positive integer.

## 0. Results and types

| Item | Statement | Type |
|---|---|---|
| [THM-4600](../../01-canon/theorems/THM-4600-anchored-debts-are-run-transparent-the-borel-torus-picture-of-pair-chains.md) | Pair-chain states anchored at the letter-a fixed point `c_a = 1/(2^a−3)` pass unchanged through common runs of a. They are the affine maps commuting with the run map: a torus in the Borel group. `a = 1` is THM-4555's collision at −1; `a = 2` (the trivial cycle +1) is new. | PROVED |
| [THM-4601](../../01-canon/theorems/THM-4601-two-anchor-reduction-of-the-residual-first-reset-two-branch.md) (i) | The residual branch is exactly the lag-one translation pair (x, x−1) with x ≡ 1 mod 8. | PROVED |
| THM-4601 (ii) | **+1 barrier.** For no deletion child h_D with D ≤ 3000 does the class-decided collision with a residual source complete before the end of the source's two-run. | PROVED for D ≤ 3000 (exact limit chains); HYP-9240 for all D |
| THM-4601 (iii) | **Universal residual state.** If the two-run has J ≥ 3 letters, the K → K−3 and K → K−4 chains end the run exactly in `(3, 1−27)`, i.e. `X − 1 = 27(Y − 1)`, whatever K, t and J are. | PROVED |
| THM-4601 (iv) | **Two-anchor compiler.** Head collisions at +1 (`|v|−|u| = 3`, `Σu = Σv+2`, `F_v(1) = 4F_u(1)+1`) give certificates uniform in both run lengths K and J. Explicit pairs `(1,10)~(1,1,2,2,3)` and `(10)~(3,1,1,3)`. | PROVED (an identity of affine maps; the child's word is forced) |
| THM-4601 (v) | **The variable-depth rule.** Six run types (j = 3, 4, 5, 6 and j ≥ 7 with r = 1, 2). For each type and each D the problem is one fixed pair-chain problem in the post-run bits, and its Haar coverage tends to 1 (THM-4581). | PROVED |
| THM-4601 (vi) | Mersenne line: in every shell `v2(K−1) = m ≥ 4` the density certified by K → K−3 within post-run depth s is the same number `ρ_r(s)`. `2^1889−1 ⇝ 2^1886−1` and `2^5249−1 ⇝ 2^5246−1` use one and the same head pair. | PROVED (shell law), FINITE-EXACT (examples) |
| [THM-4602](../../01-canon/theorems/THM-4602-tournaments-as-chirotopes-pfaffian-criterion-and-legendre-frieze-strips.md) | A tournament is a sign pattern of Gr(2,n) iff every 4-Pfaffian is ±1 iff it is locally transitive (`|Pf| = 3` exactly on 4-sets with one 3-cycle). The transitive tournament is the positive (Conway–Coxeter) part. Over F_p (p ≡ 3 mod 4) the Legendre pattern of P^1(F_p) is Paley + sink up to switching, and every SL2 frieze strip over F_p gives a switched induced subtournament. | PROVED; ingredients KNOWN |
| [HYP-9241](../hypotheses/HYP-9241-karlin-mcgregor-product-law-for-collatz-translation-partners.md) | Karlin–McGregor law: `q_R(T) ≍ q_1(T)^(R(R+1)/2)` for translation partners y, …, y+R. | NUMERICAL |

## 1. The branch, exactly

**Setup.**
* `n = 2^K t − 1`, t odd, K ≥ 2.
* After `K − 1` letters 1 the orbit is at `x = 2·3^(K−1) t − 1`. Write `x = 1 + 2^j u` with u odd.
* First reset 2 means `v2(3^K t − 1) = 1`. That is equivalent to x ≡ 1 mod 8, i.e. j ≥ 3.
* The source then has `J = ⌊(j−1)/2⌋` letters 2. After them comes letter 1 if j is odd (r = 1), or a letter ≥ 3 if j is even (r = 2), where `r = j − 2J`.
* On the Mersenne line (t = 1), `j = 3 + v2(K−1)`. This matches the tiling diagnostic.

**Children.** The deletion children are `h_D = (n+1)/2^D − 1`, with run ends `y_D = (x+1)/3^D − 1`.
* At the run ends the pair chain of THM-4581 is in state `(D, 3^D − 1)`, with Terras shift D: anchored at −1.
* The reset child is D = 1.

**Lag-one form (THM-4601 (i)).**
* `x − 1 = 3z + 1`, where `z = 2·3^(K−2) t − 1` is the last point of the reset child's ones-run. So x − 1 lies on the reset child's orbit, one Terras step before `(x−1)/2 = T^(K−1)(h_1)`.
* So the reset rule is the translation pair (x, x−1) merging at once (x ≡ 5 mod 8, `U(4u+1) = U(u)`).
* The residual branch is that same lag-one pair with x ≡ 1 mod 8. There the j forced halvings of x − 1 build the debt.

## 2. Anchors (THM-4600)

**The lemma.**
* For a letter a put `c_a = 1/(2^a − 3)` (c_1 = −1, c_2 = 1, c_3 = 1/5, …). Then `3c_a + 1 = 2^a c_a`.
* If odd u, v satisfy `u − c_a = 3^k (v − c_a)`, then:
  * v has U-letter a ⟺ `v2(v − c_a) ≥ a + 1` ⟺ u has U-letter a;
  * and then `U(u) − c_a = 3^k (U(v) − c_a)`.
* Hence common runs of a of any length are transparent for anchored states, and both runs end at the same step.

**General form.** For a periodic word w with rational cycle point `c_w = B_w/(2^(Σw) − 3^|w|)`, anchored states are carried around the cycle as long as `v2(v − c) > Σw`.

**Borel–torus reading.**
* A pair-chain state is the affine map `A(x) = 3^k x + e` with u = A(v). A common step τ acts by conjugation.
* The states fixed under conjugation by the run map `τ_a(x) = c_a + (3/2^a)(x − c_a)` form exactly its centralizer: the dilations about c_a, a torus.
* Relative to the anchor, the unipotent coordinate `ε = e − c_a(1 − 3^k)` is multiplied by `3/2^a` at each common letter:
  * archimedean-expanding for a = 1 (near −1), contracting for a ≥ 2 (near +1);
  * always 2-adically expanding, by `2^a` per letter. That is why non-anchored states leave a run after about `v2(ε)/a` letters.
* This is the "PSL(2,7) > Borel > torus" structure of THM-4552/4553, promoted to a mechanism.

## 3. The +1 barrier (THM-4601 (ii), HYP-9240)

**The limit chain.** As j → ∞ the run end x tends 2-adically to 1, and `y_D → y_D* = (2 − 3^D)/3^D`. For s < j the actual parities of x's and y_D's orbits agree with those of 1 and y_D*. So the actual chain *is* the limit chain up to time j.

**Exact computation of the limit chain, D ≤ 3000** (twoanchor_core.py, part B):
* D = 1, 2: the child's limit goes to 0, so the debt grows by one every two steps.
* D ≥ 3:
  * the child's limit makes exactly D denominator-clearing odd steps;
  * it lands on a positive integer `N_D ≤ 880` and enters the trivial cycle with debt `k_∞(D) ≥ 3`, where k_∞/D → 1;
  * 2288 values of D are in phase, 710 out of phase.
* **No D absorbs.**

**Consequence.** No D-chain is absorbed, i.e. no class-decided collision with h_D completes, before Terras time j + 1 after the run end, i.e. before the end of the two-run.
* The two-run is never skipped. The certificate depth of any K-uniform rule is at least 2J; the rule below attains 2J + O(1).
* This is the positive-cycle twin of the sign barrier of THM-4594.

**Heuristic for all D.**
* A landing debt of about D/2 or more would need an odd-step excess of that size in the orbit of a small positive integer `N_D`.
* `N_D ≤ 880` here, with no growth trend.

## 4. Universal state and two-anchor compiler (THM-4601 (iii)–(iv))

**Clearing identity.**

    F_(2,2,2)(x) − 1 = 27 (F_(4,1,1)(y_3) − 1),   y_3 = (x+1)/27 − 1,

* This is an identity of affine maps; its value is `(x+63)/64` for the child.
* The child's letters (4,1,1) are actual iff x ≡ 1 mod 128, i.e. J ≥ 3.
* The reset partner y_4 reaches the same value with letters (2,2,1,1).
* Then THM-4600 with a = 2 carries `X − 1 = 27(Y − 1)` through the remaining J − 3 letters 2.
* Verified on 640 sources with j from 16 to 64; similarly `(6, 1−3^6)` for D = 5, 6.

**Compiler.** Take heads u, v with

    |v| = |u| + 3,   Σu = Σv + 2,   F_v(1) = 4 F_u(1) + 1.

Then for all K ≥ 4, J ≥ 3, c ≥ 1, as affine maps,

    F_(1^(K−4) 4 1 1 2^(J−3) v (c+2)) ((n+1)/8 − 1)  =  F_(1^(K−1) 2^J u c)(n).

**Why the child's word is actual.** If the source word is actual (a congruence on n), the common endpoint is an odd integer. That forces every intermediate value of the child to be an odd integer with exactly the stated valuations:
* a value with negative 2-adic valuation never becomes integral again;
* an even integer value would make the next one non-integral.

This is the at-+1 version of the compiler in [reset-two rules](collatz_reset2_rules_20261007.md) (§2), composed with the clearing identity.

**Explicit heads** (`F_u(1)` → `F_v(1)`):
* r = 1: `(1,10) ~ (1,1,2,2,3)`, 7/1024 → 263/256. Post-run depth 11.
* r = 2: `(10) ~ (3,1,1,3)`, 1/256 → 65/64. Post-run depth 9.

Each certified every one of 240 random residual sources tested, with K up to 1001 and J up to 150, in its native class. The class count is 8 (r = 1) or 16 (r = 2) of 8192 odd classes of w′ mod 2^14.

The shortest absorbing patterns from the universal state have lengths 9, 10, 11, … with counts 1, 1, 3, 7, 15, 30, 67, 147, 301, 658, 1357, 2783 (lengths 9–20).

## 5. The variable-depth rule (THM-4601 (v)) and its numbers

**The rule.**
* Classify a residual source by its run type: j ∈ {3, 4, 5, 6}, or j ≥ 7 with r ∈ {1, 2}.
* For every D, the D-chain state at the end of the two-run depends only on (type, D):
  * for D ≤ 2 and j ≥ 7 it also depends on J, with debt about J;
  * for D ∈ {3, 4} and j ≥ 7 it is the universal `(3, 1−27)`.
* The certificate is then the first absorption of one of these fixed chains, driven by the post-run bits. It is uniform in K and t, and in J for D ≥ 3.
* Every such chain absorbs Haar-almost surely (THM-4581).

**Haar coverage within post-run depth 400** (child_compare.py, 300 sources per row, 400-bit t):

| J | D = 1 (reset child) | D = 3 (= D = 4) | any D ≤ 8 |
|---|---|---|---|
| 1 | 0.45 / 0.55 | 0.35 / 0.32 | 0.67 / 0.72 |
| 3 | 0.35 / 0.38 | 0.34 / 0.38 | 0.56 / 0.53 |
| 10 | 0.16 / 0.20 | 0.37 / 0.42 | 0.59 / 0.58 |
| 16 | 0.08 / 0.10 | 0.37 / 0.41 | 0.57 / 0.56 |
| 32 | 0.003 / 0.007 | 0.35 / 0.40 | 0.52 / 0.56 |

(each cell: r = 1 / r = 2)

* The reset child collapses like a random walk started at height about J.
* h_3 and h_4 are flat in J, as THM-4601 (iii) predicts.
* From the universal state (Monte Carlo, 20,000 fair-bit runs), P(absorbed by s) is 0.054 (s = 30), 0.18 (100), 0.34 (300), 0.53 (1000), 0.72 (4000). Its partner state `(3, 2−54)` behaves the same.

## 6. Mersenne line

**Shells.** For odd K with `m = v2(K−1) ≥ 4`:
* `J = 1 + ⌊m/2⌋` and `r = 1 + (m mod 2)`.
* The post-run class `w′ = 3^(J−3)(3^(K−1) − 1)/2^(m+2) mod 2^L` runs bijectively over the odd classes as `(K−1)/2^m` does.
* So the density of K certified by K → K−3 within post-run depth s is the same `ρ_r(s)` in every shell.

**Observed** (K < 12000, depth 200): m = 4..8 give 0.283, 0.283, 0.255, 0.319, 0.304. Pooled: 0.277 (r = 1) and 0.292 (r = 2).

**Two exponents, one head pair.**
* `2^1889 − 1 ⇝ 2^1886 − 1`: m = 5, J = 3, words `1^1888 2^3 (10) 3` and `1^1885 (4,1,1) (3,1,1,3) 5`.
* `2^5249 − 1 ⇝ 2^5246 − 1`: m = 7, J = 4, words `1^5248 2^4 (10) 1` and `1^5245 (4,1,1) 2 (3,1,1,3) 3`.

**Coverage of all shells.**
* Shells m = 1, 2, 3 (J = 1, 2, 2; 2^1459 − 1 has m = 1) are the three fixed short types.
* Together with the two universal types, the whole Mersenne residual branch is five fixed chain problems.

## 7. Limits (what the rule does not do)

* **Heavy tail.** Under Haar, coverage tends to 1 but slowly (0.72 at depth 4000 from the universal state).
* **Positive integers.** A finite pair whose orbits both reach 1 with unequal odd-step counts never merges at equal time. So a single D certifies only about 35–42% of large random residual sources within 400 post-run steps, and about 28% of Mersenne shells within 200.
* **Orphans.** Some sources have no deletion partner at all (THM-4556 (vi): 50 of 2500 Mersenne exponents in [10^3, 6·10^3]). For them no K-uniform certificate exists, and a complete rule must fall back to certificates of depth ∝ K (descent).
* **No re-anchoring after the two-run.** Positive limit points never reach negative anchors (T preserves positivity), so the chain cannot pass from anchor +1 to anchor −1. A long ones-run after the two-run therefore raises the debt again. Only negative → positive transitions occur; −1 → +1 is the clearing identity.

## 8. Cluster algebras, total positivity, friezes and the tournament work

**8.1 Borel and torus** (section 2). Pair chains are conjugation dynamics in the Borel subgroup. Anchored states are tori and runs are diagonal flows. This is the concrete mechanism behind "Borel > torus".

**8.2 Friezes and Collatz.**
* The tiling package showed that the prefix carries of a word give a positive configuration in Gr(2, r+1). Its Plücker coordinates are subword carries (in-repo).
* Here: the head condition `F_v(1) = 4F_u(1) + 1` is the compiler of the reset-two rules moved from the anchor −1 to +1. The clearing identity is a transition between the two anchors.
* An honest frieze identity of these heads was not found (the ear-move obstruction of the tiling note stands). Type: ANALOGY at most.

**8.3 Tournaments as chirotopes (THM-4602).**
* A tournament is a skew sign matrix. It comes from vectors in R^2 iff it is locally transitive, iff all 4-Pfaffians `b_ij b_kl − b_ik b_jl + b_il b_jk` are ±1.
  * `|Pf| = 3` exactly on the 4-sets with one 3-cycle.
  * The transitive tournament is the totally positive part, which is where positive friezes live (Conway–Coxeter).
* Over F_p, with "positive" read as "quadratic residue", the configuration of all points of P^1(F_p) gives the Paley tournament plus a sink, up to switching. Scaling a vector by a non-residue switches its vertex.
* An SL2 frieze strip over F_p gives a switched induced subtournament of it. Singer strips through all p+1 points realise the whole class (checked for 13 primes up to 83).
* So the finite-field analogue of "positive frieze ↔ transitive tournament" is "frieze over F_p ↔ Paley tournament": the Paley tournament is the finite-field totally positive tournament.
* Correction made during the session: constant-quiddity (elliptic Chebyshev) friezes cover only (p+1)/2 points, since the elliptic torus of SL2 has image of order (p+1)/2 in PSL2. Covering all p+1 points needs the PGL2 Singer cycle.

**8.4 Catalan dictionary.**
* THM-438's leading Paley cluster-integral patterns are Euler tours of plane trees, counted by C_k.
* The same Catalan numbers count triangulations, hence Conway–Coxeter friezes, via the standard bijections. Type: DICTIONARY (KNOWN bijections), not a new identity.

**8.5 Karlin–McGregor (HYP-9241).**
* Let q_R(T) be the probability that y, y+1, …, y+R have pairwise distinct Terras orbits up to time T (Haar y).
* Measured over 20,000 samples: `q_2/q_1^3` = 1.11–1.30 and `q_3/q_1^6` = 1.7–1.9 for T = 16..512, while q_3 falls by a factor 35.
* Local slopes are in ratio about 1 : 2.9 : 5.6, against the vicious-walker 1 : 3 : 6.
* With THM-4593's α = 1/2 this predicts exponents `R(R+1)/4`: total-positivity (Karlin–McGregor) structure in Collatz coalescence. NUMERICAL; the asymptotic regime has not been reached at T ≤ 2048 (q_1's local slope is still 0.34).

**8.6 Rejected or not pursued.**
* Real-rootedness of the tournament F-polynomial is not universal (the repo's Worpitzky synthesis has about 89%, and F depends on the labelling). It gives no TP law.
* Quiver mutation does not preserve tournaments (2-cycles cancel, double arrows appear).
* The frieze ear move is not a Collatz rewrite (tiling note).

## 9. The four attached papers

* **Primitive sextic torus packets** (equidistribution in SL6(Z)\SL6(R), §§6–7):
  * long diagonal segments correspond to long runs;
  * the cross-root group Z_J (normal in VA, h = zd) corresponds to the translation coordinate ε relative to the anchor;
  * an inactive cut corresponds to an anchored state, the only kind that survives arbitrarily long runs;
  * their arithmetic tube bound `e^(−ηt)` corresponds to the 2-adic loss of precision `2^(−a)` per letter.

  Type: ANALOGY. It shaped section 2.
* **Exact orders and invariant anticanonical systems** ("why finite generation attains the limit"): the certificate family is finitely generated over the free run parameters (K, J): finitely many heads per post-run depth, with transparent runs. Exact attainment for every integer, the analogue of their theorem, is precisely what remains open. Type: ANALOGY.
* **Kakeya maximal (3D) and Kakeya (4D)**: multiscale tube incidence. No usable transfer was found, so not used.

## 10. Next obligations

1. HYP-9240 for all D: bound the landing `N_D` or its odd-step excess.
2. Compile the four short run types (j = 3..6) into fixed-state certificate trees (the tiling's DP at −1). Then give the exact Haar coverage of the whole residual branch per depth.
3. A positive-integer statement for the universal state: the probability of absorption before both orbits reach 1, as a function of size.
4. Prove the Karlin–McGregor product law for R = 2 from THM-4581's weight.
5. Anchors at −5 and −17 (3-adic, contracting inverse cycles) and transitions into them. The tiling's ternary turns `G(m) = (8m−5)/9` are inverse −5-cycles.

## 11. Reproduction

All scripts are in `04-computation/experiments/twoanchor_20261007/`:
* `twoanchor_core.py [Dmax]`: parts A–E; ALL CHECKS PASSED (119,914 assertions, about 9 s).
* `tournament_chirotope.py`: the 4-Pfaffian criterion and Legendre frieze strips; ALL CHECKS PASSED.
* `child_compare.py` (+ .out): the J-table.
* `mersenne_shells.py`: Mersenne shells.
* `km_exponent.py` (+ `_n3000.out`, `_n20000.out`) and `km_analysis.py` (+ .out): the Karlin–McGregor law.
* Ad hoc exploration that produced the explicit heads (limit chains, BFS of absorbing patterns) is reproduced inside twoanchor_core.py parts B and D.
