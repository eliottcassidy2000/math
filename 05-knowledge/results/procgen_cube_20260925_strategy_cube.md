# The strategy cube of Althöfer's 3n±1 game: provability is decided by extremal cycle densities, not by drift

**Status: PROVED** (Theorems A–E: the exact criteria for bounded-lookahead convergence and for residue divergence certificates, in terms of the odd densities of the cycles of a finite graph; the sandwich ρ_min ≤ π_odd ≤ ρ_max that puts the drift strictly inside the provability window; the exceptional-dimension trichotomy; negation symmetry and its transport corollary) **+ FINITE-EXACT** (the classification of all 65,814 strategies of levels 1–5 into (i)–(iv), with an exact Python reference re-deriving every class, certificate and cycle; the exact distances δ_k, ε_k from Collatz to the two provable classes for k ≤ 9 and k ≤ 8, solver-certified by exact MaxSAT and cross-checked by CP-SAT and exhaustive search) **+ EMPIRICAL** (uniform samples at levels 6–12; the k = ∞ random-sign model). **Collatz remains OPEN.** No HYP or THM file was created; no novelty is claimed for the elementary lemmas.

Session `collatz-procgen-20260922`, strategy-cube lane, 2026-09-25. Scripts `04-computation/experiments/procgen_cube_20260925_{run,core,macro,boundary}.py` and the C engine `procgen_cube_20260925_engine.c`; output [procgen_cube_20260925.out](procgen_cube_20260925.out).

## 0. The answer in brief

Setting: a level-k strategy is σ : (odd residues mod 2^k) → {±1}, and T_σ(n) = n/2 (n even), (3n + σ(n mod 2^k))/2 (n odd). Level 1 is {Collatz, 3n−1}; level 2 is the Kuratowski lane's square.

1. **Classification** (every strategy, levels 1–5; (iii) = a positive cycle not through 1 found among starts ≤ 131072, cap 2^60):

   | k | strategies | (i) PROVED convergent | (ii) PROVED divergent | (iii) extra cycle | (iv) OPEN |
   |---|---|---|---|---|---|
   | 1 | 2 | 0 | 0 | 1 | 1 |
   | 2 | 4 | 1 | 1 | 1 | 1 |
   | 3 | 16 | 1 | 1 | 5 | 9 |
   | 4 | 256 | 16 | 3 (1 also iii) | 115 | 123 |
   | 5 | 65,536 | 1,052 | 32 (16 also iii) | 35,399 | 29,069 |

   * Every (i) strategy is **transitive**: its complete cycle list, from the finite check below the certificate threshold, is {1, 2}. Certificates need lookahead L ≤ 20 and threshold n_0 ≤ 22 at level 5.
   * Budget: the Terras DP runs to L = 40 for every strategy, with a zero-search to L = 400. By Theorem A the budget is immaterial: no strategy outside (i) has a certificate at *any* L.
   * No (i) strategy has an extra cycle at levels ≤ 5; (ii) and (iii) overlap (17 strategies).
2. **The exact provability boundary.** Let G_σ be the parity-transition graph on Z/2^k (§1) and c = log_3 2 = 0.6309 the critical odd density.
   * **Theorem A.** T_σ has a bounded-lookahead descent certificate ⟺ its exceptional set is empty at some finite level ⟺ **every cycle of G_σ has odd density < c**.
   * **Theorem B.** T_σ has a residue-class divergence certificate ⟺ **some closed (bottom) strongly connected class of G_σ has every cycle of odd density > c**. "A positive-drift trap" is necessary but *not* sufficient.
   * **Theorem C.** On every closed class, ρ_min ≤ π_odd ≤ ρ_max: the drift (π_odd log 3 − log 2) is an average of cycle densities, while provability needs the extreme ones on one side of c.
   * So provability is decided by two finite max/min-density cycle computations, and the drift sits strictly between them. Collatz has the widest possible window, [ρ_min, ρ_max] = [0, 1], at every level.
3. **The OPEN set.** Every OPEN strategy has a nonempty exceptional set (Theorem A). Most are Collatz-like (all closed classes of negative drift), but not all:
   * level 5: 28,064 Collatz-like and 1,005 "5n+1-like" (a positive-drift trap containing a contracting cycle; in 409 of them that cycle is the trivial cycle {1, 2} itself);
   * the drift ranges of (i) and of Collatz-like OPEN strategies overlap, so drift is not the boundary;
   * the exceptional dimension D ∈ [0.363, 1) on Collatz-like OPEN strategies; the thinnest is 0.363, Collatz's is h(log_3 2) = 0.950.
4. **Symmetry.** Negation ν exchanges Collatz and 3n−1 and fixes the mixed corners (checked).
   * Level k has 2^(2^(k−2)) self-dual strategies and the rest in dual pairs: 2+1, 4+6, 16+120, 256+32,640 for k = 2..5.
   * Everything residue-level (classes (i), (ii), densities, drift, D, |Bad_L|) is ν-invariant, verified on every pair.
   * **The OPEN set is not ν-closed**: it is carved out of the ν-invariant set "not (i), not (ii)" by the order-dependent predicate "no extra *positive* cycle". At level 5, 14,808 of 29,069 OPEN strategies are **sheet-obstructed** like Collatz (their mirror has a positive extra cycle), so the transport theorem forces any proof of them to be side-aware.
5. **Collatz's position.** Collatz is one flip away from (iii) at every level ≥ 3, but reaching (i) takes δ_k = 1, 2, 2, 4, 5, 9, 14, 23 flips for k = 2..9, and reaching (ii) takes ε_k = 1, 2, 3, 5, 10, 18, 30 for k = 2..8 (exact).
   * In Haar measure these are 0.500, 0.500, 0.250, 0.250, 0.156, 0.141, 0.109, 0.090 for (i) and 0.500, 0.500, 0.375, 0.313, 0.313, 0.281, 0.234 for (ii), both decreasing.
   * So Collatz is on the provable/unprovable boundary only at level 2 in the hypercube sense. In the Haar sense it approaches both provable classes, (i) faster, and plausibly lies on the boundary in the limit (OPEN).
6. **Microcosm → macrocosm.** The OPEN set **stabilises**; it neither shrinks nor grows.
   * The OPEN fraction is 0.444 at level 5 and 0.415–0.443 for sampled levels 6–12 (±0.008 to ±0.025). The k = ∞ model, an independent random sign for each odd integer, gives P(no extra cycle) = 0.435 ± 0.016.
   * The provable classes vanish: level 5 has 1.6% and 0.05%; samples have 3/4000 at k = 6 and none beyond.
   * Every residue-level statistic converges to Collatz's: π_odd → 1/2 with s.d. ≈ 0.54·2^(−k/2); D → 0.95; ρ_max → 1. So the drift boundary becomes "sharp" only in the degenerate sense that no strategy is near it.
   * What survives is the microcosm. OPEN versus (iii) is decided by whether the smallest integers close into a cycle: in 1000 random fields every extra cycle has minimum ≤ 1000, and 97.8% of level-5 witnesses have minimum ≤ 100.

## 1. Setting and the parity-transition graph

* **Strategies.** σ : {1, 3, ..., 2^k − 1} → {±1}; there are 2^(2^(k−1)) at level k.
  * A level-k strategy lifts to level k+1 (σ(r) := σ(r mod 2^k)) with the same map; counts "new at level k" exclude lifts.
  * Encoding: bit i of the mask ↔ residue 2i+1, bit 1 = minus. Mask 0 is Collatz; the all-ones mask is 3n−1.
* **u/d reading** (Althöfer lane): at an odd residue r, the step is *d* if 4 | 3r + σ(r) (at least two halvings follow) and *u* if exactly one. Collatz is d on 1 (mod 4), u on 3 (mod 4); 3n−1 is the reverse; (+,−) is all-d; (−,+) is all-u.
* **The graph G_σ.** Nodes Z/2^k. From s, the two edges go to the two lifts mod 2^k of T(s) mod 2^(k−1). Node weight w(s) = log(3/2) if s is odd, −log 2 if even.
  * A cycle with a odd nodes and length p is **expanding** iff 3^a > 2^p, i.e. iff its odd density a/p exceeds c = log_3 2, and **contracting** otherwise (3^a ≠ 2^p always).
  * ρ_max(X), ρ_min(X) denote the extreme odd densities of cycles inside a node set X.
  * Bottom SCCs (no edge leaves) are the **closed classes**. π_C is the stationary law of the uniform-edge chain on a closed class C; the drift is π_C(odd) log 3 − log 2 per step.
* **Exceptional set.** Bad_L = classes mod 2^(k+L) whose first L prefix multipliers all exceed 1 (3^(a_j) > 2^j for j = 1..L), and Bad_∞ = ⋂ Bad_L ⊂ Z_2.

**Lemma 1 (Terras bijection, PROVED).** For m ≥ k, T maps each class r + 2^m Z_2 affinely and bijectively onto T(r) + 2^(m−1) Z_2. Hence the classes mod 2^(k+L) correspond bijectively to the paths of length L in G_σ (the itineraries mod 2^k), and on each class T^j(x) = (3^(a_j) x + c_j)/2^j = (3^(a_j)/2^j)(x + Off_j) for j ≤ L, with a_j, c_j, Off_j = Σ_{odd i<j} σ_i 2^i / 3^(a_(i+1)) constant.

*Proof.* On r + 2^m Z_2 ⊂ r + 2^k Z_2 the sign is constant, so T is x ↦ (3x + σ)/2 or x/2, of slope 3/2 or 1/2, and the image is T(r) + 3·2^(m−1) Z_2 = T(r) + 2^(m−1) Z_2. Induct on L: the x ∈ s_0 + 2^k Z_2 whose image lies in the class of a tail path form the preimage, under the affine bijection T, of a class mod 2^(k+L−1), hence a class mod 2^(k+L). The recursion c_(i+1) = 3c_i + σ_i 2^i (odd step), c_i (even step) gives the closed forms. ∎

**Lemma 2 (periodic points, PROVED).** Each cycle γ of G_σ (length p, a odd nodes) carries exactly one x_γ ∈ Z_2 whose itinerary is γ repeated. It is the rational x_γ = c_γ/(2^p − 3^a). If γ is expanding, some rotation of γ has all prefix sums positive and the corresponding periodic point lies in Bad_∞.

*Proof.* T^p maps the class of γ (mod 2^(k+p)) affinely onto s_0 + 2^k Z_2, which contains it, with slope of 2-adic norm 2^p. The inverse is a contraction with one fixed point. The cycle lemma: with P_i the prefix sums and i* the last index attaining min P_i, starting at i* makes every partial sum positive when the total is positive. ∎

## 2. The theorems

**Theorem A (bounded-lookahead convergence, PROVED).** The following are equivalent.
* (a) Every cycle of G_σ is contracting: ρ_max(G_σ) < log_3 2.
* (b) Bad_L = ∅ for some L.
* (c) There are L and n_0 such that every integer n > n_0 has T^j(n) < n for some 1 ≤ j ≤ L.

Under (a), the least such L satisfies L_min ≤ 2^k (1 + log(3/2)/|μ|) + 1, with μ = ρ_max log 3 − log 2 < 0. Every positive orbit then enters a cycle meeting [1, n_0], and a finite computation lists all cycles.

*Proof.*
* (a) ⇒ (b). Decompose a walk of L steps into cycles plus a simple path of at most 2^k nodes. The cycles weigh at most μ·(their length), so S_L ≤ 2^k log(3/2) + μ(L − 2^k) < 0 for L beyond the bound: no class survives.
* (b) ⇒ (c). Since Bad_L is empty, each class mod 2^(k+L) has a first j ≤ L with 3^(a_j) < 2^j, and T^j(n) < n iff n > c_j/(2^j − 3^(a_j)). Take n_0 the largest of these finitely many thresholds.
* (c) ⇒ (a). If γ is expanding, the class of x_γ mod 2^(k+L) is in Bad_L (Lemma 2). On that class T^j(n) − n = ((3^(a_j) − 2^j)n + c_j)/2^j > 0 for all j ≤ L once n is large, contradicting (c).
* The last claim: strong induction gives that every orbit visits [1, n_0] infinitely often, hence is eventually periodic through a cycle meeting [1, n_0]. ∎

So "bounded-lookahead provable ⟺ exceptional set empty at some level" (the task's first example) is TRUE, and it is equivalent to a max-density-cycle computation. The density is computed exactly by Karp's algorithm on integer odd counts: all walks of length j have the same number of steps, so maximising the weight maximises a.

**Theorem B (residue divergence certificates, PROVED).** The following are equivalent.
* (a) Some closed class C of G_σ has every cycle expanding: ρ_min(C) > log_3 2.
* (b) There are a nonempty node set N closed under the edges of G_σ and an M ≥ 1 such that every path of length M inside N has 3^(a_M) > 2^M.

Then let θ > 1 be the least multiplier 3^(a_M)/2^M over M-paths in N, W = max over M-paths of max(0, −Off_M), and B = θW/(θ − 1). Every integer n in a class of N with n > B satisfies T^(jM)(n) − B ≥ θ^j (n − B), so its orbit diverges. The runner uses M with θ ≥ 2, so B ≤ 2W.
* Any T-invariant union of classes at any modulus 2^m reduces to (b), because by Lemma 1 T^(m−k) blows each class up to a full class mod 2^k.

*Proof.*
* (a) ⇒ (b). Take N = C, λ = min cycle mean > 0, and M with λ(M − 2^k) − 2^k log 2 ≥ log 2, so that θ ≥ 2.
* Divergence. By Lemma 1, T^M(n) = (3^(a_M)/2^M)(n + Off_M) ≥ θ(n − W) for n ≥ W, and N is closed. Hence e_j = T^(jM)(n) − B satisfies e_(j+1) ≥ θ e_j.
* (b) ⇒ (a). N contains a closed class C. If C had a contracting cycle γ, then along γ repeated M times, the M-blocks would each multiply by > 1 while their product, a power of the multiplier of γ, is < 1. ∎

**Corollary B′.** The task's second example, "a divergence proof iff a positive-drift trap", is FALSE as stated.
* Positive drift on the trap is necessary: Theorem C with ρ_min > c.
* It is not sufficient. At level 3, σ = (+,+,−,+) (duuu) has the single closed class {1,2,3,5,6,7} with drift +0.222 per step, yet it contains the contracting cycle 1 → 2 → 1, which is the trivial cycle {1, 2}. There are 1, 11 and 1,005 such OPEN strategies at levels 3, 4, 5; they are the 5n+1 phenomenon inside the cube.
* In particular, **a strategy with σ(1) = +1 can be proved divergent only through a trap avoiding the residues of 1 and 2**, because the trivial cycle contracts.

**Theorem C (sandwich, PROVED).** For every closed class C, ρ_min(C) ≤ π_C(odd) ≤ ρ_max(C).
* Hence (i) implies every drift is negative, and (ii) implies the trap's drift is positive.
* Neither converse holds: Collatz has π_odd = 1/2 but window [0, 1].

*Proof.* The stationary edge flow f(u→v) = π(u)/2 is a circulation, so it is a positive combination Σ λ_γ 1_γ of simple cycles. Then π(odd) = Σ λ_γ a_γ / Σ λ_γ p_γ, a mediant of the cycle densities. ∎

**Theorem D (exceptional dimension, PROVED; the formula is standard).** Let D(σ) = limsup (1/L) log_2 |Bad_L| be the growth exponent of the exceptional counts. For Collatz it is the Hausdorff dimension h(log_3 2) of Bad_∞ ([choice ladder](collatz_procgen_20260922_choice_ladder.md) §2). Numerically, D = inf over β ≥ 0 of log_2 ρ(A(β)), where A(β) is the adjacency matrix of G_σ with row s scaled by e^(β w(s)). Its value is:
* D = −∞ (the counts vanish) ⟺ (i);
* D = 1 ⟺ some closed class has positive drift ⟺ Bad_∞ has positive Haar measure;
* otherwise 0 ≤ D < 1.

*Proof sketch.*
* The first line is Theorem A.
* If a closed class has positive drift, it contains an expanding cycle (Theorem C). Prefixing that cycle, rotated and repeated, to a Markov continuation that stays above −B with positive probability (ergodic theorem) gives a positive-measure subset of Bad_∞.
* If every closed class has negative drift, then d/dβ log ρ at β = 0 is the drift of the closed classes by Perron–Frobenius perturbation, and transient classes have ρ < 2. So ρ(A(β)) < 2 for small β > 0, and the Chernoff bound |{S_L > 0}| ≤ 2^k ρ(A(β))^L gives D < 1.
* The equality D = inf_β log_2 ρ(A(β)) is the standard Chernoff/Cramér law for finite Markov chains: the upper bound is the Chernoff bound above, and the lower bound (tilting plus a ballot argument) is standard and not re-proved here.
* Numerically the formula matches the exact counts to L = 600 (§4, F3). ∎

**Theorem E (negation, PROVED).** For ν(x) = −x, T_(νσ) = ν T_σ ν with (νσ)(m) = −σ(−m mod 2^k), and s ↦ −s is an isomorphism G_σ ≅ G_(νσ) preserving parity.
* Hence classes (i) and (ii), ρ_min, ρ_max, the stationary laws, the drifts, D and every |Bad_L| are ν-invariant.
* The positive cycles of νσ are the negatives of the negative cycles of σ.
* **Corollary (transport).** Let σ be OPEN with νσ ∈ (iii). Then no side-blind (ν-invariant) argument proves that every positive orbit of σ reaches the cycle through 1, since it would prove the same for νσ, which is false. This is Corollary 7 of the [mod-192 note](collatz_procgen_20260924_inverse_tree_mod192.md) transferred to the cube.

**Proposition F (the two fixed-point obstructions, PROVED).** (i) forces σ(1) = +1 and σ(2^k − 1) = −1, written σ(−1) = −1 below. With σ(1) = −1, 1 is an expanding fixed point (a self-loop at node 1). With σ(−1) = +1, the 2-adic point −1 is one (the Collatz loop at node 2^k − 1).
* At level 3 the further condition "not (σ(3) = + and σ(5) = −)" kills the expanding 2-cycle {−1/5, 1/5}.
* This is the cube's version of THM-4470(5)'s "pair 0". There, both choices of a pairing near −1 rise, so **no** periodic pairing is provable. Here the sign pair (σ(1), σ(−1)) = (+, −) makes both fixed points ±1 contracting ({1,2} and {−1,−2}), which is why class (i) is nonempty at every level ≥ 2.
* At level 5, 3/4 of all strategies carry at least one of the two obstructions (49,152 = 3·16,384).

## 3. The classification, levels 1–5 (FINITE-EXACT)

**Method.** Two independent code paths:
* the C engine computes Tarjan SCCs, float Karp, Gaussian stationary laws, the Terras DP for |Bad_L| (L ≤ 40, zero-search to 400), a Perron-formula dimension, and a memoised cycle search with Brent detection;
* the Python reference re-derives, for every strategy at levels 1–5, the classes (i)/(ii) by exact integer Karp on odd densities.

For levels 1–4 it also re-derives:
* the simple-cycle enumeration;
* exact Fraction stationary laws and exact drift signs;
* the DP counts, including a brute force over every residue class for L ≤ 9 (L ≤ 7 on 16 masks at level 4);
* the full cycle search (same budget, same cap).

Every certificate is recomputed in Python with integers:
* for (i), the first L with Bad_L empty, the max-plus DP for the maximal threshold, and the finite check below it;
* for (ii), the least block length M with every M-path satisfying 3^(a_M) ≥ 2^(M+1), and the exact W.

Every reported cycle (114,320 at level 5) is re-verified. Budget for (iii): starts ≤ 131,072, cap 2^60; no start was unresolved.

**Level 3, all 16 strategies** (σ(1)σ(3)σ(5)σ(7); one closed class each; D = exceptional dimension):

| σ | u/d | class | ρ_max(G) | closed class: size, ρ_min, π_odd, ρ_max | drift/step | D | extra cycles (min) | ν-partner |
|---|---|---|---|---|---|---|---|---|
| ++++ | dudu | **IV** (Collatz) | 1 | 8, 0, 0.500, 1 | −0.144 | 0.950 | — | −−−− |
| −+++ | uudu | III | 1 | 8, 0, 0.500, 1 | −0.144 | 0.950 | 7 | −−−+ |
| +−++ | dddu | IV | 1 | 8, 0, 0.375, 1 | −0.281 | 0.804 | — | −−+− |
| −−++ | uddu | IV | 1 | 8, 0, 0.400, 1 | −0.254 | 0.846 | — | self |
| ++−+ | duuu | IV | 1 | 6, 1/2, 0.833, 1 | +0.222 | 1 | — | −+−− |
| −+−+ | uuuu | **II** (all-u) | 1 | 4, 1, 1, 1 | +0.405 | 1 | — | self |
| +−−+ | dduu | IV | 1 | 8, 0, 0.444, 1 | −0.205 | 0.876 | — | −++− |
| −−−+ | uduu | IV | 1 | 8, 0, 0.500, 1 | −0.144 | 0.950 | — | −+++ |
| +++− | dudd | IV | 2/3 | 8, 0, 0.417, 2/3 | −0.235 | 0.440 | — | +−−− |
| −++− | uudd | IV | 1 | 8, 0, 0.444, 1 | −0.205 | 0.876 | — | +−−+ |
| +−+− | dddd | **I** (all-d) | 1/2 | 8, 0, 0.333, 1/2 | −0.327 | −∞ | — | self |
| −−+− | uddd | IV | 1 | 8, 0, 0.375, 1 | −0.281 | 0.804 | — | +−++ |
| ++−− | duud | III | 1 | 6, 1/2, 0.667, 1 | +0.039 | 1 | 5, 13, 53 | self |
| −+−− | uuud | III | 1 | 6, 1/2, 0.833, 1 | +0.222 | 1 | 5 | ++−+ |
| +−−− | ddud | III | 2/3 | 8, 0, 0.417, 2/3 | −0.235 | 0.440 | 5 | +++− |
| −−−− | udud | III (3n−1) | 1 | 8, 0, 0.500, 1 | −0.144 | 0.950 | 5, 17 | ++++ |

* The level-3 microcosm already shows every phenomenon. The provable corners are the lifted all-d and all-u, both self-dual.
* The OPEN strategy +++− (Collatz with σ(7) = −) has the thinnest exceptional set, D = 0.440, and the smallest window, ρ_max = 2/3. Its only expanding simple cycle is 3 → 1 → 6 → 3, the residues of Collatz's cycle −5 → −7 → −10. Its mirror +−−− carries the 3n−1 cycle {5, 7, 10}, so it is sheet-obstructed exactly like Collatz.
* The OPEN −−−+ has Collatz's drift and D, and its mirror −+++ has a cycle through 7.

**Level 4.**
* (i): 16 strategies, all transitive. L_min histogram {2: 1, 4: 7, 5: 4, 7: 4}; thresholds n_0 ≤ 4; ρ_max ∈ {1/2, 3/5}; 2 self-dual and 7 dual pairs.
* (ii): all-u; uuuuduuu (d at 9); and uuuduuuu (d at 7, which also carries the cycle {5, 7, 10}). Their traps have ρ_min = 2/3 and 1.
* OPEN: 123 = 112 Collatz-like + 11 positive-drift.

**Level 5.**
* (i): 1,052 strategies, all transitive. L_min histogram {2: 1, 4: 67, 5: 254, 7: 368, 8: 164, 10: 132, 12: 34, 18: 24, 20: 8}; n_0 ≤ 22; all L_min ≤ 40, as Theorem A's bound guarantees.
* (ii): 32 strategies. Block lengths M ∈ {2, 8, 18, 29}; divergence thresholds 2W ≤ 18.8.
* (iii): 35,399. OPEN: 29,069.

## 4. The OPEN set (FINITE-EXACT at levels ≤ 5)

| level | OPEN | Collatz-like (all drifts < 0) | positive-drift (5n+1-like) | σ(1) = − or σ(−1) = + | neither | D range (Collatz-like) |
|---|---|---|---|---|---|---|
| 3 | 9 | 8 | 1 | 8 | 1 | 0.440–0.950 |
| 4 | 123 | 112 | 11 | 101 | 22 | 0.414–0.9998 |
| 5 | 29,069 | 28,064 | 1,005 | 22,440 | 6,629 | 0.363–1− |

* **Answers to §2 of the task.**
  * Every OPEN strategy has a nonempty exceptional set (Theorem A).
  * Not every OPEN strategy has negative drift: the positive-drift ones are 5n+1-like.
  * Every provably convergent strategy has an empty exceptional set at a finite level, by definition and Theorem A.
  * Divergence proofs exist exactly for uniformly expanding traps (Theorem B).
* **Drift is not the boundary.** At level 5 the per-step drift of (i) lies in [−0.327, −0.241], and that of Collatz-like OPEN strategies in [−0.320, −0.0002]; the ranges overlap. The largest trap drift of (ii) lies in [0.222, 0.406], and the positive-drift OPEN strategies reach 0.360.
* **Windows.** ρ_max(G) over OPEN strategies is mostly 1 (an all-odd cycle; 25,492 of 29,069), then 3/4, 2/3, 4/5, 5/6, 5/7, ... The least is 2/3.
* **Near-critical strategies and their giant cycles.**
  * Negative-drift strategies with drift within 0.02 of zero let orbits wander past 2^60: up to 122,423 of 131,072 starts escape the cap.
  * Their extra cycles are giants at the convergents and semiconvergents of log_2 3. At level 5 there are 2,094 cycles of length ≥ 300 in 1,742 strategies. The most frequent (p, a) are:
    * (1054, 665): 1,052 cycles in 902 strategies;
    * (485, 306);
    * (569, 359) = (485 + 84, 306 + 53);
    * (401, 253) = (485 − 84, 306 − 53);
    * 2 × (1054, 665), (317, 200), 4 × (84, 53), (653, 412), (1539, 971) = (1054 + 485, 665 + 306), ...
  * Every one of them has |a/p − log_3 2| < 3·10^(−4), and 1,834 have < 10^(−5). The reason is the exact identity a log 3 − p log 2 = −Σ log(1 + σ(x)/(3x)) over the odd points x of a cycle: a long cycle through large numbers is forced to the critical density. The same law pins the fragile pairs of THM-4470(6).
  * Budget sensitivity. Re-searching all 29,069 OPEN strategies with starts ≤ 2^20 (8 times the budget) finds **no** new extra cycle. 2,248 of them have some start escaping past 2^60: the near-critical and the positive-drift ones. The (iii)/OPEN split is stable in the start budget; only the value cap limits it, for near-critical strategies.
* **Where the witnesses live.** At level 5 the smallest extra cycle of a (iii) strategy has minimum ≤ 10 for 25,508 of the 35,399, in (10, 100] for 9,104, in (100, 1000] for 700, and above 10^4 for only 23. The most common minima are 5 (18,966 strategies), 7, 11, 17, 29, 43, 41, 53, 13.
* **Exceptional dimension** (F3). The Perron formula gives D = 0.949956 (Collatz), 0.804263 (+−++) and 0.439822 (+++−). The slopes of log_2 |Bad_L| over [400, 600] are 0.9459, 0.8005 and 0.4376, approaching D with the usual O(log L / L) lag.

## 5. The negation involution (FINITE-EXACT; Theorem E)

| level | self-dual | dual pairs | pairs (OPEN, OPEN) | pairs (OPEN, iii) = sheet-obstructed | self-dual OPEN |
|---|---|---|---|---|---|
| 1 | 0 | 1 | 0 | 1 (Collatz, 3n−1) | 0 |
| 2 | 2 | 1 | 0 | 1 | 0 |
| 3 | 4 | 6 | 2 | 4 | 1 |
| 4 | 16 | 120 | 41 | 33 | 8 |
| 5 | 256 | 32,640 | 7,074 | 14,808 | 113 |

* The square is exactly "dual pair + two self-dual": {Collatz, 3n−1} and the fixed (+,−), (−,+). **Both provable corners are the self-dual ones.**
* At levels 4–5, (i) splits into self-dual strategies and (I, I) pairs (2 + 7 pairs; 12 + 520 pairs).
* The residue-level invariants are identical on every pair (checked). The OPEN set is not ν-closed at any level.
* The transport theorem acts as follows. Residue data cannot prefer either member of a pair. In a mixed pair the OPEN member is the one whose extra cycles lie on its negative side. For the 14,808 sheet-obstructed OPEN strategies at level 5, like Collatz, any proof must use the order.
* The 7,074 (OPEN, OPEN) pairs and the 113 self-dual OPEN strategies are **sheet-free**: the transport theorem does not obstruct them, yet Theorem A excludes residue proofs. Their difficulty is purely INTEGRAL/DEFECT: a nonempty exceptional set that may contain no integer.
* In the macrocosm (§7) the sheet-obstructed share of OPEN is 0.52–0.59 (k = 6–12). This matches the k = ∞ prediction 1 − 0.435 = 0.565: there the two sides carry independent signs, so the mirror has an extra cycle with the same probability 0.565, independently.

## 6. Collatz's neighbourhood (FINITE-EXACT, solver-certified)

| k | residues | single flips → (iii) | δ_k (to (i)) | Haar | ε_k (to (ii)) | Haar |
|---|---|---|---|---|---|---|
| 2 | 2 | 0 | 1 | 0.500 | 1 | 0.500 |
| 3 | 4 | 1 | 2 | 0.500 | 2 | 0.500 |
| 4 | 8 | 3 | 2 | 0.250 | 3 | 0.375 |
| 5 | 16 | 5 | 4 | 0.250 | 5 | 0.313 |
| 6 | 32 | 6 | 5 | 0.156 | 10 | 0.313 |
| 7 | 64 | 11 | 9 | 0.141 | 18 | 0.281 |
| 8 | 128 | 12 | 14 | 0.109 | 30 | 0.234 |
| 9 | 256 | 12 | 23 | 0.090 | — | — |

* **Method.** Lazy-constraint MaxSAT (pysat RC2), cross-checked by CP-SAT for k ≤ 7.
  * Variables: x_r (flip residue r). For (ii), also y_s (node s in the trap), with hard clauses making the trap closed.
  * Lazy no-goods: an offending cycle (expanding for (i), contracting inside the trap for (ii)) exists whenever σ agrees on its odd nodes, so each clause is a necessary condition and the optimum is a lower bound. The first optimum with no offending cycle is optimal.
  * Every optimum is re-verified exactly: (i) witnesses by exact Karp, the descent certificate and the finite check (all transitive); (ii) witnesses by exact Karp on their closed classes.
  * Cross-checks: the exhaustive classification for k ≤ 5, and exhaustive Hamming-ball searches (δ_6 = 5 found at radius 5; ε_6 ≥ 10 from radius 9; δ_7 ≥ 7 from radius 6).
* **Reading.**
  * In the hypercube Collatz is adjacent to both provable classes only at level 2. At levels ≥ 3 it is adjacent only to (iii), where a single flip creates a cycle.
  * In the density picture it is as far from provability as possible: ρ_min = 0 (the loop at 0) and ρ_max = 1 (the loop at −1) at every level.
  * In Haar measure its distance to the provable classes shrinks monotonically, since lifting preserves the map: to (i) roughly like 2^(−k/3), and more slowly to (ii).
  * Whether δ_k/2^(k−1) → 0 is **OPEN**: it would put Collatz on the boundary of the provable region in the limit. It is the residue-periodic analogue of HYP-9136, the pairing family's "price of provable descent tends to zero". The contrast is real: in the pairing family THM-4470(5) forbids *every* periodic member, whereas here class (i) is nonempty and approaches Collatz.
  * The optimal flip sets are mostly u → d conversions at residues ≡ 3 (mod 4). For example, level 6 flips {7, 15, 27, 59, 63}: 63 = −1 is forced by the −1 loop, and 59 = −5 breaks the −5 cycle. Some optima also convert a d → u (level 5: {7, 15, 25, 31}, with 25 ≡ 1 mod 4).
  * The optimal (i) witnesses at levels 6–9 all have ρ_max = 5/8, just below log_3 2: the provable strategies nearest Collatz are near-critical.

## 7. Microcosm → macrocosm

Levels 1–5 are exhaustive; levels 6–12 are uniform samples (seed 20260925), each run together with its ν-partner under the same certificates and the same budget (starts ≤ 131,072, cap 2^60). Sections I, I3 and I4 of the output.

| k | strategies | (i) | (ii) | (iii) | OPEN (± 1 s.e.) | sheet-obstructed share of OPEN | positive drift | ρ_max = 1 | mean π_odd | s.d.(π_odd)·2^(k/2) | mean D (negative drift) |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | 2 | 0 | 0 | 0.500 | 0.500 | 1.000 | 0 | 1.000 | 0.5000 | 0 | 0.950 |
| 2 | 4 | 0.250 | 0.250 | 0.250 | 0.250 | 1.000 | 0.250 | 0.750 | 0.5833 | 0.500 | 0.950 |
| 3 | 16 | 0.0625 | 0.0625 | 0.313 | 0.563 | 0.444 | 0.250 | 0.813 | 0.5337 | 0.534 | 0.808 |
| 4 | 256 | 0.0625 | 0.0117 | 0.449 | 0.480 | 0.268 | 0.184 | 0.848 | 0.5202 | 0.515 | 0.890 |
| 5 | 65,536 | 0.0161 | 0.0005 | 0.540 | 0.444 | 0.509 | 0.105 | 0.878 | 0.5092 | 0.516 | 0.908 |
| 6 | 4,000 | 0.0008 | 0 | 0.563 | 0.436 ± 0.008 | 0.519 | 0.040 | 0.912 | 0.5051 | 0.523 | 0.927 |
| 7 | 3,000 | 0 | 0 | 0.572 | 0.428 ± 0.009 | 0.571 | 0.006 | 0.916 | 0.5031 | 0.536 | 0.939 |
| 8 | 2,000 | 0 | 0 | 0.577 | 0.423 ± 0.011 | 0.534 | 0 | 0.928 | 0.5002 | 0.520 | 0.944 |
| 9 | 1,500 | 0 | 0 | 0.563 | 0.437 ± 0.013 | 0.592 | 0 | 0.949 | 0.5004 | 0.542 | — |
| 10 | 1,000 | 0 | 0 | 0.574 | 0.426 ± 0.016 | 0.566 | 0 | 0.949 | 0.5003 | 0.565 | — |
| 11 | 600 | 0 | 0 | 0.585 | 0.415 ± 0.020 | 0.574 | 0 | — | 0.4988 | 0.568 | — |
| 12 | 400 | 0 | 0 | 0.558 | 0.443 ± 0.025 | 0.582 | 0 | — | 0.5004 | 0.551 | — |
| ∞ | 1,000 fields | — | — | 0.565 | 0.435 ± 0.016 | 0.565 (independence) | — | — | — | — | — |

* **The k = ∞ model.** For 1000 fields, P(no extra cycle with minimum ≤ X) is 0.613, 0.496, 0.436 and 0.435 at X = 10, 30, 100 and 1000, and stays at 0.435 up to X = 10^6. The number of extra cycles per field is distributed {0: 435, 1: 476, 2: 85, 3: 4}, mean 0.658. No start escaped 2^60.
* **The OPEN fraction stabilises** at the k = ∞ value. Levels 1–4 are too small to be random; from level 5 on the fraction sits at 0.42–0.44, within sampling error of 0.435.
  * This is not a shrinking to 0. A heuristic "log-many random cycles" would predict P(no extra cycle) → 0. It fails because cycles through large integers need 2^p ≈ 3^a (the identity in §4), so they are rare. P(no extra cycle with minimum ≤ X) is flat from X = 100 to 10^6.
  * It is not a growth either: every residue-level obstruction to OPEN, i.e. membership in (i) or (ii), disappears.
* **Provability vanishes** (EMPIRICAL, with a HEURISTIC reason).
  * (i) needs every cycle of G_σ contracting. A uniform σ realises, in expectation, about 2^(0.95p)/p expanding 2-adic periodic orbits of each period p ≲ k: the signed words of odd density > log_3 2, each realised with probability 2^(−a).
  * (ii) needs a closed class free of contracting cycles, while the positive-drift fraction itself vanishes.
  * Both fractions are 0 in the samples from k = 7 on.
* **The drift boundary.**
  * π_odd concentrates at 1/2 with s.d. ≈ 0.54·2^(−k/2), a CLT over the 2^(k−1) independent signs. The positive-drift fraction goes 0.25, 0.18, 0.10, 0.04, 0.006, 0 for k = 2..8, and the mean exceptional dimension rises to Collatz's 0.95.
  * So the drift boundary becomes sharp only in the degenerate sense that almost no strategy is near it. The provability boundary does not converge to it: the window [ρ_min, ρ_max] widens to [0, 1] for almost every strategy (ρ_max = 1 in 95% at k = 9–10).
* **The owner's theme, precisely.** The macrocosm (drift, exceptional dimension, cycle-density window) of a random strategy converges to Collatz's, and it is the same on both sheets. What distinguishes one strategy from another, and decides OPEN versus (iii), is the microcosm: the signs at the few smallest odd integers.
  * The single-flip census (H4) says the same: from level 7 on, about 12 residues of Collatz are "fragile", i.e. flipping one creates a small cycle (minima 5, 11, 31, 47, 61, 91). This is the cube's analogue of THM-4470(6)'s 24 fragile pairs.
  * Collatz itself is a typical macrocosm with a clean microcosm on the positive side: no small positive extra cycle, while its mirror 3n−1 has the extra cycles through 5 and 17.

## 8. What this does and does not say

* **What it does.**
  * The strategy cube has an exact, finite provability test on both sides (Theorems A, B). The drift, the quantity every heuristic uses, provably sits strictly inside the window [ρ_min, ρ_max] that decides provability (Theorem C).
  * Collatz has the maximal window at every level. The square's picture, "Collatz is one flip from both provable corners", is a level-2 accident: at higher levels Collatz is adjacent to intransitive strategies (iii), not to provable ones.
* **What it does not do.**
  * Nothing here bears on Collatz itself: Collatz is in (iv) at every level, and Theorem A says no bounded-lookahead residue argument can move it.
  * The OPEN strategies are open for the same reasons Collatz is. The sheet-obstructed ones need the order (Theorem E). The sheet-free ones have a nonempty exceptional set whose integer points are unknown (DEFECT/INTEGRAL). The positive-drift ones are 5n+1-type divergence problems.
* **Relation to the session.**
  * The square (Kuratowski lane) is the level-2 slice; its PROVED corners are the self-dual all-d and all-u.
  * The choice ladder's "choice collapses the exceptional set" is realised by the all-d strategy, i.e. choice frozen into residues; Theorem A says exactly which frozen choices keep the collapse.
  * The pair-0 obstruction of THM-4470(5) becomes the two fixed-point obstructions of Proposition F, which a sign strategy can avoid.

## 9. Reproduction

```bash
python3 04-computation/experiments/procgen_cube_20260925_run.py > 05-knowledge/results/procgen_cube_20260925.out
```

* **Requirements.** The runner compiles the C engine into `scratch/procgen_cube/strategy_cube/` (caches there are not to be committed) and needs `ortools` and `pysat`.
* **Cost.** One process at a time (CP-SAT with 2 workers); peak memory 514 MB; wall time 31 min (1849 s). The cost is dominated by the level-5 census, the 2^20 re-search, the Hamming-ball searches and the level-9/level-8 MaxSAT runs. The level-5 caches (`lvl5_*.txt.gz`) are reused when present.
* **Checks.** Every claim in the output is a `check(...)` that raises on failure.
* **Standalone solvers.** `python3 04-computation/experiments/procgen_cube_20260925_boundary.py i 2 3 4 5 6 7 8 9` and `... ii 2 3 4 5 6 7 8`.
* **Inputs used.**
  * [Kuratowski/Tait note](procgen_kuratowski_20260925_tait_kempe_triples.md) §1.7 (the strategy square);
  * [choice ladder](collatz_procgen_20260922_choice_ladder.md) (exceptional sets, the additive-choice zoo);
  * THM-4470 (pairing ladder, §5 and §6);
  * [mod-192 note](collatz_procgen_20260924_inverse_tree_mod192.md) (transport theorem, Corollary 7);
  * [Althöfer game note](collatz_procgen_20260922_althofer_game.md) (u/d moves, L1);
  * HYP-9136.
