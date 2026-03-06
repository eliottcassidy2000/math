# Deep Analysis of parity_tournaments_fixed.tex

**Author:** kind-pasteur-2026-03-05-S9
**Purpose:** Thorough audit for errors, extraction of geometric insights, catalog of referenced papers and investigation leads.

---

## A. Errors and Issues Found

### ISSUE-1: DR mod-4 proof (Thm 7.4) has broken arithmetic

**Location:** Lines 1800-1848
**Severity:** MEDIUM (result is likely correct but proof is invalid)

The proof attempts to determine the parity of |C_3| for doubly-regular tournaments from Moon's formula. The computation explicitly produces v_2(|C_3|) = -2, which the paper acknowledges as "impossible." The proof falls back to direct verification for n=3,7,11 only.

**Problem:** Moon's formula for DR tournaments gives |C_3| = n(n-1)(n-3)/24 (note: the paper writes n(n-1)(n-2)/24 on line 1813, which may itself be wrong — Moon's formula for the number of directed 3-cycles in a doubly regular tournament is actually C(n,3)/4 = n(n-1)(n-2)/24, but the 2-adic analysis then breaks as shown).

**Fix needed:** Either:
1. Use Kummer's theorem to properly compute v_2(C(n,3)) and v_2(C(n,3)/4)
2. Use the OCF directly: alpha_1 counts directed odd cycles (not just 3-cycles), so for DR_n, alpha_1 includes all odd cycle lengths. The mod-4 result should come from alpha_1 mod 2, not just |C_3| mod 2.
3. Downgrade from "Theorem" to "Verified Conjecture"

### ISSUE-2: SE-SYT formula gives non-integer (Thm 7.3)

**Location:** Lines 1739-1790
**Severity:** LOW (paper is honest about it, result verified for small cases)

The cited classical formula 2^((|lambda| - r(lambda))/2) gives 2^(3/2) for m=2 (n=5), which is not an integer. The paper says "the classical formula cited above requires a more careful statement."

**Root cause:** The formula the paper attempts to cite is likely from Stembridge's work on self-evacuating tableaux, but the correct formula for 2-core shapes is different. For staircase delta_{n-2} with n=2m+1, the correct count of SE-SYT is 2^(m^2) — this is verifiable directly but the paper doesn't have a clean proof.

**Fix:** Find the correct classical reference (likely Stembridge 1996, "Canonical bases and self-evacuating tableaux" or similar). Alternatively, give a direct computational proof using the promotion/evacuation theory for 2-cores.

### ISSUE-3: Transitive tournament uniqueness proof is incomplete (Prop 2.1)

**Location:** Lines 463-476
**Severity:** LOW (well-known correct result, just sketchy proof)

The proof claims: if T has a 3-cycle a→b→c→a, "inserting the remaining vertices in their unique consistent order gives at least two distinct Hamiltonian paths." This glosses over the fact that there may not be a "unique consistent order" — the remaining vertices may interact with {a,b,c} in complex ways.

**Better proof:** Use induction on n. If T is not transitive, there exist i,j with T(i,j)=1 and score(i) < score(j). The standard Rédei induction argument gives this. Alternatively, just cite Rédei (1934) directly for this well-known fact.

### ISSUE-4: Verification record outdated (lines 361-391)

**Severity:** LOW (just needs updating)

The table shows exhaustive verification only up to n≤6 and n=7 sampled. Our later work proved OCF as a polynomial identity at n≤8 (2^27 = 134M configs, all passing). Also proved at n≤7 as SymPy symbolic identity.

### ISSUE-5: Rajkumar et al. (arXiv:2110.05188) missing from bibliography

**Severity:** LOW (informational)

This paper is central to several of our investigation leads (flip classes, R-cones, sign rank) but is not cited in the tex bibliography.

---

## B. Geometric Insights Worth Preserving

### B1: Pin grid = staircase Young diagram delta_{n-2}

The bijection: tournament T with base path P_0 ↔ tiling t in {0,1}^m where m = C(n-1,2). Coordinates (r,c) encode arc (a,b) = (r+c+1, c). This makes the S_3 action on the triangular grid transparent as vertex-relabeling operations.

**Why this matters for OCF:** The subset convolution sum_S f_i(S)*g_j(R) decomposes over the Boolean lattice 2^{others}. The pin grid gives a GEOMETRIC interpretation of this lattice — subsets correspond to bit patterns on grid tiles, and the S_3 symmetry constrains which patterns are equivalent.

### B2: All hook lengths are odd (Thm 7.1)

For (r,c) in Grid(n): hook(r,c) = 2k'+1 where k' = n-1-r-c. This means delta_{n-2} is a 2-core (no hook of length 2). Consequences:
- The hook-length formula f^lambda = |lambda|! / prod(hooks) has only odd factors in the denominator (for this specific shape)
- This constrains the representation theory of S_m acting on tournaments
- The 2-adic valuation v_2(f^{delta_{n-2}}) is determined entirely by v_2(m!) where m = C(n-1,2)

### B3: The orphan concept as structural barrier

An "orphan" path P in Ham(T) satisfies P\v not in Ham(T-v). The orphan count equals H(T) - sum insact = H(T) - (B_v + S_v). The existence of orphans is WHY the insertion argument fails to prove Claim A directly. Understanding orphan structure is key to any bijective proof.

**Connection to OCF:** The Type-II positions (10 transitions in signature) exactly correspond to 3-cycles through v. Orphans arise from these 3-cycles. For n≥6, 5-cycles through v create additional orphan structure that the per-path identity can't capture.

### B4: The Q-function involution iota

The involution iota(L, P_L, P_R) = (V\L, bar{P_R}, bar{P_L}) is fixed-point-free because a fixed point requires L = V\L (impossible for distinct u,w). This is an instance of path reversal / complement duality in the pin grid. Under the tiling encoding, iota exchanges the "left" and "right" halves while bit-complementing.

### B5: Double-Burnside formula (Thm 6.1)

The orbit count under combined position-permutation and label-permutation actions gives the isomorphism class count. This is a COMPUTATIONAL tool: for each n, it reduces the space of tournaments to check from 2^m to roughly 2^m / n!.

---

## C. Papers Referenced — Full Investigation Leads

### PRIORITY HIGH

**Forcade 1973** — "Parity of paths and circuits in tournaments"
- Discrete Math 6 (1973), 115-118
- Original GF proof of F_2-invariance for k-block decompositions
- Our paper gives new combinatorial proof but Forcade's GF machinery may encode the polynomial identity we need for OCF
- **Action:** Read the 4-page paper carefully. Extract the generating function. Check if it lifts from F_2 to Z_2 (2-adic tower connection).

**El Sahili & Ghazo Hanna 2023** — "About the number of oriented Hamiltonian paths and cycles in tournaments"
- J. Graph Theory 102 (2023), 684-701
- Studies H(T) directly — bounds, formulas, structural results
- **Action:** Read for any results constraining H(T) that connect to OCF.

**Striker 2011** — "A unifying poset perspective on alternating sign matrices, plane partitions, Catalan objects, tournaments, and tableaux"
- Adv. Appl. Math. 46 (2011), 583-609
- The S_3-equivariance question (Open Problem 5) is COMPLETELY UNEXPLORED
- **Action:** Read the paper. Test equivariance computationally at n=3,4,5.

**Chapman 2001** — "Alternating sign matrices and tournaments"
- Adv. Appl. Math. 27 (2001), 290-298
- Direct ASM-tournament bijection. ASMs have determinantal formulas.
- **Action:** Read the paper. If H(T) maps to an ASM statistic with a known formula, this could give OCF.

**Rajkumar et al. 2021** (arXiv:2110.05188) — "A Theory of Tournament Representations"
- Not in bibliography! Must add.
- Flip classes, R-cones, locally transitive, sign rank
- **Action:** Add to bibliography. The flip-class proof strategy (INV-004) and mu(T) induction (INV-005) both come from this paper.

### PRIORITY MEDIUM

**El Sahili & Abi Aad 2020** — "Parity of paths in tournaments"
- Discrete Math 343 (2020), Art. 111695
- Decisive/concordant terminology comes from here
- **Action:** Read for mod-4 connections to our alpha_1 parity question.

**Schweser-Stiebitz-Toft 2025** (arXiv:2510.10659) — "The tournament theorem of Rédei revisited"
- Stronger Rédei for mixed graphs
- **Action:** Read for Q-Lemma generalization possibilities.

**Feng 2025** (arXiv:2510.25202) — "The dual Burnside process"
- Already partially explored by opus-S4b
- Transfer matrix symmetry connection (detailed balance / reversibility)
- **Action:** Deepen connection to our transfer matrix M[a,b].

### PRIORITY LOW

**Eplett 1979** — Self-converse tournament counts (specialized, less relevant to OCF)
**Frame-Robinson-Thrall 1954** — Hook-length formula (classical, well-understood)
**Stanley EC2** — Reference text (used for SE-SYT and general combinatorics)
**Rédei 1934** — Original theorem (starting point, nothing new to extract)

---

## D. Connections Between Paper and Recent Discoveries

| Paper concept | Recent discovery | Connection |
|---------------|-----------------|------------|
| Strategy 4: adjacency formula (line 1307) | THM-013/THM-014/THM-015 | Direct evolution — Strategy 4 became our arc-flip induction |
| Strategy 2: direct bijection | INV-007 (odd-cycle bijection) | Same problem, still open |
| Per-path identity (Thm 5.2) | Even-Odd Split Lemma | The per-path identity sums to give even-odd structure; the split is the aggregate form |
| Q-function / Forcade (Strategy 5) | Transfer matrix symmetry | Q-function counts 2-block cuts; transfer matrix M generalizes this to alternating-signed subset decomposition |
| S_3 x Z_2 symmetry | Pin grid encoding | Not yet exploited for OCF — could the S_3 symmetry constrain the polynomial identity? |
| Orphan paths (Def 3.6) | Swap involution unmatched paths (THM-014) | Orphans ≈ unmatched paths. The swap involution is the computational version of the orphan concept |

---

## E. Most Promising Unexploited Directions

1. **Forcade's generating function** (INV-023): A 4-page paper that proved F_2-invariance. Our whole project is about lifting this to Z. The GF itself has not been examined.

2. **Chapman's ASM bijection** (INV-021): ASMs have the ASM theorem (determinantal formula). If H(T) translates to an ASM statistic, we get a formula.

3. **Striker's unifying poset** (INV-008): The S_3-equivariance question connects our pin grid symmetry to ASMs and plane partitions. Completely unexplored.

4. **SE-SYT formula gap** (ISSUE-2): Tractable problem that would clean up Section 7 of the paper. Needs the right classical reference.

5. **Pin grid symmetry for OCF**: The S_3 action on barycentric coordinates has NOT been used for the OCF identity. The polynomial identity delta_H = delta_I lives on the Boolean lattice 2^{others}, which is a sublattice of the pin grid. Could the S_3 symmetry reduce the number of monomials to check?

---

Source: kind-pasteur-2026-03-05-S9
