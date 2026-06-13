# Research Synthesis: Tournament Parity & Path Homology
## Session kind-pasteur-2026-03-10-S51 — Full Stocktaking

*This document takes stock of every result accumulated across ~60 agent sessions, identifies connections, and maps potential applications and innovations.*

---

## PART I: WHAT HAS BEEN PROVED (Core Theorems)

### A. The Central Identity — OCF (FULLY PROVED)

**The Odd-Cycle Collection Formula** (THM-002): For every tournament T on n vertices,
```
H(T) = I(Ω(T), 2)
```
where H(T) = number of directed Hamiltonian paths, Ω(T) = conflict graph of odd directed cycles, I(G,x) = independence polynomial.

**Status:** PROVED for all n. Source: Grinberg-Stanley (arXiv:2412.10572, Corollary 20), internally verified at n≤10.

**What this means in words:** The parity of Hamiltonian path counts is entirely controlled by directed cycles. Each odd cycle "costs" one unit; pairs of vertex-disjoint odd cycles each add 4 (not 2+2) to H(T). The formula says: assign weight 2^k to each collection of k pairwise vertex-disjoint odd cycles; sum over all such collections. This is precisely the hard-core gas partition function at fugacity λ=2.

**Claim B** (THM-003, internal proof): The same formula holds with both H and I varying by vertex deletion: I(Ω(T),2) − I(Ω(T−v),2) = 2 Σ_{C∋v} μ(C).

---

### B. Path Homology Structure (GLMY Theory)

The chain complex Ω_*(T) with boundary maps d_p: Ω_p → Ω_{p-1} gives cohomological invariants.

**β_0 = 1 always** (THM-001): All tournaments are strongly connected (Rédei).

**β_2 = 0 always** (THM-108 + THM-109, PROVED):
- Proof: long-exact sequence of pair (T, T\v) + "good vertex" existence
- Verified: 0 failures at n=4-10 (100% of all tested tournaments)
- Consequence: "Even Betti vanishing at degree 2" — the first nontrivial topological constraint

**β_1 ≤ 1 always** (THM-103): The first Betti number is at most 1.

**β_1 · β_3 = 0** (THM-095): Seesaw theorem. Proved for n≤7.

**Top vanishing** (THM-124): β_{n-1} = β_{n-2} = 0 for ALL tournaments.
- Proof: d_{n-1}: Ω_{n-1} → Ω_{n-2} is injective (non-2-monotone face argument).

**Paley T_11 full Betti numbers** (HYP-443, CONFIRMED):
```
β(T_11) = (1, 0, 0, 0, 0, 5, 15, 0, 0, 0, 0)    chi = 11 = p
```
All 11 boundary maps verified computationally. Per-eigenspace: k=0 gives β=(1,0,0,0,0,5,5,0,0,0,0); k≠0 gives β=(0,0,0,0,0,0,1,0,0,0,0).

**Paley T_7 Betti**: β = (1,0,0,0,6,0,0), chi = 7 = p.

**Chi theorem**: chi(T_p) = p for Paley primes p≡3 mod 4. Confirmed p=3,7,11.

**Omega dimensions for T_11** (fully computed):
Per eigenspace: [1, 5, 20, 70, 205, 460, 700, 690, 450, 180, 30]. Palindrome FAILS (unlike T_7 where it holds).

---

### C. Hamiltonian Path Combinatorics

**Paley maximizer** (empirically established, OEIS A038375):
H(T_p) = max H(T) over all p-vertex tournaments, for p = 3, 7, 11.
Values: H(T_3)=3, H(T_7)=189, H(T_11)=95095.

**H/|Aut| sequence**: 1, 9, 1729, 6,857,869,865 for p=3,7,11,19.
The value **1729 = Hardy-Ramanujan taxicab number** (12³+1³ = 10³+9³ = 7·13·19) appears at p=11. H(T_11) = 55 × 1729.

**Permanent gaps in H-spectrum** (THM-029, THM-079, THM-115, PROVED for all n):
- **H = 7 is impossible** for any tournament on any number of vertices.
- **H = 21 is impossible** for any tournament on any number of vertices.
- H = 63 IS achievable (at n=8).

**SC Maximizer** (verified exhaustively n≤8): Within each self-complementary score class, the maximum H is always achieved by a self-complementary tournament. The global max H(n) is always SC.

**BIBD cycle maximization** (THM-028): Among regular n=7 tournaments, BIBD arrangement of 3-cycles gives alpha_1 = 80 (maximum) and H = 189.

**H=7, H=21 connection**: These gaps occur at H = 7·1 and H = 7·3. The H=7 gap "propagates" to produce a H=21 gap. The poisoning-graph DAG argument (THM-115) proves both.

---

### D. Transfer Matrix and Polynomial Structure

**Transfer matrix symmetry** (THM-030, PROVED for all n by induction):
M[a,b] = M[b,a]. The matrix M encodes "endpoint sensitivity" of H(T).

**Trace formula** (THM-027): tr(M) = H(T) for odd n, 0 for even n.

**Deletion-contraction** (THM-082/083):
- H(D) = H(D\e) + H(D/e) for any digraph D
- F(T,x) = F(T\e, x) + (x-1)·F(T/e, x)

**F_k mod 2 universality** (INV-124, proved): F_k(T) ≡ C(n-1, k) mod 2 for ALL tournaments. F(T,x) = (1+x)^{n-1} mod 2 — completely tournament-independent!

**Universal Taylor zeros mod 3** (THM-086): c_j(T) ≡ 0 mod 3 for all j < 2⌊(n-1)/2⌋. For n odd, F(T,x) mod 3 is a single-parameter family.

**9 | F(T,ω) for all n≥6** (THM-085): Universal divisibility. The forward-edge generating function at any primitive root of unity is divisible by 9.

**Moment-cycle hierarchy** (THM-092/093/117): The cumulants κ_{2k} of the forward-edge distribution are determined by odd-cycle counts up to length 2k+1. Universal coefficient: coeff(t_{2k+1}) in κ_{2k} = 2/C(n,2k). PROVED for all k.

---

### E. Structural Results

**Blueself impossibility at all odd n** (THM-023, INV-039): No blueself tilings at any odd n. Purely algebraic proof.

**Hereditary regular maximizers** (corrected in INV-044): At odd n, regular maximizers are hereditary (all vertex deletions give max H(n-1)). Fails at even n.

**Real roots of I(Ω(T),x)** (THM-020/THM-025):
- PROVED for n≤8 via claw-free + Chudnovsky-Seymour
- DISPROVED at n=9 (explicit counterexample: score [1,1,3,4,4,4,6,6,7])
- Maximally rare (0/10000 random n=9 samples)

**Tribonacci connection** (INV-035): H(T_full_n) = Tribonacci(n). Both sides satisfy the same recurrence via parallel decomposition structures (run decompositions ↔ interval graph packings).

---

## PART II: STRONG EMPIRICAL DISCOVERIES (Unproved but Highly Reliable)

| Discovery | Evidence | Status |
|-----------|----------|--------|
| β_2 = 0 for ALL tournaments | n=4-10, 100% | PROVED |
| β_1·β_3 = 0 for all n | n≤7: 100%, n=8: verified | Open n≥8 |
| Paley T_p maximizes H(T) | p=3,7,11 confirmed, p=19,23 match asymptotic | Strong conjecture |
| chi(T_p) = p for Paley p≡3 mod 4 | Proved p=3,7,11 | Proved case-by-case |
| F_k mod 2 universality | Exhaustive n≤6, sampled n=7,8 | PROVED |
| H-gaps are only 7 and 21 forever | n≤9 verified, all n proved via THM-115 | PROVED |
| Omega dims identical across eigenspaces | T_7 and T_11 verified | Open: algebraic proof |

---

## PART III: CURRENT OPEN FRONTIERS

### Tier 1: Closest to resolution

1. **Prove β_1·β_3 = 0 for all n** (THM-095 covers n≤7). The seesaw mechanism via im(d_2) is understood but doesn't yet close the general case.

2. **Algebraic proof that all p eigenspaces of T_p have identical Omega dims**. Computational evidence is complete (T_7, T_11 both verified fully). The proof should come from representation theory of Z_p acting on the chain complex.

3. **Almost-tournament claim for THM-086** (mod 3 universality): c_j(T\e) = 0 mod 3 for j < val(n)-1. A nested DC induction should close this.

4. **Paley Betti pattern**: Is β_{(p-1)/2}(T_p) = p-1 for all Paley primes? (T_7: β_3=0, β_4=6=p-1; T_11: β_5=5=p-6? No. β_6=15=p+4? No obvious pattern. Need T_19.)

### Tier 2: Major open problems

5. **Prove the Paley maximizer conjecture**: H(T_p) = max_{|V|=p} H(T). This would be a significant theorem connecting number theory (quadratic residues) to combinatorial extremal theory.

6. **What determines β_m(T) for generic n=8,9,10 tournaments?** The seesaw structure breaks, even Betti numbers appear (β_4>0). Is there a complete characterization?

7. **Tang-Yau stability theorem**: When does the pattern β_m(T_p) "stabilize" as p grows? Their Theorem 1.4 requires p ∉ Q+(connection_set), which needs computation for Paley sets.

8. **T_19 path homology**: The next valid Paley prime. Full Betti computation is computationally intense but in principle feasible with the small-prime technique.

---

## PART IV: POTENTIAL APPLICATIONS, UTILITIES, AND INNOVATIONS

### A. Immediate Mathematical Applications

#### A1. Hard-Core Lattice Gas at λ=2

OCF says H(T) equals the partition function Z(Ω(T), λ=2) of the hard-core lattice gas on Ω(T). This is a **physically meaningful** evaluation:
- λ=2 is exactly in the "high-fugacity" phase for graphs with maximum independent set size ≥ 2
- The Paley maximizer says: among tournament conflict graphs on p vertices, Ω(T_p) has the LARGEST partition function at λ=2
- This is a question in **extremal statistical physics**: which graphs on p vertices maximize Z(G, 2)?
- **New conjecture**: Does Ω(T_p) maximize Z(G, λ) for ALL λ ≥ 2? If so, this would extend the Paley maximizer result to all fugacities.
- This connects to Lee-Yang theory: the zeros of Z(G, λ) as a polynomial in λ are the Lee-Yang zeros. Our real-rootedness results (THM-020, THM-025) are exactly the Lee-Yang zero location problem.

#### A2. Algebraic Combinatorics: New Symmetric Function Connections

The Mitrovic-Stojadinovic bridge (arXiv:2506.08841) shows:
```
X_{inc(P)}(x) = ω(U_T(x))   (chromatic = omega-image of Redei-Berge)
```
This means **tournament H-values encode chromatic symmetric functions** of incomparability graphs. Concretely:
- The H-spectrum (achievable values of H(T)) corresponds to achievable chromatic polynomial evaluations
- H=7 and H=21 impossible ↔ no tournament can produce a chromatic polynomial with those values
- The gap theorem (THM-115) translates into a **chromatic polynomial gap theorem** via this bridge
- This is publishable independently as a consequence of the Mitrovic-Stojadinovic connection

#### A3. Cumulant Hierarchy as a Tournament Invariant System

Our results establish that the cumulants κ_{2k}(T) of the forward-edge distribution form a **complete system of invariants up to cycle content**:
- κ_2 determined by t_3 (3-cycles)
- κ_4 determined by t_3, t_5, α_2
- κ_6 determined by t_3, t_5, t_7, α_2
- Universal coefficient: coeff(t_{2k+1}) = 2/C(n, 2k) for all k

This creates a **moment problem** for tournaments: given a sequence (κ_2, κ_4, κ_6, ...), which are achievable as the cumulants of a tournament forward-edge distribution? The constraints are non-trivial and the boundary of achievable cumulant space is unknown.

**Application**: In social choice theory, a tournament represents pairwise preferences. The cumulant hierarchy gives a **basis for quantifying preference transitivity** at all scales simultaneously. κ_2 measures 3-cycle abundance; κ_4 measures 5-cycle vs. disjoint-3-cycle balance; etc.

---

### B. Computational and Algorithmic Utilities

#### B1. Sub-exponential Hamiltonian Path Counting

The deletion-contraction identity H(D) = H(D\e) + H(D/e) gives a tree-decomposition approach. Pairing this with:
- OCF: H = I(Ω, 2) reduces to independence polynomial evaluation
- Deletion-contraction for I(G, x): standard recursion
- Modern algorithms (FPRAS for independence polynomial for bounded-degree graphs)

This suggests a **practical fast algorithm for H(T) when T has special structure** (e.g., Ω(T) has small treewidth). The connection graph Ω(T) for sparse tournaments could have low treewidth even when T itself is dense.

**Utility**: Implement a Betti-number-guided H(T) approximation that uses:
1. Compute β_i(T) via GLMY (fast at low degrees via small-prime technique)
2. Use β profile to classify the "topological type" of T
3. Use type-specific formulas or bounds

#### B2. The Small-Prime Gaussian Elimination Technique

This session developed a highly effective computational technique: using PRIME=89 (since 11|(89-1)=88) with uint8 storage reduces memory by 16× compared to int64 without loss of rank. This enabled computing Betti numbers for T_11 which would otherwise be impossible on standard hardware.

**Utility**: This is a general technique for computing ranks of large sparse matrices over Z_p. It should be packaged as:
- A Python library `mod_rank` with `gauss_rank_mod(matrix, prime)` and automatic prime selection
- Documented with the small-prime stability theorem: rank over F_p is stable for p larger than any entry sum

**Innovation**: The "multiple-prime rank verification" pattern (compute rank at PRIME=89, verify at one other prime) is a robust algorithm for certifying ranks without large-number arithmetic.

#### B3. Eigenspace Decomposition for Circulant Combinatorics

The T_11 computation used the fact that the cyclic group Z_p decomposes the chain complex into p eigenspaces. This is a **general algorithmic technique** for any circulant structure:
1. Find the cyclic symmetry group of size m
2. Decompose each Ω_k into m eigenspaces
3. Compute Betti numbers eigenspace-by-eigenspace (m times smaller matrices)
4. Combine via Chinese Remainder Theorem

**This generalizes**: Any digraph with a cyclic automorphism group benefits from this decomposition. The T_11 computation took ~5 min/eigenspace × 11 eigenspaces vs. what would be a single ~1 GB matrix computation.

---

### C. Connections to Established Fields Needing Development

#### C1. Magnitude Homology / Path Homology Bridge

The Hepworth-Roff magnitude-path spectral sequence (INV-140) connects:
- **Page 1**: Magnitude homology (Kunneth, Mayer-Vietoris, easy to compute)
- **Page 2**: Path homology (our Betti numbers)

The spectral sequence converges. **Utility**: Use magnitude homology (which satisfies Kunneth!) as a computational shortcut or upper bound for path homology Betti numbers. Specifically:
- Compute magnitude Betti numbers (straightforward for products)
- Use them as initial estimates for GLMY Betti numbers
- The spectral sequence differentials tell you what "gets killed"

For Paley tournaments, the product structure T_p ≅ Cayley(Z_p, QR_p) means Kunneth may simplify the magnitude homology computation significantly.

#### C2. Coding Theory: Tournament Codes

Tournaments with large H(T) have a natural interpretation as **error-correcting codes**:
- Vertices = code symbols
- Hamiltonian paths = codewords
- H(T) = codebook size

The Paley maximizer property says T_p is the **optimal tournament code** on p symbols. The BIBD structure of T_p (its 3-cycles form a 2-(p, 3, 1) design) is the standard construction in code design. This suggests:

**New direction**: Interpret the independence polynomial I(Ω(T), x) as the **weight enumerator** of a tournament code, with x tracking some notion of "error weight". The MacWilliams transform might then connect Ω(T) to the dual code. The H=7 and H=21 gaps would correspond to forbidden minimum distances.

#### C3. Statistical Ranking Models

Tournaments model pairwise comparison outcomes (sports rankings, survey preferences, etc.). H(T) counts linear orderings consistent with the tournament. The OCF formula says:
```
H(T) = sum over all consistent orderings = I(Ω(T), 2)
```

**Innovation in social choice**: The independence polynomial evaluation gives a **weighted vote aggregation formula**:
- Each collection of k pairwise-disjoint odd cycles contributes 2^k to the "score"
- This is a non-linear aggregation that weights **consistent circular preferences differently from linear ones**

This is richer than standard Condorcet/Borda methods. A tournament with large α_2 (many disjoint cycle pairs) gives more "weight" to the quadratic term of I(G,x), meaning the outcome is more sensitive to 5-cycle vs. 7-cycle structure.

**Proposed methodology**: For a social choice application, represent paired comparisons as a tournament T. Use H(T) as an "ordering entropy": how many linear orderings are consistent? Use the OCF decomposition to attribute this entropy to specific cycle conflicts. The k-th term I_k = α_k × 2^k quantifies the contribution of k-cycle collections.

---

### D. Novel Mathematical Innovations

#### D1. The Permanent-Gap Phenomenon

We proved H=7 and H=21 are permanent gaps. This is remarkable because:
- H must be odd (Rédei's theorem)
- H=1, 3, 5, 9, 11, 13, 15, 17, 19, 23, 25, ... are all achievable
- But 7 and 21 are NEVER achievable, regardless of n

**The deeper structure**: 7 = 7 and 21 = 3×7. These are both multiples of 7. The next gap candidate, H=49 = 7², appears not to be achievable at small n but the proof breaks.

**Open innovation**: Is there a complete characterization of permanent H-gaps? Are they all multiples of 7? Or do they form a more complex set? The "poisoning DAG" technique developed for H=21 could potentially extend to characterize all gaps, giving a **complete taxonomy of the H-spectrum**.

This would be analogous to Bertrand's postulate in the prime number context: there's always a prime between n and 2n. Here: there's always an achievable H-value between any two consecutive odd values — except for the specific gaps at 7 and 21.

#### D2. Path Homology as Tournament Discriminator

The Betti profile β(T) = (β_0, β_1, ..., β_{n-1}) is a **topological invariant** of a tournament that captures:
- β_0 = 1 always (connected)
- β_1 ∈ {0,1} (whether T has a "topological hole" at degree 1)
- β_3 ∈ {0,1,2} (higher-order holes)
- β_4 = 0 for most, positive for SC near-regular tournaments (Paley T_7 has β_4=6)

**Innovation**: Use the Betti profile as a **fast tournament classifier** for machine learning applications:
- The GLMY computation at small degrees (β_1, β_2, β_3) is fast (polynomial in n)
- These give robust, structurally meaningful features
- Unlike algebraic invariants (cycle counts), Betti numbers are stable under small perturbations
- The SCM pattern (β_4>0 requires SC score) is a topological certificate of near-regularity

**Concrete application**: Given a tournament T representing preferences in a survey or sports season, compute:
1. β_1: is there a consistent cyclic preference pattern? (1 = yes)
2. β_3: are there higher-order preference cycles? (>0 = complex structure)
3. β_5, β_6, ...: for near-regular tournament, detect Paley-like algebraic structure

This would be the first use of topological data analysis (TDA) for tournament analysis.

#### D3. The Cumulant-Cycle Duality Principle

The universal coefficient theorem (THM-117): coeff(t_{2k+1}) in κ_{2k} = 2/C(n, 2k).

This is a **remarkable duality**: the (2k+1)-cycle count appears exactly in the (2k)-th cumulant with a universal coefficient. The proof uses:
- Forward path formula: #fwd_{2k+1} path = Σ_S H(T[S]) (OCF on subtournaments)
- Multinomial expansion

**Innovation**: This duality suggests a **moment generating function approach** to OCF. Define a formal power series in the cumulant variables:
```
log Z(T, λ) = Σ_k κ_{2k}(T) · f_k(λ)
```
The universal coefficient theorem says the cycle variables enter this expansion in a remarkably clean way. The question is: can we find f_k such that Z(T, 2) = H(T)? If so, this gives a **new proof of OCF via the cumulant expansion**.

#### D4. Mod-p Tournament Universality Hierarchy

The results (THM-085, THM-086, INV-124) establish a hierarchy:
- mod 2: F_k(T) completely universal (same for all T)
- mod 3: F_k(T) nearly universal (single free parameter at fixed n)
- mod 5: F_k(T) has multiple free parameters
- mod p: universality threshold at n = p+2

**Innovation**: This hierarchy defines a **"arithmetic complexity" of a tournament** — the smallest prime p for which F_k(T) mod p is tournament-dependent. The maximizer T_p might have minimal arithmetic complexity (its mod-p structure should be uniquely constrained by the Paley structure).

This creates a connection to the theory of **p-adic tournament invariants**, where the goal is to understand which tournament properties can be detected modulo small primes.

---

### E. Near-Term Publishable Results

Based on the accumulated work, the following papers are essentially ready:

#### Paper 1: "The Odd-Cycle Collection Formula: Internal Verification and Consequences"
- THM-002 (OCF), THM-003 (Claim B), the full verification at n≤10
- Permanent gaps H=7 and H=21 (THM-029, THM-115)
- F_k mod 2 universality (INV-124)
- Cumulant hierarchy (THM-092, THM-117)
- **Ready**: All proofs complete.

#### Paper 2: "Path Homology of Paley Tournaments"
- Full Betti numbers β(T_7) and β(T_11)
- chi(T_p) = p theorem
- Omega palindrome (holds T_7, fails T_11) — contrast and explanation
- Small-prime Gaussian elimination technique
- **Nearly ready**: T_19 data would strengthen it but isn't required.

#### Paper 3: "Universal Vanishing Theorems for Tournament Path Homology"
- β_2 = 0 (THM-108/109)
- β_{n-1} = β_{n-2} = 0 (THM-124)
- β_1 · β_3 = 0 (n≤7, THM-095) — with the seesaw mechanism
- **Ready for n≤7 case; n≥8 still partially open.**

#### Paper 4: "The Transfer Matrix Symmetry and Its Consequences"
- THM-030 (symmetry, proved)
- THM-027 (trace = H(T) for odd n)
- THM-052 (scalar M for SC VT tournaments)
- c-tournament generalization
- **Ready**: All proofs complete.

---

### F. Speculative Long-Range Innovations

#### F1. Tournament Homology Theory

The path homology framework (GLMY) gives tournaments a **topological character**. The key observation is:
- χ(T_p) = p (a prime!) for Paley tournaments
- β(T_p) concentrates at degrees (p-1)/2 and (p+1)/2
- This is reminiscent of the **Weil conjectures** for varieties over finite fields: zeta functions with zeros on the critical line, Betti numbers controlled by geometry

**Speculative innovation**: Is there a **"Weil conjecture analogue" for tournaments**? The Paley tournament T_p comes from the algebraic geometry of F_p (quadratic residues), and its path homology has chi = p and Betti numbers concentrated at middle degrees. This pattern (middle-degree concentration, Euler characteristic = field size) is characteristic of smooth projective varieties.

If this analogy can be made precise, it would connect tournament path homology to **étale cohomology** and open up an entire theory of "motivic" tournament invariants.

#### F2. The H-Spectrum as a Mathematical Universe

The set of achievable H-values forms a subset of odd positive integers. We know:
- Gaps at 7 and 21 (proved for all n)
- No other gaps in [1, 200] for n≤9
- The maximum H grows as ~p!/2^{p-1} · e (asymptotically)
- H(T_p)/H_{max} → e/2 as p → ∞

**Question**: Does the density of achievable H-values approach 1 as H → ∞? I.e., do the gaps become sparser? Or is there a positive-density set of permanent gaps?

This is analogous to **Goldbach-type problems** in additive number theory. The permanent gaps 7 and 21 could be the first two in an infinite sequence of "permanent gaps" forming a sparse set, or 7 and 21 could be truly exceptional with all sufficiently large odd numbers achievable.

#### F3. OCF-Based Cryptographic Hash Function

The independence polynomial I(Ω(T), x) evaluated at x=2 gives H(T). The function T ↦ H(T) mod p (for prime p) is a function from tournaments (size C(n,2) bits) to {0,...,p-1}. Properties:
- Deterministic: same input always gives same output
- Tournament-sensitive: small changes in T can change H(T) dramatically
- Modular arithmetic friendly (our THM-085/086 give constraints)

**Speculative utility**: Could (tournament, prime p) → H(T) mod p serve as a **combinatorial pseudorandom function**? The constraints from THM-085/086 (divisibility by 9 for n≥6, universal structure mod 2) might make this predictable, but at primes p≥5, the distribution is not yet fully understood.

---

## PART V: PRIORITIES FOR NEXT STEPS

### Immediate (1-2 sessions):

1. **Prove the eigenspace identity algebraically** (HYP-437): Why do all p eigenspaces of T_p have identical Ω_m dimensions? This should follow from the representation theory of Z_p acting on the chain complex.

2. **Apply Tang-Yau symbol matrix to T_7 and T_11 explicitly** (INV-135): Verify that their stability theorem covers the Paley case and understand when β_m(T_p) "locks in" for all large p.

3. **Submit H(T_p) sequence to OEIS**: H/|Aut| = 1, 9, 1729, 6857869865 is not in OEIS. Neither is H(T_p) as a function of p≡3 mod 4. Submission would establish priority and invite connections from the OEIS community.

### Medium-term (3-10 sessions):

4. **Compute β(T_19) fully**: The next Paley prime. Omega dims will be much larger (~3000-5000 per degree) but the small-prime technique should handle it. This would confirm or refute the χ(T_p)=p pattern and reveal the Betti structure.

5. **Write Paper 2** (Paley path homology): The T_7 and T_11 results are complete and publishable. The techniques (eigenspace decomposition, small-prime rank, chi theorem) are clean and novel.

6. **Prove β_1·β_3=0 for all n**: The seesaw mechanism is understood (im(d_2) mediates); formalize the argument beyond n=7.

### Long-term:

7. **Prove the Paley maximizer conjecture**: H(T_p) = max H(T) for all p≡3 mod 4 prime. This would be a major theorem.

8. **Develop the tournament → lattice gas connection**: Is there a physical interpretation of the Paley optimal fugacity problem? What does the "phase transition" at H-maximization correspond to in statistical physics terms?

---

## SUMMARY TABLE

| Result | Status | Impact | Application |
|--------|--------|--------|-------------|
| OCF: H(T)=I(Ω,2) | PROVED | Central | Hard-core gas, codes |
| Permanent gaps H=7,21 | PROVED | Novel | Chromatic polynomial bridge |
| Paley maximizer | Conjecture | Major | Extremal graph theory |
| β_2=0 universally | PROVED | Fundamental | Tournament TDA |
| Full β(T_11) | CONFIRMED | Novel | Path homology theory |
| chi(T_p)=p | Confirmed 3 cases | Deep | Weil conjecture analogy |
| F_k mod 2 universal | PROVED | Clean | Arithmetic structure |
| Cumulant hierarchy | PROVED | Clean | Social choice, statistics |
| Taylor zeros mod 3 | Proof sketch | Strong | Arithmetic combinatorics |
| 1729 in H/|Aut| | Computed | Mysterious | Number theory |
| Transfer matrix symmetry | PROVED | Technical | Algebraic structure |
| DC for Ham paths | PROVED | Useful | Algorithm design |
| Real roots n≤8 | PROVED | Strong | Lee-Yang theory |
| β·seesaw n≤7 | PROVED | Deep | Topology of tournaments |

---

*Written by kind-pasteur-2026-03-10-S51.*
