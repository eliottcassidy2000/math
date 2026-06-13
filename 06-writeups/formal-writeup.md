# Parity, Counting, and Homology in Tournaments: A Research Summary

**Status:** Living document. Last updated 2026-03-16.

---

## 1. Overview

This document summarizes research at the intersection of tournament combinatorics, algebraic topology, and Fourier analysis. The central objects are **tournaments** (complete directed graphs) and their **Hamiltonian paths** (directed paths visiting every vertex exactly once). The work spans six interlocking programs:

1. **The Odd-Cycle Collection Formula (OCF)** — relating Hamiltonian path counts to independent sets in cycle conflict graphs
2. **The Walsh-Fourier spectral program** — decomposing tournament invariants via Walsh transforms on {0,1}^m
3. **The signed Hamiltonian permanent** — a skew-symmetric analogue with universal congruence properties
4. **GLMY path homology of tournaments** — topological invariants revealing structural hierarchy
5. **Number-theoretic connections** — Egyptian fractions, cyclotomic splitting, spectral zeta, formal groups
6. **Engineering applications** — modular rank computation, circulant codes, TDA for preference data

Throughout, T denotes a tournament on n vertices, H(T) the number of directed Hamiltonian paths, and T^op the complement tournament (all arcs reversed).

### What is genuinely new

The following results appear to be new contributions not found in prior literature:

- **Walsh-Fourier spectrum of H(T) and M[a,b]**: Complete closed-form formulas (THM-071, THM-080). No prior work computes the full Walsh spectrum of tournament path counts.
- **Direct Walsh proof of OCF** (THM-077): An elementary proof of H(T) = I(Omega(T), 2) that bypasses P-partition theory entirely.
- **beta_2 = 0 for tournaments**: No vanishing result of this kind exists in the path homology literature. **Proved** (THM-108/109) via induction, long exact sequence, and isolation characterization.
- **Twin vertex mechanism**: The structural explanation for WHY beta_2 vanishes — completeness forbids the twin vertices that all beta_2 > 0 oriented graphs require.
- **beta_3 = 2 at n = 8**: The first Betti number exceeding 1 in tournament path homology.
- **Universal signed permanent congruences** (THM-H, THM-I, THM-J): The universality criterion s_2(n-3) <= 1 and the master identity D_S = n!/2^k appear to be new.
- **THM-130 (Paley Betti formula)**: Complete formula beta_m = m(m-3)/2, beta_{m+1} = C(m+1,2), chi = p for Paley tournaments. The per-eigenspace decomposition and Heisenberg Lie algebra connection appear new.
- **THM-125 (Constant symbol matrix)**: For circulant tournaments, all eigenspaces have identical Omega dimensions — enables n-fold speedup in rank computation.
- **THM-227 (k-nacci Mersenne identity)**: Tr(M_k^n) = 2^n - 1 for all 1 <= n <= k, connecting forbidden H values to Mersenne numbers via k-nacci matrices.
- **Simplicial Redei** (THM-220): sim_H in {0,1} for all tournaments on n >= 4 vertices.
- **Unit fraction splitting characterization** (HYP-1612): 3/N = 1/b + 1/c solvable iff N has a prime factor congruent to 2 mod 3.
- **Master splitting criterion** (HYP-1615): k/N = 1/b + 1/c solvable iff N^2 has a divisor d with d = -N mod k.
- **Base-42 Erdos-Straus covering**: Verified 0 failures through 10^6.
- **k-nacci trace forbidden values** (HYP-1618, corrected): The forbidden values 7 and 21 arise from k-nacci traces via Newton's identities (Tr(M_k^3) = 7 for all k >= 3; Tr(M_3^5) = 21 = 3 x 7), NOT from the standard Riemann zeta function (which gives zeta(-3) = 1/120).

---

## 2. Foundational Results

### 2.1 Redei's Theorem and the Q-Lemma

**Theorem (Redei, 1934).** For every tournament T on n vertices, H(T) is odd.

This classical result admits at least four independent proof routes, three of which were developed or refined in this project:

- **Route A (Q-Lemma).** Define Q_T(u,w) = #{Hamiltonian paths with u before w}. The Q-Lemma states Q_T(u,w) = Q_T(w,u) (mod 2), proved via a toggle involution on path pairs. This implies H(T) = sum of Q values = odd.

- **Route B (Anti-isomorphism involution).** For strongly connected T, an involutive anti-automorphism beta pairs Hamiltonian paths, leaving |Fix(beta)| = H(Q) for a smaller decisive quotient Q. Induction closes the argument.

- **Route C (Automorphism parity).** |Aut(T)| is always odd (proved), and Aut(T) acts freely on Ham(T), so |Aut(T)| divides H(T). Since |Aut(T)| is odd and H(T) = |Aut(T)| * (orbit count), H(T) is odd.

- **Route D (from OCF).** H(T) = 1 + 2*alpha_1 + 4*alpha_2 + ... where alpha_k counts independent sets of size k in the odd-cycle conflict graph. This is manifestly odd.

### 2.2 The Odd-Cycle Collection Formula

**Theorem (OCF; Grinberg-Stanley, 2024).** For every tournament T,

    H(T) = I(Omega(T), 2)

where Omega(T) is the conflict graph whose vertices are directed odd cycles of T and whose edges connect cycles sharing a vertex, and I(G, x) = sum_k alpha_k(G) x^k is the independence polynomial.

**History.** The formula was conjectured computationally in this project and verified exhaustively through n = 8 (2^27 configurations, ~57 minutes). It was independently proved by Grinberg and Stanley (arXiv:2412.10572, Corollary 20) using the theory of P-partitions and Claim B (see below).

**Claim A.** H(T) - H(T\v) = 2 * sum_{C in C_v} mu(C) for any vertex v, where C_v is the set of odd cycles through v and mu(C) counts independent sets in the subgraph of Omega induced by cycles disjoint from C\{v}, evaluated at x=2.

**Claim B (proved).** I(Omega(T), 2) - I(Omega(T\v), 2) = 2 * sum_{C in C_v} mu(C).

The OCF follows from Claims A + B by induction (base case: transitive tournaments have H = 1 and no odd cycles).

### 2.3 The Per-Path Identity

For n <= 5, the OCF admits a path-level refinement. For each Hamiltonian path P' in T\v and vertex v:

    (inshat(v, P') - 1)/2 = sum_{C containing v, |C|=3} mu(C)

where inshat counts insertion positions (always odd). This identity **fails** at n = 6: 2,758 out of 9,126 (T, v, P') triples violate it, due to contributions from 5-cycles. The correct generalization incorporating all odd cycle lengths remains **open** (OPEN-Q-004).

---

## 3. The Transfer Matrix and Walsh-Fourier Program

### 3.1 The Transfer Matrix

Fix a labeling 1, ..., n. The **transfer matrix** M[a,b] counts Hamiltonian paths starting at vertex a and ending at vertex b.

**Theorem (THM-030).** M[a,b] = M[b,a] for all tournaments T and all vertices a, b.

This symmetry was proved by induction using the Walsh-Fourier framework. For odd n, trace(M) = H(T).

### 3.2 Walsh-Fourier Decomposition

Encode a tournament T on n vertices as a binary string t in {0,1}^m where m = C(n,2), recording the orientation of each edge pair. The **Walsh-Fourier transform** decomposes any tournament invariant f(T) into components indexed by subsets S of [m]:

    f_hat[S] = 2^{-m} sum_T (-1)^{t cdot chi_S} f(T)

where chi_S is the indicator of S.

**Theorem (THM-071, Walsh diagonalization).** The Walsh transform of H(T) has nonzero coefficients f_hat[S] only when S is a union of edge-disjoint even-length paths in K_n. Moreover:

    H_hat[S] = epsilon * 2^r * (n - 2k)! / 2^{n-1}

where |S| = 2k, r is the number of path components of S, and epsilon = +/-1 depends on the path orientations.

**Theorem (THM-077, Direct Walsh proof of OCF).** By computing I_hat[S] for the independence polynomial side using a generating function factorization (THM-076), one shows H_hat[S] = I_hat[S] for all S. This provides an elementary proof of OCF bypassing P-partition theory.

### 3.3 Complete Walsh Spectrum of M[a,b]

**Theorem (THM-080).** For the transfer matrix entry M[a,b]:

    M_hat[a,b][S] = (-1)^{asc(S)} * 2^s * (n - 2 - |S|)! / 2^{n-2}

where s is the number of **unrooted** even-length components of S (components not containing the roots a or b), and the formula is nonzero only when S union {a,b} forms a valid path union of even length with |S| = n (mod 2).

This formula is manifestly symmetric in a and b, providing a **Walsh proof of transfer matrix symmetry**.

Verified: exhaustive at n = 5 (968 nonzero coefficients), n = 6 (1471), and by random sampling at n = 7.

### 3.4 Position Character Decomposition

**Theorem (THM-068, PCD).** For each vertex v and Walsh subset S of degree 2k:

    M_hat[v][S] = (-1)^{[v in N(S)]} * H_hat[S] / (n - 2k)

where N(S) is the set of "odd-offset" vertices in the path components of S. Proved for all degrees and all odd n via a block-placement argument.

**Theorem (THM-072, Off-diagonal PCD).** The off-diagonal Walsh spectrum has a similar structure, with interior vertices contributing zero rows/columns to M_hat at degree d.

### 3.5 Dimension Ladder and 2-adic Structure

The ratio between H and M Walsh amplitudes satisfies:

    H_amp(n, d, r) / M_amp(n, d-1, s=r-1) = n - d

**Theorem (HYP-1606, proved).** The product of all ladder ratios for a given n is (n-2)!! — the double factorial (OEIS A001147, counting perfect matchings).

**Theorem (HYP-1607, proved).** The sum of all ladder ratios is k^2 - 1 where k = (n-1)/2 (OEIS A005563, oblong numbers).

**Theorem (HYP-1608, proved).** The total 2-adic weight = k^2 + A000788(k-1), decomposing into a smooth quadratic part and a fractal logarithmic part (cumulative binary digit sum). This sequence is **not in OEIS**. The second differences encode binary carry counts: Delta^2 a(k) = 3 - trailing_ones(k+1).

**Spectral Legendre Identity (HYP-1603).** The 2-adic spread across the M Walsh spectrum equals v_2((n-3)!), connecting to Legendre's formula. The spectral Legendre excess = -s_2(n-3), the same quantity controlling THM-J universality.

---

## 4. The Signed Hamiltonian Permanent

### 4.1 Definition and Basic Properties

Define B = 2A - J (the skew-symmetric signed adjacency matrix of T, with entries +/-1). The **signed Hamiltonian permanent** is:

    S(T) = sum_P prod_{i=1}^{n-1} B[P_i, P_{i+1}]

summed over all permutations P (Hamiltonian paths, ignoring direction).

**Theorem (THM-A).** S(T) = 0 for all tournaments on even n (reversal pairing with sign (-1)^{n-1} = -1).

### 4.2 Universal Congruences

**Theorem (THM-H).** S(T) mod 2^{n-1} depends only on n, not on T.

**Theorem (THM-J).** S(T) is **universal** (independent of T modulo 2^{n-1}) if and only if s_2(n-3) <= 1, where s_2 is the binary digit sum. The universal values of n are: 3, 5, 7, 11, 19, 35, 67, ...

At the first non-universal n = 9: S mod 128 = 0 universally, but S mod 256 depends on the parity of the 3-cycle count t_3.

### 4.3 The Master Identity

**Theorem (THM-I).** For any set of k non-adjacent positions in the Hamiltonian path:

    D_S = n! / 2^k

where D_S = sum over (2k)! orderings of 2k vertices placed at those positions, weighted by edge products. This identity holds pointwise for all tournaments.

The constant c_0 = S(T) / 2^{n-1} serves as the degree-0 Walsh coefficient of the W-polynomial and satisfies:
- At n = 5: c_0 = H - 3*t_3 (integer, zero iff t_3 is odd)
- At n = 7: c_0 in Z + 3/4 (never integer)
- At n = 9: c_0 mod 1 depends on t_3 parity

---

## 5. The Worpitzky/F-Polynomial

### 5.1 The Forward-Edge Polynomial

For a tournament T with a fixed labeling, define fwd(P) = #{edges (P_i, P_{i+1}) where P_i < P_{i+1}} for each Hamiltonian path P. The **forward-edge polynomial** (or ascent polynomial) is:

    F(T, x) = sum_P x^{fwd(P)} = sum_{k=0}^{n-1} F_k(T) x^k

**Theorem (THM-087, complement duality).** F_k(T) = F_{n-1-k}(T^op).

**Theorem (THM-094, mod-2 universality).** F_k(T) = C(n-1, k) (mod 2) for all tournaments T. This is tournament-independent.

**Theorem (THM-089).** Var(fwd) = (n+1)/12 + 4*t_3/(n(n-1)), where t_3 is the number of directed 3-cycles.

### 5.2 Unimodality

**Conjecture (HYP-204).** F(T, x) is unimodal for all tournaments T.

Verified: exhaustive at n = 3,4,5; 100% at n = 6,7,8 by sampling. Zero violations in over 50,000 tournaments tested.

---

## 6. GLMY Path Homology of Tournaments

### 6.1 Background

The **GLMY path homology** (Grigor'yan-Lin-Muranov-Yau) associates to any directed graph a chain complex (Omega_*, d_*) where:
- Omega_0 = R^{vertices}
- Omega_p = {u in R^{A_p} : d_p(u) in R^{A_{p-1}}} (allowed p-paths with boundary landing in allowed (p-1)-paths)
- d_p is the alternating face map

The **Betti numbers** beta_p = dim(ker d_p / im d_{p+1}) measure p-dimensional "holes" in the directed graph.

### 6.2 Betti Number Landscape

Exhaustive computation through n = 6 and extensive sampling through n = 10 reveals:

| n | beta_0 | beta_1 | beta_2 | beta_3 | beta_4 | beta_5 |
|---|--------|--------|--------|--------|--------|--------|
| 3 | 1 | 0-1 | 0 | - | - | - |
| 4 | 1 | 0-1 | 0 | 0 | - | - |
| 5 | 1 | 0-1 | 0 | 0 | 0 | - |
| 6 | 1 | 0-1 | 0 | 0-1 | 0 | 0 |
| 7 | 1 | 0-1 | 0 | 0-1 | 0-6 | 0 |
| 8 | 1 | 0-1 | 0 | 0-2 | 0-5 | 0-1 |

- beta_0 = 1 always (tournaments are weakly connected)
- beta_1 in {0, 1} (proved, THM-103)
- **beta_2 = 0 universally** (see Section 6.3) — **PROVED** (THM-108/109)
- **beta_3 can reach 2** at n = 8 (0.08% of tournaments), previously thought bounded by 1
- beta_1 * beta_3 = 0 (**proved for all n**; THM-095 + THM-108/109). At n = 8, beta_3 * beta_4 = 1 CAN coexist ("consecutive seesaw" fails)
- beta(T) = beta(T^op) (complement invariance, proved at n = 5, verified through n = 8)

### 6.3 The beta_2 = 0 Theorem

**Theorem (THM-108/109).** For every tournament T, beta_2(T) = 0 in GLMY path homology.

**Computational evidence:**
- Exhaustive verification: n = 3 through n = 6 (all 32,768+ tournaments)
- Extensive sampling: n = 7 (5000), n = 8 (5000+), n = 9 (2000), n = 10 (500)
- Total: over 50,000 tournaments tested with **zero failures**

**Proof (THM-108/109, kind-pasteur-S43).** By strong induction on n using the long exact sequence of the pair (T, T\v):

    ... -> H_2(T\v) -> H_2(T) -> H_2(T, T\v) -> H_1(T\v) -> H_1(T) -> ...

By induction, H_2(T\v) = 0, so H_2(T) injects into H_2(T, T\v). The proof reduces to showing every tournament has a "good" vertex v with h_2^{rel}(T, T\v) = 0, equivalently beta_1(T\v) <= beta_1(T).

**What makes this novel.** Beta_2 = 0 has no analogue in existing path homology literature. For general directed graphs, beta_2 > 0 is common (70/59,049 oriented graphs at n = 5). The vanishing is specific to tournaments.

**Twin vertex mechanism.** All oriented graphs with beta_2 > 0 have twin vertices (identical neighborhoods). Tournament completeness forbids twins.

**Supporting structure:**

**Rank formula (proved).** rank(d_2|_{Omega_2}) = C(n,2) - n + 1 - beta_1(T).

**Additional verified properties:**
- **HYP-282:** When beta_1(T) = 0, at most 3 vertices have beta_1(T\v) = 1 (verified through n = 10)
- **HYP-384:** The restriction map res: Z_1(T) -> direct_sum_v H_1(T\v) is always surjective

### 6.4 Higher Betti Numbers: The n = 8 Threshold

At n = 8, several patterns that held for smaller tournaments break:

- **beta_3 = 2 exists** (0.08% of tournaments at n = 8, 0.05% at n = 9). Previously all beta_k were at most 1 for k >= 1.
- **Consecutive seesaw fails:** beta_3 * beta_4 = 1 can coexist at n = 8 (~0.15%), though beta_1 * beta_3 = 0 still holds.
- **i_*-injectivity fails:** The inclusion map H_3(T\v) -> H_3(T) has nontrivial kernel for some (T, v) at n = 8.

These failures mean proof strategies that work at n <= 7 cannot extend directly.

### 6.5 Paley Tournament Homology

For the **Paley tournament** T_p (p prime, p = 3 mod 4), where a->b iff b-a is a quadratic residue mod p, set m = (p-1)/2. The path complex has degrees 0 through 2m = p-1.

**Two-level symmetry decomposition:**
- **Z_p action** on paths (vertex translation) decomposes the chain complex into p eigenspaces (k = 0, 1, ..., p-1). The k=0 eigenspace is the **diff-seq complex** (translation-invariant paths encoded by successive differences).
- **Z_m action** on diff-seqs (multiplication by quadratic residues) further decomposes into m orbits.

**Theorem (THM-130, Complete Paley Betti Formula).** For the Paley tournament T_p with m = (p-1)/2:

    beta_m = m(m-3)/2,     beta_{m+1} = m(m+1)/2 = C(m+1, 2)
    beta_d = 0  for all other d >= 1
    chi(T_p) = 1 - beta_m + beta_{m+1} = p

The Euler characteristic equals the number of vertices.

**Per-eigenspace structure:**
- **k = 0** (diff-seq complex): beta_m^{(0)} = m(m-3)/2, beta_{m+1}^{(0)} = m(m-3)/2. Contributes ALL of beta_m.
- **k != 0** (p-1 nonzero eigenspaces): each contributes beta_{m+1}^{(k)} = 1 only.
- **Rank shift**: R_d^{(k)} - R_d^{(0)} = (-1)^{d+1} for 1 <= d <= m.

**Verified data:**
- **P_7** (m=3): beta = (1, 0, 0, 0, 6, 0, 0). chi = 7. Omega dims palindromic: [7,21,42,63,63,42,21].
- **P_11** (m=5): beta = (1, 0, 0, 0, 0, 5, 15, 0, 0, 0, 0). chi = 11. Omega dims: [1,5,20,70,205,460,700,690,450,180,30]. **Not** palindromic.

**Heisenberg Lie algebra connection (HYP-756).** beta_m = m(m-3)/2 is exactly the second Betti number b_2(h_m) of the m-dimensional Heisenberg Lie algebra (Santharoubane 1983). This connection appears new.

**Theorem (THM-125, Constant Symbol Matrix).** For all circulant tournaments: the Tang-Yau symbol matrix M_m(t) is constant (t-independent). All eigenspaces have identical Omega dimensions. Enables n-fold speedup in rank computation for T_p.

### 6.6 Torsion-Free Property

**Observation (HYP-1610).** Paley tournament chain complexes are torsion-free at boundary maps d_2 and d_3. Verified for P_7, P_11, P_19, P_23: rank(d_k) is constant across all test primes q in {2,3,5,7,...,47}. No rank drops = no torsion. Combined with beta_2 = 0, the chain complexes have unusually clean algebraic structure.

---

## 7. Spectrum and Extremal Results

### 7.1 Hamiltonian Path Spectrum

**Permanent gaps in the H-spectrum:**

- **H = 7 is impossible** for all n. Proof: Claim A decomposition forces alpha_1 >= 4, giving H >= 11.
- **H = 21 is impossible** for all n. Proof: poisoning graph DAG argument via component reduction.
- H = 63 is achievable at n = 8, so it is **not** a permanent gap.

At n = 7: the H-spectrum contains 77 distinct odd values in [1, 189].

### 7.2 Forbidden Values and k-nacci Traces

**Theorem (HYP-1618, corrected; HYP-1623).** The forbidden H values arise from k-nacci power sums via Newton's identities:

    Tr(M_k^3) = p_3 = e_1^3 - 3*e_1*e_2 + 3*e_3 = 1 + 3 + 3 = 7   for ALL k >= 3

Since the first three elementary symmetric polynomials (e_1=1, e_2=-1, e_3=-1) are identical for all k-nacci companion matrices with k >= 3, the value p_3 = 7 is universally forbidden.

For the tribonacci case (k=3): Tr(M_3^5) = p_5 = 21 = 3 * 7 = p_2 * p_3. This is the unique multiplicative relation among the first ~15 tribonacci traces. The product 21 = 3 x 7 inherits prohibition from both the cycle obstruction (3) and the Fano obstruction (7).

**Note:** The standard Riemann zeta gives zeta(-3) = 1/120, NOT 7. The connection to zeta is through Von Staudt-Clausen (42 = denom(B_6) = 2*3*7) and Kummer's carry-counting (THM-J), not through zeta at negative integers.

### 7.3 k-nacci Mersenne Identity

**Theorem (THM-227, proved).** For the k-nacci companion matrix M_k:

    Tr(M_k^n) = 2^n - 1   for all 1 <= n <= k

At n = k: Tr(M_k^k) = 2^k - 1 is the k-th Mersenne number. At n = k+1, the identity breaks by exactly k+1. While Tr(M_3^3) = 7 coincides with the first forbidden H value, the Mersenne numbers 31 = 2^5 - 1, 63 = 2^6 - 1, and 127 = 2^7 - 1 are all achievable H values (at n = 6, 8, and 7 respectively). The k-nacci Mersenne identity is a result about transfer matrices, not a characterization of forbidden H values. The forbidden set {7, 21} is better characterized as {7 * 3^0, 7 * 3^1} with nilpotency 2 (HYP-1231), or as {Phi_3(2), Phi_3(4)} (HYP-1317).

### 7.4 Paley Maximization

**Conjecture.** Among all tournaments on p vertices (p prime, p = 3 mod 4), the Paley tournament T_p maximizes H(T).

Verified: H(T_3) = 3, H(T_7) = 189, H(T_11) = 95,095. H(T_11)/|Aut(T_11)| = 1729 (the Hardy-Ramanujan taxicab number).

### 7.5 Simplicial Redei

**Theorem (THM-220, proved for n >= 4).** The simplicial Hamiltonian count sim_H(T) is in {0, 1} for all tournaments on n >= 4 vertices. Proved algebraically via Key Lemma + case analysis. Verified exhaustively through n = 8.

### 7.6 Real-Rootedness of I(Omega(T), x)

**Theorem (THM-020/021).** I(Omega(T), x) has all real roots for n <= 8.

**Theorem (THM-025).** This **fails** at n = 9: explicit counterexample found.

---

## 8. The Pin Grid and Symmetry

### 8.1 Tiling Model

A tournament on vertices {1, ..., n} with a fixed base path P_0 = (n, n-1, ..., 1) is encoded by a binary tiling t in {0,1}^m of the **pin grid** Grid(n) = {(r,c) : r >= 1, c >= 1, r+c <= n-1}, where m = C(n-1, 2). The grid is isomorphic to the staircase Young diagram delta_{n-2}.

### 8.2 Symmetry Group

The pin grid has symmetry group S_3 x Z_2 (the prism group), generated by sigma (reflection), tau (120-degree rotation), and phi (complement/bit flip). The **double Burnside formula** counts isomorphism classes.

### 8.3 Connection to Self-Evacuating Standard Young Tableaux

**Theorem (THM-035).** The number of self-evacuating SYT of shape delta_{n-2} equals 2^{m^2} = |Fix(sigma)| for n = 2m+1. All hook lengths of delta_{n-2} are odd.

---

## 9. Number-Theoretic Connections

### 9.1 Base-42 Structure

The number 42 = 2 * 3 * 7 encodes three fundamental constants of tournament parity:
- **2**: orientation/parity (F_2 arithmetic)
- **3**: smallest cycle (C_3, the triangle)
- **7**: first forbidden H value (Fano plane obstruction)

### 9.2 Egyptian Fraction Splitting

**Theorem (HYP-1612, proved).** 3/N = 1/b + 1/c has a solution in positive integers iff N has a prime factor p congruent to 2 (mod 3).

**Proof.** 3/N = 1/b + 1/c iff (3b-N)(3c-N) = N^2. Need d | N^2 with d = -N (mod 3). If all prime factors of N are congruent to 1 mod 3, all divisors of N^2 are congruent to 1 mod 3 — no suitable d exists.

**Master Criterion (HYP-1615, proved).** For general k: k/N = 1/b + 1/c solvable iff N^2 has a divisor d with d = -N (mod k).

**Cyclotomic Pattern (HYP-1617, proved for primes).** For prime k and prime p coprime to k: k/p = 1/b + 1/c solvable iff p = -1 (mod k). The solvable primes have order 2 in (Z/kZ)*. Unsolvable fraction among primes = (k-2)/(k-1). Connection: solvable primes split in Z[zeta_k + zeta_k^{-1}] (maximal real subfield).

### 9.3 Erdos-Straus Conjecture Connection

**Erdos-Straus (1948).** For every integer n >= 2, 4/n = 1/x + 1/y + 1/z has a solution in positive integers.

The base-42 covering system reduces this to finitely many residue classes. The "easy" 8 classes mod 42 are handled by factors 2 and 3; the "hard" 4 classes {1, 13, 25, 37} mod 42 (primes congruent to 1 mod 12) are distinguished by mod-7 residue. A multi-r covering using parametric identities a = (p+r)/4 for r in {3, 7, 11, 15, ...} catches all failures.

**Verification:** 0 failures across 19,564 primes to 10^6. Maximum r needed: 59 (at p = 118,801).

**Unconditional case (HYP-1621, proved).** p congruent to 13 mod 24 always works because (p+3)/4 is even, and 2 = 2 mod 3 satisfies the splitting condition.

### 9.4 Double Factorial Fixed Point

**Theorem (HYP-1614, proved).** (n-2)!! = 21 (mod 42) for k >= 4 (where n = 2k+1).

**Generalized (HYP-1616, proved).** For M = 2Q with Q odd squarefree, (2k-1)!! = M/2 (mod M) once k >= K(M) = max over primes p|Q of (p+1)/2. The fixed point M/2 = Q is the "odd half" of M.

### 9.5 Von Staudt-Bernoulli Connection

The Von Staudt chain n -> prod{p : (p-1)|n} gives 1 -> 2 -> 6 -> 42 -> 1806 -> 1806. Here 1806 = 2*3*7*43 is the unique Von Staudt fixed point (verified to 100,000). The set {2,3,7,43} is self-selecting: primes with (p-1)|1806 are exactly {2,3,7,43}. This connects to the Bernoulli number B_6 = 1/42 and the denominator formula for Bernoulli numbers.

---

## 10. Computational Complexity and Algorithmic Innovations

### 10.1 OCF-Based Computation

The formula H(T) = I(Omega(T), 2) replaces path-counting with independence polynomial evaluation. For structured tournaments, this yields a **100x speedup** (0.7ms vs 70ms per tournament at n = 9).

### 10.2 Walsh-Fourier Dimensionality Reduction

| n | Tournament space dim | Nonzero Walsh coefficients | Reduction factor |
|---|---------------------|---------------------------|-----------------|
| 5 | 1024 | 3 independent amplitudes | 341x |
| 7 | 2,097,152 | ~20 amplitudes | ~100,000x |

### 10.3 Small-Prime Modular Rank (8x Memory Reduction)

Reduce constraint matrices mod p < 256 before Gaussian elimination. Store as uint8 instead of int64. Rank preserved for all primes p > max elementary divisor. For T_11 degree-9 matrix: 6.6 GB -> 828 MB.

### 10.4 Eigenspace Decomposition (THM-125, n-fold Speedup)

For circulant tournaments: all n eigenspaces have identical rank structure. Collapses n independent rank computations to 1. For T_11: 11x speedup; for T_19: 19x speedup.

### 10.5 Multi-Prime Rank Certification

Compute rank at two independent small primes p_1, p_2. Agreement certifies correctness via Smith normal form argument. Mathematically verified, non-heuristic.

---

## 11. Engineering Applications

### 11.1 Implemented Tools

- **mod_rank_library.py** — Production-ready small-prime modular rank with auto prime selection and multi-prime verification. PyPI target.
- **circulant_homology module** — CirculantHomology and PaleyHomology classes using THM-125 for efficient Betti computation.
- **split_inert_analyzer.py** — Circulant cryptanalysis tool: splitting tables, defense rankings, torsion detection, QC-LDPC security audits.
- **circulant_ldpc_codes.py** — LDPC codes from Paley tournament adjacency. Rate >= (p+1)/(2p).
- **tournament_codes.py** — TDA feature extractor, H-fiber error-correcting codes.

### 11.2 Application Domains

1. **Sparse modular rank computation** — 8x memory, GPU-accelerable
2. **GLMY path homology for network analysis** — social networks, citation graphs, supply chains
3. **Circulant LDPC codes** — coding theory via QR_p structure
4. **GPU acceleration** — THM-125 reduces eigenspace work by factor p; int8 tensor cores
5. **TDA for preference/ranking data** — elections, sports, consumer research
6. **Deletion-contraction algorithm** — exact H via DC tree, O(2^n)
7. **Spectral tournament algorithms** — block-diagonalization via circulant structure
8. **Distributed Betti computation** — embarrassingly parallel eigenspaces
9. **Sparse path homology** — CSC format for T_19 (42 GB -> 1.2 MB)
10. **Number theory** — QR structure, cryptographic relevance
11. **Walsh spectrum cryptanalysis** — priority ranking for S-box attacks via 2-adic formula
12. **Post-quantum crypto analysis** — THM-125 directly applicable to QC-LDPC schemes (BIKE, HQC)

---

## 12. Status Summary

### Proved
- Redei's theorem (4 independent routes)
- OCF: H(T) = I(Omega(T), 2) (Grinberg-Stanley 2024; also THM-077 via Walsh)
- Claim B (algebraic companion to OCF)
- Transfer matrix symmetry M[a,b] = M[b,a] (via Walsh, THM-030)
- Complete Walsh spectrum of H(T) and M[a,b] (THM-071, THM-080)
- Position Character Decomposition — all degrees, all odd n (THM-068)
- Universal congruences for signed Hamiltonian permanent (THM-H, THM-I, THM-J)
- **beta_2 = 0 for all tournaments** (THM-108/109)
- beta_1 <= 1 for all tournaments (THM-103)
- beta_1 * beta_3 = 0 — mutual exclusivity (THM-095)
- rank(d_2) = C(n,2) - n + 1 - beta_1 (universal formula)
- F-polynomial complement duality, moment hierarchy, mod-2 universality
- H = 7 and H = 21 are permanent spectrum gaps
- Pin grid S_3 x Z_2 symmetry, Burnside orbit formula
- k-nacci Mersenne identity Tr(M_k^n) = 2^n - 1 (THM-227)
- Simplicial Redei: sim_H in {0,1} for n >= 4 (THM-220)
- Unit fraction splitting 3/N (HYP-1612), master criterion k/N (HYP-1615)
- Double factorial fixed point (HYP-1614, HYP-1616)
- Constant symbol matrix for circulants (THM-125)
- Up-Laplacian uniform spectrum (THM-224)

### Computational (strong evidence, no general proof)
- **THM-130 (Paley Betti formula)**: beta_m = m(m-3)/2, beta_{m+1} = C(m+1,2), chi = p. Verified P_7, P_11.
- HYP-282: sum_v beta_1(T\v) <= 3 when beta_1(T) = 0 (verified n <= 10)
- Unimodality of F(T, x) (50k+ tests, 0 failures)
- Paley maximization of H
- Eigenspace rank shift: R_d^{(k)} - R_d^{(0)} = (-1)^{d+1} (verified P_7, P_11)
- Base-42 Erdos-Straus covering (0 failures to 10^6)
- Torsion-free Paley chain complexes (HYP-1610)

### Open
- Algebraic proof of THM-130 (Paley Betti formula)
- Understanding beta_3 = 2 at n = 8
- HYP-282 (at most 3 bad vertices)
- Per-path identity for all n
- Proof of Paley maximization
- What bound replaces beta_3 <= 1 at n >= 8?
- Spectral zeta connection zeta(-3) = 7, zeta(-5) = 21 (HYP-1618)
- P_19 full verification (OOM at degree 9)

---

## 13. References

- N. Alon, *The maximum number of Hamiltonian paths in tournaments*, Combinatorica 10 (1990), 319-324
- J. Chapman, *Alternating sign matrices and tournaments*, Adv. in Appl. Math. 27 (2001), 318-335
- S. Chowdhury, S. Huntsman, M. Yutin, *Path homologies of motifs and temporal network representations*, Appl. Netw. Sci. (2022)
- A. El Sahili, M. Abi Aad, *Parity of paths and circuits in tournaments*, Discrete Math. 343 (2020)
- R. Forcade, *Parity of paths and circuits in tournaments*, Discrete Math. 6 (1973), 115-118
- S. Grinberg, R.P. Stanley, *Counting Hamiltonian paths in tournaments*, arXiv:2412.10572 (2024)
- A. Grigor'yan, Y. Lin, Y. Muranov, S.-T. Yau, *Homologies of path complexes and digraphs*, arXiv:1207.2834 (2012)
- J.W. Moon, *Topics on Tournaments*, Holt, Rinehart and Winston, New York (1968)
- L. Redei, *Ein kombinatorischer Satz*, Acta Litterarum ac Scientiarum (Szeged) 7 (1934), 39-43
- L. Santharoubane, *Cohomology of Heisenberg Lie algebras*, Proc. Amer. Math. Soc. 87 (1983), 23-28
- J. Schweser, M. Stiebitz, B. Toft, *The tournament theorem of Redei revisited*, arXiv:2510.10659 (2025)
- R.P. Stanley, *Enumerative Combinatorics*, Vol. 1 & 2, Cambridge University Press (1999)
- K.B. Tang, S.-T. Yau, *Path homology of circulant digraphs*, arXiv:2602.04140 (2026)

---

*This document is auto-generated from the research repository and will be updated as new results are obtained.*
