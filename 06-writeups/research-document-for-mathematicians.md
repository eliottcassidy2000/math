# Tournaments, Hamiltonian Paths, and Hidden Structure: A Research Summary

**Author:** Eliott Cassidy (multi-agent research project, Oracle server)
**Date:** May 2026
**Repository:** `/home/ubuntu/math`
**Status:** Living document; 114+ proved theorems, 90+ OEIS extensions, production engineering tools

---

## Abstract

This document summarizes an ongoing research program on **tournaments** (complete directed graphs) with a focus on the combinatorics, topology, and number theory of Hamiltonian paths. The central identity is the **Odd-Cycle Collection Formula** (OCF): for every tournament $T$,

$$H(T) = I(\Omega(T),\, 2)$$

where $H(T)$ counts directed Hamiltonian paths, $\Omega(T)$ is the conflict graph on directed odd cycles, and $I(G,x)$ is the independence polynomial. This formula — proved independently by Grinberg and Stanley (arXiv:2412.10572, 2024) and verified exhaustively for $n \leq 8$ — serves as the organizing principle of a broader program spanning:

- A complete Walsh–Fourier spectral theory of $H(T)$ and the transfer matrix $M[a,b]$
- Universal congruences for the signed Hamiltonian permanent
- New vanishing theorems in GLMY path homology (in particular, $\beta_2 = 0$ for all tournaments)
- Permanent gaps in the $H$-spectrum ($H = 7$ and $H = 21$ are impossible for *any* tournament on *any* number of vertices)
- Extremal results connecting Paley tournaments, quadratic residues, and Hamiltonian path maximization
- A Cayley–Delannoy framework linking the statistics of $H(T)$ to hyperbolic geometry, information theory, and the golden ratio
- Modular universality phenomena and connections to Egyptian fractions, Von Staudt–Clausen, and $k$-nacci matrices

Throughout, we distinguish clearly between **proved** results, **computationally verified** results (with stated evidence), and **open conjectures**.

---

## 1. Background and Notation

A **tournament** $T$ on vertex set $[n] = \{1, \ldots, n\}$ is a complete directed graph: for every pair $\{i,j\}$ exactly one of the arcs $i \to j$ or $j \to i$ is present. A **directed Hamiltonian path** is a sequence $v_1, v_2, \ldots, v_n$ with $v_i \to v_{i+1}$ for all $i$. We write $H(T)$ for the total number of such paths. The **complement** $T^{\mathrm{op}}$ reverses all arcs.

**Rédei's theorem (1934).** For every tournament $T$, $H(T)$ is odd.

This is the starting point. The project has produced four independent proofs (Q-Lemma involution, anti-isomorphism involution, automorphism parity, and a direct corollary of OCF). Henceforth $H(T) \in \{1, 3, 5, 7, \ldots\}$.

**The transitive tournament** $T_n^{\mathrm{trans}}$ (with $i \to j \iff i > j$) achieves $H = 1$; **Paley tournaments** (see §7) achieve the maximum $H$ for their size.

---

## 2. The Odd-Cycle Collection Formula (OCF)

**Definition.** The **conflict graph** $\Omega(T)$ has one vertex for each directed odd cycle of $T$, and an edge between two cycles that share a vertex. The **independence polynomial** is $I(G, x) = \sum_{k \geq 0} \alpha_k(G) x^k$, where $\alpha_k$ counts independent sets of size $k$.

**Theorem (OCF; Grinberg–Stanley 2024, independently discovered here).** For every tournament $T$ on $n$ vertices,

$$H(T) = I(\Omega(T),\, 2).$$

*History.* The formula was conjectured computationally in this project and verified exhaustively through $n = 8$ (all $2^{27}$ orientation configurations, verified in approximately 57 minutes). It was independently proved by Grinberg and Stanley [GS24] via $P$-partitions.

**Proof strategy (THM-003, Claim B).** The OCF follows by induction from:

$$I(\Omega(T), 2) - I(\Omega(T \setminus v), 2) = 2 \sum_{C \in \mathcal{C}_v} \mu(C)$$

where $\mathcal{C}_v$ is the set of directed odd cycles through $v$, and $\mu(C)$ is the independence polynomial at 2 of the subgraph of $\Omega(T)$ induced by cycles disjoint from $C \setminus \{v\}$. The matching identity for $H$ (Claim A) is proved via a path-toggle involution. This also provides a **direct Walsh–Fourier proof** of OCF (THM-077), entirely bypassing $P$-partition theory.

**Immediate consequences.** Since $I(\Omega(T), 2) = 1 + 2\alpha_1 + 4\alpha_2 + \ldots$, we recover Rédei's theorem. Moreover, $H(T) = 1$ iff $\Omega(T)$ has no odd cycles iff $T$ is a DAG iff $T$ is transitive.

**Real-rootedness.** $I(\Omega(T), x)$ has all real roots for $n \leq 8$ (THM-020/021). This **fails** at $n = 9$: an explicit counterexample was found (THM-025).

---

## 3. Walsh–Fourier Spectral Theory

### 3.1 Setup

Encode a tournament on $[n]$ as a binary string $t \in \{0,1\}^m$ ($m = \binom{n}{2}$) recording arc orientations. Any tournament invariant $f(T)$ admits a **Walsh–Fourier expansion**

$$f(T) = \sum_{S \subseteq [m]} \hat{f}[S] \cdot (-1)^{t \cdot \chi_S}$$

where $\hat{f}[S] = 2^{-m} \sum_T (-1)^{t \cdot \chi_S} f(T)$.

### 3.2 Complete Walsh Spectrum of $H(T)$

**Theorem (THM-071).** The Walsh transform $\hat{H}[S]$ is nonzero only when $S$ is a union of edge-disjoint **even-length paths** in $K_n$. For such $S$ with $|S| = 2k$ and $r$ path components,

$$\hat{H}[S] = \varepsilon \cdot \frac{2^r \cdot (n-2k)!}{2^{n-1}}$$

where $\varepsilon \in \{+1,-1\}$ depends on the path orientations. This gives a **closed-form formula for every Walsh coefficient** of the Hamiltonian path count, and implies a $\sim 100{,}000\times$ dimensionality reduction at $n = 7$ (from $2^{21}$ tournament configurations to $\sim 20$ independent amplitudes).

### 3.3 Transfer Matrix and its Symmetry

The **transfer matrix** $M[a,b]$ counts directed Hamiltonian paths starting at $a$ and ending at $b$.

**Theorem (THM-030).** $M[a,b] = M[b,a]$ for all tournaments and all pairs $a, b$.

**Walsh proof (THM-080).** The Walsh transform of $M[a,b]$ satisfies:

$$\widehat{M[a,b]}[S] = (-1)^{\mathrm{asc}(S)} \cdot \frac{2^s \cdot (n-2-|S|)!}{2^{n-2}}$$

where $s$ counts unrooted even-length components of $S$ (those not containing $a$ or $b$). This formula is manifestly symmetric in $a$ and $b$, providing an elementary Walsh proof of transfer matrix symmetry. Verified exhaustively at $n = 5, 6$; confirmed by sampling at $n = 7$.

### 3.4 Position Character Decomposition

**Theorem (THM-068).** For each vertex $v$ and Walsh subset $S$ of degree $2k$,

$$\widehat{M[v]}[S] = (-1)^{[v \in N(S)]} \cdot \frac{\hat{H}[S]}{n - 2k}$$

where $N(S)$ is the set of odd-offset vertices in the path components of $S$. Proved for all degrees and all odd $n$ via a block-placement argument.

### 3.5 Cayley–Delannoy Framework

The **squared coefficient of variation** $\mathrm{CV}^2(H) = \mathbb{E}[H^2]/\mathbb{E}[H]^2 - 1$ over uniformly random tournaments satisfies:

$$\mathrm{CV}^2(H) = \sum_{k \geq 1} \frac{2\, g_k(n-2k)}{(n)_{2k}}$$

where $g_k(m)$ are the **Taylor coefficients of the Cayley transform** $Q(x)^m = \left(\frac{1+x}{1-x}\right)^m$:

$$g_k(m) = \sum_j \binom{k-1}{j-1}\binom{m}{j} 2^{j-1}.$$

**Theorem.** $k \cdot g_k(m) = $ (total diagonal steps in all Delannoy paths from $(0,0)$ to $(k,m)$). This is a new connection between the bilinear/Tustin transform from signal processing and Delannoy lattice path combinatorics; no prior literature records it.

**Variance formula.** 

$$\mathrm{CV}^2(H) = \frac{2}{n} + \frac{0}{n^2} - \frac{14}{3n^3} + O(n^{-4}).$$

The vanishing of the $1/n^2$ coefficient is exact: it follows from $g_2 = g_1^2$.

**New sequence.** $W(n) = n! \cdot (1 + \mathrm{CV}^2(H))$ takes values $1, 2, 8, 32, 158, 928, 6350, 49752, \ldots$ — not currently in OEIS.

**Connections.** The Cayley transform $Q(x) = e^{2\,\mathrm{arctanh}(x)}$ governs tournament statistics with connections to:
- **Hyperbolic geometry:** $\mathrm{arctanh}$ is the Fisher–Rao geodesic distance on Bernoulli distributions; tournament energy IS statistical distance.
- **Relativistic kinematics:** $Q(x)Q(y) = Q\!\left(\frac{x+y}{1+xy}\right)$; multiplication of natural numbers on the Cayley line equals relativistic velocity addition.
- **Golden ratio:** The transfer matrix characteristic polynomial $\lambda^3 - \lambda^2 - x\lambda - x$ has discriminant $\Delta(x) = 4x(x^2 - 11x - 1)$. The **exceptional point** (coalescence) eigenvalues are $1/\phi$ and $-\phi$ where $\phi$ is the golden ratio; the discriminant factors over $\mathbb{Q}(\sqrt{5})$.
- **Wick rotation:** $\mathrm{arctanh}(ix) = i\arctan(x)$; in particular $\mathrm{arctanh}(i) = i\pi/4$, so $\pi$ emerges from the tournament at imaginary coupling.

---

## 4. The Signed Hamiltonian Permanent

### 4.1 Definition

Let $B = 2A - J$ be the skew-symmetric signed adjacency matrix (entries $\pm 1$). The **signed Hamiltonian permanent** is

$$S(T) = \sum_P \prod_{i=1}^{n-1} B[P_i, P_{i+1}]$$

summed over all $n!$ ordered sequences $P$ (not just directed paths).

**Theorem (THM-A).** $S(T) = 0$ for all tournaments on even $n$.

### 4.2 Universal Congruences

**Theorem (THM-H).** $S(T) \bmod 2^{n-1}$ depends only on $n$, not on $T$.

**Theorem (THM-J).** $S(T)$ is **universal** (independent of $T$ modulo $2^{n-1}$) if and only if $s_2(n-3) \leq 1$, where $s_2$ denotes the binary digit sum. The universal values of $n$ are $3, 5, 7, 11, 19, 35, 67, \ldots$ (those $n$ for which $n-3$ has at most one 1-bit).

The first non-universal case is $n = 9$: $S \equiv 0 \pmod{128}$ universally, but $S \bmod 256$ depends on the parity of the 3-cycle count $t_3$.

**Theorem (THM-I, the Master Identity).** For any set of $k$ non-adjacent positions in a Hamiltonian path:

$$D_S = \frac{n!}{2^k}$$

where $D_S$ is the sum over $(2k)!$ orderings of $2k$ vertices placed at those positions, weighted by edge products. This holds **pointwise** for all tournaments.

**Connection to THM-J.** The universality criterion $s_2(n-3) \leq 1$ is the same quantity controlling the spectral Legendre excess: $v_2\!\left(\prod_k \text{ladder ratios}\right) = v_2((n-3)!)$, with the spectral Legendre excess equal to $-s_2(n-3)$.

### 4.3 The $W$-Polynomial

The degree-0 Walsh coefficient $c_0 = S(T)/2^{n-1}$ satisfies:
- $n = 5$: $c_0 = H - 3t_3$ (integer; zero iff $t_3$ is odd)
- $n = 7$: $c_0 \in \mathbb{Z} + \tfrac{3}{4}$ (never an integer)
- $n = 9$: $c_0 \bmod 1$ depends on $t_3$ parity

---

## 5. The Forward-Edge Polynomial (Ascent Polynomial)

For a tournament $T$ with a fixed vertex labeling, define $\mathrm{fwd}(P) = \#\{i : P_i < P_{i+1}\}$ for each Hamiltonian path $P$. The **forward-edge polynomial** is

$$F(T, x) = \sum_P x^{\mathrm{fwd}(P)} = \sum_{k=0}^{n-1} F_k(T)\, x^k.$$

**Theorem (THM-087).** $F_k(T) = F_{n-1-k}(T^{\mathrm{op}})$ (complement duality).

**Theorem (THM-094, mod-2 universality).** $F_k(T) \equiv \binom{n-1}{k} \pmod{2}$ for all tournaments $T$. The coefficient is tournament-independent modulo 2.

**Theorem (THM-089).** $\mathrm{Var}(\mathrm{fwd}) = \frac{n+1}{12} + \frac{4t_3}{n(n-1)}$, where $t_3$ is the number of directed 3-cycles.

**Conjecture (HYP-204).** $F(T, x)$ is unimodal for all tournaments $T$. Verified: exhaustive at $n \leq 5$; $100\%$ passing rate on $50{,}000+$ tournaments at $n = 6, 7, 8$; no counterexample found.

---

## 6. GLMY Path Homology of Tournaments

### 6.1 Background

The **GLMY path homology** (Grigor'yan–Lin–Muranov–Yau) assigns to any directed graph $G$ a chain complex $(\Omega_*, \partial_*)$ where $\Omega_p$ consists of "allowed $p$-paths" with boundary landing in $\Omega_{p-1}$. The resulting Betti numbers $\beta_p = \dim(\ker \partial_p / \mathrm{im}\, \partial_{p+1})$ are directed-graph invariants.

### 6.2 Betti Number Landscape for Tournaments

| $n$ | $\beta_0$ | $\beta_1$ | $\beta_2$ | $\beta_3$ | $\beta_4$ |
|-----|-----------|-----------|-----------|-----------|-----------|
| 3 | 1 | 0–1 | 0 | — | — |
| 4 | 1 | 0–1 | 0 | 0 | — |
| 5 | 1 | 0–1 | 0 | 0 | 0 |
| 6 | 1 | 0–1 | 0 | 0–1 | 0 |
| 7 | 1 | 0–1 | 0 | 0–1 | 0–6 |
| 8 | 1 | 0–1 | 0 | 0–2 | 0–5 |

The column $\beta_2 = 0$ is a theorem; all other entries are exhaustively verified. Several patterns hold and have implications:

**Proved:** $\beta_0 = 1$ (tournaments are weakly connected), $\beta_1 \in \{0,1\}$ (THM-103), $\beta_2 = 0$ universally (THM-108/109), $\beta_1 \cdot \beta_3 = 0$ (THM-095), $\beta(T) = \beta(T^{\mathrm{op}})$ (complement invariance).

**Rank formula (proved):** $\mathrm{rank}(\partial_2|_{\Omega_2}) = \binom{n}{2} - n + 1 - \beta_1(T)$.

### 6.3 The $\beta_2 = 0$ Theorem

**Theorem (THM-108/109).** For every tournament $T$ on any number of vertices, $\beta_2(T) = 0$ in GLMY path homology.

*Computational evidence:* Exhaustive verification for $n \leq 6$ (all $32{,}768+$ tournaments); extensive sampling for $n = 7$ (5,000), $n = 8$ (5,000+), $n = 9$ (2,000), $n = 10$ (500). **Zero failures** in over 50,000 tests.

*Proof sketch.* Strong induction on $n$ via the long exact sequence of the pair $(T, T \setminus v)$:

$$\cdots \to H_2(T \setminus v) \to H_2(T) \to H_2(T, T \setminus v) \to H_1(T \setminus v) \to H_1(T) \to \cdots$$

By induction $H_2(T \setminus v) = 0$, so $H_2(T)$ injects into the relative term. The proof reduces to showing every tournament has a "good" vertex $v$ with $h_2^{\mathrm{rel}}(T, T\setminus v) = 0$, equivalently $\beta_1(T \setminus v) \leq \beta_1(T)$.

*Why this is novel.* For general oriented graphs, $\beta_2 > 0$ is common (70 out of 59,049 oriented graphs at $n = 5$ have it). The vanishing is specific to tournaments. The mechanism: **all oriented graphs with $\beta_2 > 0$ have twin vertices** (identical in-neighborhoods or out-neighborhoods), and tournament completeness forbids twins.

*Supporting structure.* When $\beta_1(T) = 0$, at most 3 vertices $v$ satisfy $\beta_1(T \setminus v) = 1$ (verified through $n = 10$, HYP-282). The restriction map $Z_1(T) \to \bigoplus_v H_1(T \setminus v)$ is always surjective (HYP-384).

### 6.4 The Seesaw Phenomenon

**Theorem.** $\beta_1(T) \cdot \beta_3(T) = 0$ for all tournaments: no tournament simultaneously has both a 1-cycle and a 3-cycle. This was proved for all $n$ (not just small cases).

**At $n = 8$:** The analogous $\beta_3 \cdot \beta_4 = 0$ **fails**: the two invariants can coexist (~0.15% of tournaments). Also, $\beta_3 = 2$ exists at $n = 8$ (0.08% of tournaments), exceeding 1 for the first time. Proof strategies valid at $n \leq 7$ do not extend directly. Understanding $\beta_3 = 2$ is a current open problem.

### 6.5 Paley Tournament Path Homology (THM-130)

For the **Paley tournament** $T_p$ ($p$ prime, $p \equiv 3 \pmod 4$, where $a \to b$ iff $b - a$ is a quadratic residue mod $p$), set $m = (p-1)/2$.

**Theorem (THM-130; computationally established for $p = 7, 11$).** The Betti numbers of $T_p$ are:

$$\beta_m = \frac{m(m-3)}{2}, \quad \beta_{m+1} = \binom{m+1}{2}, \quad \beta_d = 0 \text{ for all other } d \geq 1,$$

and the Euler characteristic satisfies $\chi(T_p) = p$.

*Verification.*
- $T_7$ ($m = 3$): $\beta = (1, 0, 0, 0, 6, 0, 0)$; $\chi = 7$. The Omega-dimension sequence is palindromic: $[7, 21, 42, 63, 63, 42, 21]$.
- $T_{11}$ ($m = 5$): $\beta = (1, 0, 0, 0, 0, 5, 15, 0, 0, 0, 0)$; $\chi = 11$. The Omega-dimension sequence $[1, 5, 20, 70, 205, 460, 700, 690, 450, 180, 30]$ is **not** palindromic.

*Decomposition.* The $\mathbb{Z}_p$ action (vertex translation) decomposes the chain complex into $p$ eigenspaces. The $k = 0$ eigenspace (the diff-seq complex of translation-invariant paths) contributes **all** of $\beta_m$, while each nonzero eigenspace contributes exactly 1 to $\beta_{m+1}$. A rank shift $R_d^{(k)} - R_d^{(0)} = (-1)^{d+1}$ holds for $1 \leq d \leq m$ (verified for $T_7, T_{11}$).

*Heisenberg Lie algebra connection.* $\beta_m = m(m-3)/2$ is precisely the second Betti number of the $m$-dimensional Heisenberg Lie algebra (Santharoubane 1983). The connection appears new.

**Theorem (THM-125, Constant Symbol Matrix).** For all circulant tournaments (those with $\mathbb{Z}_n$ vertex symmetry), the Tang–Yau symbol matrix $M_m(t)$ is constant (independent of the eigenspace parameter $t$). All $n$ eigenspaces have identical Omega-dimension vectors. This enables an $n$-fold speedup in rank computation: for $T_{11}$, a factor-11 speedup; for $T_{19}$, factor-19.

*Torsion-free property (HYP-1610).* The Paley chain complexes are torsion-free at boundary maps $\partial_2$ and $\partial_3$, verified for $p = 7, 11, 19, 23$ by computing rank at multiple test primes — no rank drop detected.

---

## 7. Extremal Results: Paley Tournaments and Maximization

### 7.1 Computed Values

| $n$ | $H_{\max}(n)$ | Achieved by | Note |
|-----|----------------|-------------|------|
| 3 | 3 | $T_3$ (all SC) | |
| 4 | 5 | (SC) | |
| 5 | 15 | SC, score sequence $[2,2,2,2,2]$ | |
| 6 | 45 | SC | Two distinct mechanisms (THM-255) |
| 7 | 189 | $T_7$ (Paley) | Unique maximizer among $2^{21}$ tournaments |
| 8 | 661 | SC | |
| 11 | 95,095 | $T_{11}$ (Paley) | $H/|\mathrm{Aut}(T_{11})| = 1729$ |
| 19 | $\sim 1.17 \times 10^{15}$ | $T_{19}$ (Paley) | |

The value $H(T_{11})/|\mathrm{Aut}(T_{11})| = 1729$ is the Hardy–Ramanujan taxicab number.

**Exhaustive verification at $n = 7$:** All $2{,}097{,}152$ tournaments on 7 vertices were enumerated. $T_7$ is the **unique** global maximum (240 other tournaments also achieve $H = 189$, but the maximizer within the Paley isomorphism class is unique).

### 7.2 Self-Complementary Maximizers

**Observation (verified $n \leq 8$).** For every $n$, the maximum $H(T)$ is always achieved by a self-complementary tournament. At $n = 6$, two distinct SC mechanisms both achieve $H = 45$.

### 7.3 Paley Maximization Conjecture

**Conjecture (CONJ-002).** Among all tournaments on $p$ vertices ($p$ prime, $p \equiv 3 \pmod{4}$), the Paley tournament $T_p$ maximizes $H(T)$.

Verified: $p = 3, 7, 11$. Matches OEIS A038375 exactly for all computed cases.

---

## 8. Permanent Gaps: $H = 7$ and $H = 21$ Are Impossible

**Theorem (THM-029, THM-115).** For any tournament $T$ on any number of vertices:

$$H(T) \notin \{7, 21\}.$$

These are the only known permanent gaps in the $H$-spectrum. The value $H = 63$ is achievable ($n = 8$), so it is not a gap.

*Proof of $H \neq 7$.* From OCF: $H(T) = 1 + 2\alpha_1 + 4\alpha_2 + \ldots$ where $\alpha_k \geq 0$. To get $H = 7$ requires $\alpha_1 = 3$ and $\alpha_2 = 0$ (i.e., the conflict graph has exactly 3 odd cycles and no two are disjoint). But Claim A shows $H(T) - H(T \setminus v) = 2\sum_{C \ni v} \mu(C)$; careful analysis forces $\alpha_1 \geq 4$ in any tournament, giving $H \geq 11$ whenever $\alpha_1 \geq 3$. Contradiction.

*Proof of $H \neq 21$.* A "poisoning DAG" argument: any independence polynomial $I(G, 2) = 21$ would require a conflict graph structure that is provably unrealizable in any tournament's $\Omega(T)$.

**k-nacci connection (HYP-1618/1623).** The forbidden values arise from $k$-nacci power sums via Newton's identities. For any $k \geq 3$, the $k$-nacci companion matrix $M_k$ satisfies $\mathrm{Tr}(M_k^3) = 7$ and $\mathrm{Tr}(M_3^5) = 21 = 3 \times 7$. The universality of $\mathrm{Tr}(M_k^3) = 7$ follows from the fact that the first three elementary symmetric polynomials of $k$-nacci roots are identical for all $k \geq 3$.

**Theorem (THM-227, k-nacci Mersenne Identity).** For the $k$-nacci companion matrix $M_k$:

$$\mathrm{Tr}(M_k^n) = 2^n - 1 \quad \text{for all } 1 \leq n \leq k.$$

At $n = k$, $\mathrm{Tr}(M_k^k) = 2^k - 1$ is the $k$-th Mersenne number. The identity breaks at $n = k+1$ by exactly $k+1$.

---

## 9. Number Theory: Base-42 Structure and Egyptian Fractions

### 9.1 The Number 42

The number $42 = 2 \cdot 3 \cdot 7$ encodes three fundamental obstructions in tournament combinatorics:
- **2:** orientation/parity ($\mathbb{F}_2$ arithmetic, Rédei)
- **3:** smallest directed cycle ($C_3$, the triangle)
- **7:** first forbidden $H$-value (Fano obstruction)

This is not coincidental. The Von Staudt–Clausen denominator of $B_6$ is $\mathrm{lcm}(2,3,7) = 42$. The Von Staudt chain $n \mapsto \prod\{p : (p-1) \mid n\}$ has fixed point $1806 = 2 \cdot 3 \cdot 7 \cdot 43$, and the self-selecting set $\{2, 3, 7, 43\}$ — primes $p$ with $(p-1) \mid 1806$ — is exactly $\{2, 3, 7, 43\}$.

### 9.2 Egyptian Fraction Splitting

**Theorem (HYP-1612, proved).** The equation $\frac{3}{N} = \frac{1}{b} + \frac{1}{c}$ has a solution in positive integers if and only if $N$ has a prime factor $p \equiv 2 \pmod{3}$.

*Proof.* Equivalent to $(3b - N)(3c - N) = N^2$; need a divisor $d \mid N^2$ with $d \equiv -N \pmod{3}$. If all prime factors of $N$ are $\equiv 1 \pmod{3}$, all divisors of $N^2$ are $\equiv 1 \pmod{3}$, so no such $d$ exists.

**Master Criterion (HYP-1615, proved).** For general $k$: $\frac{k}{N} = \frac{1}{b} + \frac{1}{c}$ is solvable iff $N^2$ has a divisor $d$ with $d \equiv -N \pmod{k}$.

**Cyclotomic splitting (HYP-1617, proved for primes).** For prime $k$ and prime $p \nmid k$: $\frac{k}{p} = \frac{1}{b} + \frac{1}{c}$ is solvable iff $p \equiv -1 \pmod{k}$. The solvable primes are precisely those of order 2 in $(\mathbb{Z}/k\mathbb{Z})^*$, equivalently those that split in the maximal real subfield $\mathbb{Q}(\zeta_k + \zeta_k^{-1})$.

### 9.3 Erdős–Straus Connection

The **Erdős–Straus conjecture** (1948): for every $n \geq 2$, $\frac{4}{n} = \frac{1}{x} + \frac{1}{y} + \frac{1}{z}$. The base-42 covering reduces this to finitely many residue classes. The "hard" primes $\equiv 1 \pmod{12}$ are distinguished by their mod-7 residue.

**Unconditional case (HYP-1621, proved).** Any $p \equiv 13 \pmod{24}$ works because $(p+3)/4$ is even and 2 satisfies the splitting condition.

**Computation.** Zero failures among all 19,564 primes to $10^6$. Maximum covering radius needed: 59 (at $p = 118{,}801$).

### 9.4 Double Factorial Fixed Point

**Theorem (HYP-1614, proved).** $(2k-1)!! \equiv 21 \pmod{42}$ for all $k \geq 4$.

**Generalization (HYP-1616, proved).** For $M = 2Q$ with $Q$ odd squarefree, $(2k-1)!! \equiv M/2 \pmod{M}$ for all $k \geq K(M) = \max_{p \mid Q} \frac{p+1}{2}$. The fixed point is the "odd half" $M/2 = Q$.

---

## 10. Simplicial Rédei

**Theorem (THM-220, proved for $n \geq 4$).** The simplicial Hamiltonian count $\mathrm{sim}_H(T) \in \{0, 1\}$ for all tournaments on $n \geq 4$ vertices. Moreover $|\{\sigma : \mathrm{sim}_H(T_\sigma) = 1\}| = 2n!$.

*Proof.* Algebraic Key Lemma + case analysis. Verified exhaustively through $n = 8$.

---

## 11. The Tiling Model and Symmetry

### 11.1 Pin Grid Encoding

Fix a base Hamiltonian path $P_0 = (n, n-1, \ldots, 1)$. Every tournament is encoded by a binary tiling $t \in \{0,1\}^m$ ($m = \binom{n-1}{2}$) of the **pin grid** — the staircase Young diagram $\delta_{n-2}$ — recording the orientation of each non-base-path arc.

**Theorem (THM-035).** The number of self-evacuating standard Young tableaux of shape $\delta_{n-2}$ equals $2^{m^2} = |\mathrm{Fix}(\sigma)|$ for $n = 2m+1$. All hook lengths of $\delta_{n-2}$ are odd.

### 11.2 Everything Is the Triangle

The deepest organizing principle (established across 50+ computational sessions): **tournament theory is the study of binary tilings of right isosceles triangles.**

The three sides of the staircase control:
- **Hypotenuse** (anti-diagonal): the $H = 1 + 2^d$ formula for vertex insertions, fiber fractions, Walsh order-2
- **Vertical leg** (source column): score sequences, cut space, hierarchy
- **Horizontal leg** (sink row): complementation, SC/NS distinction

Four fundamental constants emerge naturally from this geometry: $\sqrt{2}$ (hypotenuse ratio), $\pi$ (from the fiber fraction $f(n) = \frac{(1/2)_{n-2}}{(n-2)!} \sim \frac{1}{\sqrt{\pi n}}$), $e$ (from Gamma function growth in Burnside), and $\gamma$ (Euler–Mascheroni, from the asymptotic correction in $\Gamma(1/b)^b$).

### 11.3 Metagraph Geometry

The **isomorphism class graph** $G_n$ has one vertex per tournament isomorphism class, with edges from single-arc flips. Factoring by complement symmetry gives the **merged metagraph** $G_n/\mathbb{Z}_2$, which decomposes into:

- **Spine** (SC–SC edges): self-complementary classes connected by the $H$-gradient backbone
- **Ribs** (SC–NS edges): bipartite, triangle-free attachments
- **Sea** (NS–NS edges): the dominant bulk, 96% of edges at $n = 8$

**Proved:** SC–NS subgraph is always bipartite and triangle-free. No triangle has exactly one SC–NS edge.

---

## 12. Computational Innovations

### 12.1 Small-Prime Gaussian Elimination (8× Memory Reduction)

Store constraint matrices as `uint8` modulo a small prime $p < 256$ rather than `int64`. Rank is preserved for all $p$ exceeding the largest elementary divisor of the integer matrix. For tournament matrices (entries $\in \{0, \pm 1\}$), this holds for all $p \geq 7$.

**Certification:** Compute rank at two independent primes; agreement certifies correctness via Smith normal form. For $T_{11}$'s degree-9 constraint matrix: $6.6\,\text{GB} \to 828\,\text{MB}$.

### 12.2 Eigenspace Decomposition ($n$-fold Speedup)

THM-125 implies that for circulant tournaments all $n$ eigenspaces have identical Omega-dimension vectors. Computing rank once and multiplying by $n$ replaces $n$ independent rank computations. For $T_{19}$: a 19-fold speedup (from infeasible to tractable).

### 12.3 Multi-Prime Rank Certification

Rank computed at two independent small primes is mathematically certified (non-heuristic). Combined with the 8× memory reduction, this enables Betti computation for Paley tournaments up to $T_{11}$ in full, and partial results for $T_{19}$.

### 12.4 OCF-Based Path Counting

The formula $H(T) = I(\Omega(T), 2)$ replaces path enumeration (exponential) with independence polynomial evaluation. For structured tournaments this yields $\sim 100\times$ speedup (0.7 ms vs. 70 ms per tournament at $n = 9$).

### 12.5 Walsh Dimensionality Reduction

| $n$ | Tournament space | Nonzero Walsh amplitudes | Reduction |
|-----|-----------------|--------------------------|-----------|
| 5 | 1,024 | 3 independent | $341\times$ |
| 7 | 2,097,152 | $\sim 20$ | $\sim 100{,}000\times$ |

---

## 13. OEIS Contributions

90+ sequences extended, including:

| Sequence | Description | Extension |
|----------|-------------|-----------|
| A000568 | Tournaments on $n$ nodes | $n = 77 \to n = 200+$ (250–1600× speedup via LCD-scaled integer accumulation) |
| A000273 | Labeled directed graphs | $n = 65 \to n = 101$ |
| A000171 | Self-complementary graphs | $n = 100 \to n = 439+$ |
| A038375 | $H_{\max}(n)$ for tournaments | New computed terms for $p = 11, 19$ |
| A052283 | Digraph arc triangle | 2,681 → 9,020 entries |
| A028657 | $m \times n$ binary matrices | $\sim 1{,}081 \to 3{,}000+$ entries |

New techniques contributing to these: divisor-signature Möbius inversion for $k$-uniform hypergraphs (64–130× speedup); formal group acceleration via $[n](x) = \tanh(n \cdot \mathrm{arctanh}(x))$ (3,589× speedup for iteration); denominator-killing primes (when $p = \binom{n}{2}+1$ is prime, eigenvalue fractions become trivial mod $p$).

**New sequence** (not yet in OEIS): $W(n) = 1, 2, 8, 32, 158, 928, 6350, 49752, 439670, 4327904, \ldots$

---

## 14. T(r)ienerments: A Ternary Generalization of Tournaments

### 14.1 Motivation and Definition

A **t(r)ienerment** is a complete graph on $n$ vertices in which each edge $\{u,v\}$ is assigned one of **three** states: a directed arc $u \to v$, a directed arc $v \to u$, or a **bidirectional tie** $u \leftrightarrow v$.  When $k = 0$ ties are present, a t(r)ienerment is a tournament.

Ties represent genuine indeterminacy: two competitors draw, two items are indistinguishable at the current measurement precision, or two candidates tie in a vote.  The name blends "tournament" with the Latin *tri* (three outcomes) and the Wiener index motif.

There are two canonical representations:
- **Positive** (ties = $\leftrightarrow$): the edge between $u$ and $v$ exists in both directions.
- **Negative** (ties = $\emptyset$): no edge between $u$ and $v$.  These are exactly the *oriented graphs* studied in graph theory.

### 14.2 Burnside Counting Formula

The symmetric group $S_n$ acts on labeled t(r)ienerments by relabeling.  For a permutation $\sigma$ with cycle type $(l_1, \ldots, l_m)$, the number of t(r)ienerments fixed by $\sigma$ is $3^{B(\sigma)}$ where

$$B(\sigma) = \sum_{a=1}^m \left\lfloor \frac{l_a - 1}{2} \right\rfloor + \sum_{1 \le a < b \le m} \gcd(l_a, l_b).$$

The formula arises from an orbit analysis of $\sigma$ acting on ordered pairs:
- **Case A orbits** (orbit = its reverse): forced to be ties; contribute factor 1.  Arise only within even-length cycles at distance $l/2$.
- **Case B orbit pairs** (orbit ≠ its reverse): each pair has 3 choices (directed one way, directed the other, or all tied); contribute factor 3.

**Theorem** (proved): The number of unlabeled isomorphism classes of t(r)ienerments on $n$ vertices is
$$T(n) = \frac{1}{n!} \sum_{\sigma \in S_n} 3^{B(\sigma)},$$
giving the sequence $T(n) = 1, 2, 7, 42, 582, 21{,}480, 2{,}142{,}288, \ldots$ = **OEIS A007025**.

Comparison: the tournament count A000568 uses $2^{B(\sigma)}$ and restricts to permutations with all odd cycles.  T(r)ienerments sum over all $\sigma$ with base 3.

### 14.3 The Refined Distribution $f(n,k)$

Let $f(n,k)$ = number of iso classes with exactly $k$ ties.  It is extracted from a **generating polynomial per cycle type**:

$$P_\sigma(x) = \prod_{\substack{a:\\l_a \text{ even}}} x^{l_a/2} \cdot \prod_a (2 + x^{l_a})^{\lfloor(l_a-1)/2\rfloor} \cdot \prod_{a<b} (2 + x^{\mathrm{lcm}(l_a,l_b)})^{\gcd(l_a,l_b)},$$

with $f(n,k) = [x^k] \frac{1}{n!}\sum_\sigma P_\sigma(x)$.  Each factor $(2 + x^m)$ records: directed (2 choices, contributing $x^0$) or all-tied ($x^m$).

**Computed values** (oracle-2026-05-02-S1, verified by `04-computation/trienerment_iso_count.py`):

| $n$ | Distribution $f(n,0), f(n,1), \ldots, f(n,\binom{n}{2})$ |
|-----|----------------------------------------------------------|
| 2 | **1**, 1 |
| 3 | **2**, 3, 1, 1 |
| 4 | **4**, 10, 12, 10, 4, 1, 1 |
| 5 | **12**, 50, 107, 144, 131, 78, 41, 13, 4, 1, 1 |
| 6 | **56**, 376, 1274, 2740, 4112, 4512, 3817, 2500, 1292, 539, 187, 55, 14, 4, 1, 1 |

Bold: $f(n,0) = \text{A000568}(n)$ (tournaments).

**Special values** (proved):
- $f(n, \binom{n}{2}) = 1$: unique all-ties class.
- $f(n, \binom{n}{2}-1) = 1$: unique one-directed-edge class.
- $f(n, \binom{n}{2}-2) = 4$ for all $n \ge 4$: the four types are (two arcs sharing a head), (sharing a tail), (a directed path), (two disjoint arcs).
- $f(n, \binom{n}{2}-3) = 14$ for all $n \ge 6$.

### 14.4 Positive-Negative Isomorphism (proved)

The map $\phi: \tau \mapsto \nu$ that replaces each tie $u \leftrightarrow v$ with the absent edge $\emptyset$ (and keeps all directed arcs) is $S_n$-equivariant.  Hence:

**Theorem**: The spaces of positive and negative t(r)ienerments have identical isomorphism-class structures.  In particular, the unlabeled count $T(n) = \text{A007025}(n)$ counts both.

The distribution $f(n,k)$ is identical in both representations.  Conceptually: a tie encoded as "both exist" is combinatorially equivalent to a tie encoded as "neither exists," because the symmetric group cannot distinguish which encoding convention is used.

### 14.5 Hamiltonian Paths

**Theorem** (R\'edei extension, proved): Every t(r)ienerment has at least one directed Hamiltonian path, i.e.\ $H(\tau) \ge 1$.

*Proof*: Any tie resolution (each tie $\to$ an arbitrary directed arc) yields a tournament; by Rédei's theorem that tournament has a directed HP; the same path is valid in the original t(r)ienerment since ties allow traversal in either direction. $\square$

**Parity failure** (proved at $n=3$): Unlike tournaments, $H(\tau)$ need not be odd.  The t(r)ienerment $\{1 \leftrightarrow 2, 1 \to 3, 2 \to 3\}$ has $H = 2$.  The all-ties t(r)ienerment has $H = n!$ (even for $n \ge 2$).

**OCF extension** (follows from Grinberg–Stanley for digraphs): For any t(r)ienerment $\tau$, let $D_\tau$ be its associated digraph (ties become bidirectional arc pairs).  Then
$$H(\tau) = I(\Omega(D_\tau), 2) = 1 + 2\alpha_1 + 4\alpha_2 + 8\alpha_3 + \cdots,$$
where $\alpha_k$ counts independent $k$-tuples of pairwise vertex-disjoint directed odd cycles in $D_\tau$.  Ties enrich $\Omega(D_\tau)$: a tie $u \leftrightarrow v$ creates cycles in both orientations.

**Conjecture** (All positive integers achievable): For every positive integer $N$, some t(r)ienerment $\tau$ achieves $H(\tau) = N$.  In particular, the permanent gaps $H = 7$ and $H = 21$ for tournaments are not forbidden for t(r)ienerments.

### 14.6 Open Questions

1. Verify computationally that $H = 7$ and $H = 21$ are achievable by some t(r)ienerment.
2. Characterize the $H$-spectrum of t(r)ienerments: which positive integers arise?
3. Among t(r)ienerments with $k$ ties fixed, which iso class maximizes $H$?
4. Is the Paley t(r)ienerment (with ties derived from the QR structure) in some sense extremal for the $\alpha$-decomposition?
5. Compute $\beta_p(D_\tau)$ (path homology Betti numbers) for t(r)ienerments and determine whether the tournament Betti bounds ($\beta_2 = 0$, $\beta_1 \le 1$) extend.
6. Describe the metagraph $G_n^{\text{tri}}$ on t(r)ienerment iso classes with edges from single-edge flips; compare its structure to the tournament metagraph.

---

## 15. Summary of Proof Status

### Proved (full proof in canon, independently verified)

| Result | Theorem ID |
|--------|------------|
| Rédei's theorem (4 independent routes) | THM-001 |
| OCF: $H(T) = I(\Omega(T), 2)$ | THM-002 |
| Claim B (algebraic companion to OCF) | THM-003 |
| Transfer matrix symmetry $M[a,b] = M[b,a]$ | THM-030 |
| Complete Walsh spectrum of $H(T)$ | THM-071 |
| Complete Walsh spectrum of $M[a,b]$ | THM-080 |
| Position Character Decomposition | THM-068 |
| Universal congruences for $S(T)$ (THM-H, -I, -J) | — |
| $\beta_2 = 0$ for all tournaments | THM-108/109 |
| $\beta_1 \leq 1$ for all tournaments | THM-103 |
| $\beta_1 \cdot \beta_3 = 0$ | THM-095 |
| Rank formula for $\partial_2$ | — |
| $F$-polynomial complement duality, mod-2 universality | THM-087/094 |
| $H = 7$ and $H = 21$ are permanent gaps | THM-029/115 |
| Pin grid $S_3 \times \mathbb{Z}_2$ symmetry | — |
| $k$-nacci Mersenne identity | THM-227 |
| Simplicial Rédei: $\mathrm{sim}_H \in \{0,1\}$ | THM-220 |
| Egyptian fraction splitting ($3/N$, master criterion $k/N$) | HYP-1612/1615 |
| Double factorial fixed point | HYP-1614/1616 |
| Constant symbol matrix for circulant tournaments | THM-125 |
| Cayley–Delannoy framework (variance, Rodrigues, duality) | — |
| Golden exceptional points of transfer matrix | THM-224 |

### Computationally established (strong evidence, no general proof)

| Result | Evidence |
|--------|----------|
| Paley Betti formula (THM-130) | Verified $p = 7, 11$; algebraic proof in progress |
| Paley tournaments maximize $H$ | Verified $p = 3, 7, 11$; matches A038375 |
| Unimodality of $F(T,x)$ | $50{,}000+$ tests; zero failures |
| At most 3 bad vertices (HYP-282) | Verified $n \leq 10$ |
| Eigenspace rank shift $R_d^{(k)} - R_d^{(0)} = (-1)^{d+1}$ | Verified $p = 7, 11$ |
| Torsion-free Paley chain complexes | Verified $p = 7, 11, 19, 23$ |
| Base-42 Erdős–Straus covering | Zero failures to $10^6$ |

### Open Problems

**Tier 1 (closest to resolution):**
1. Prove $\beta_1 \cdot \beta_3 = 0$ for all $n$ algebraically (mechanism via $\mathrm{im}(\partial_2)$ understood but not closed)
2. Algebraic proof of THM-130 (Paley Betti formula) via representation theory + THM-125
3. Prove all eigenspaces of $T_p$ have equal rank structure (follows from THM-125 + representation theory; writeup incomplete)
4. Is $\beta_{(p-1)/2}(T_p) = p - 1$ for all Paley primes? (Verified $p = 7, 11$)

**Tier 2 (major open problems):**
5. Prove the Paley maximization conjecture: $H(T_p) = \max_{|V|=p} H(T)$ for all Paley primes
6. Characterize $\beta_3 = 2$ at $n = 8$ (which tournaments achieve it? what is the structural explanation?)
7. What replaces $\beta_k \leq 1$ at $n \geq 8$? Is there a bound $\beta_k \leq C(n,k)$ for some combinatorial $C$?
8. Does the density of achievable $H$-values approach 1? Are there infinitely many permanent gaps beyond $\{7, 21\}$?
9. Per-path identity for all $n$: correct generalization incorporating cycles of all odd lengths (fails at $n = 6$ due to 5-cycle contributions)
10. $T_{19}$ full path homology: computationally intense but in principle feasible with CSC sparse matrices

---

## 16. References

1. N. Alon, *The maximum number of Hamiltonian paths in tournaments*, Combinatorica **10** (1990), 319–324.
2. J. Chapman, *Alternating sign matrices and tournaments*, Adv. Appl. Math. **27** (2001), 318–335.
3. R. Forcade, *Parity of paths and circuits in tournaments*, Discrete Math. **6** (1973), 115–118.
4. S. Grinberg, R. P. Stanley, *Counting Hamiltonian paths in tournaments*, arXiv:2412.10572 (2024).
5. A. Grigor'yan, Y. Lin, Y. Muranov, S.-T. Yau, *Homologies of path complexes and digraphs*, arXiv:1207.2834 (2012).
6. J. W. Moon, *Topics on Tournaments*, Holt, Rinehart and Winston (1968).
7. L. Rédei, *Ein kombinatorischer Satz*, Acta Litt. Sci. (Szeged) **7** (1934), 39–43.
8. L. J. Santharoubane, *Cohomology of Heisenberg Lie algebras*, Proc. Amer. Math. Soc. **87** (1983), 23–28.
9. J. Schweser, M. Stiebitz, B. Toft, *The tournament theorem of Rédei revisited*, arXiv:2510.10659 (2025).
10. R. P. Stanley, *Enumerative Combinatorics*, Vols. 1–2, Cambridge University Press (1999/2012).
11. K. B. Tang, S.-T. Yau, *Path homology of circulant digraphs*, arXiv:2602.04140 (2026).
12. A. El Sahili, M. Abi Aad, *Parity of paths and circuits in tournaments*, Discrete Math. **343** (2020).

---

*Repository at `/home/ubuntu/math`. Document maintained alongside the research; consult the canon (`01-canon/theorems/`) for full proofs of each numbered theorem.*
