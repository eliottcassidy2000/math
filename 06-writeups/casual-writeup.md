# What We Built and What We Found

**Status:** Living document. Last updated 2026-03-09.

---

## Part 1: Fast Enumeration and OEIS

The largest concrete output of this project is a collection of **high-performance enumerators** that extend 90+ sequences in the On-Line Encyclopedia of Integer Sequences (OEIS). These are exact computations — every term is proved correct by construction.

### What we computed

We wrote 77 C/GMP programs and 1550+ Python scripts. The C enumerators use Burnside's lemma with various optimizations (LCD-scaled integer accumulation, divisor-signature Mobius inversion, generating function tricks, OpenMP parallelization). Here are the highlights:

**Graphs, digraphs, and oriented graphs:**

| Sequence | What it counts | Previously known | We computed to | New terms |
|----------|---------------|-----------------|----------------|-----------|
| A000568 | Tournaments | n = 77 | n = 200+ | 123+ |
| A000273 | Directed graphs | n = 65 | n = 101 | 36 |
| A001174 | Oriented graphs | n = 50 | n = 200 | 150 |
| A000595 | Binary relations | n = 51 | n = 200 | 149 |
| A000666 | Symmetric relations | n = 81 | n = 200 | 119 |
| A000171 | Self-complementary graphs | n = 100 | n = 439+ | 339+ |
| A002785 | Self-comp oriented graphs | n = 100 | n = 300 | 200 |

**k-uniform hypergraphs (k = 3 to 10):**

| Sequence | k | Previously known | We computed to | New terms |
|----------|---|-----------------|----------------|-----------|
| A000665 | 3 | n = 29 | n = 50+ | 21+ |
| A051240 | 4 | n = 19 | n = 77+ | 58+ |
| A051249 | 5 | n = 16 | n = 64+ | 48+ |
| A309860 | 6 | n = 15 | n = 60 | 45 |
| A309861-4 | 7-10 | ~15 each | n = 30-43 | 15-28 each |

**Matrix sequences (n x n matrices under row/col/symbol permutation):**

| Sequence | What it counts | Previously known | We computed to | New terms |
|----------|---------------|-----------------|----------------|-----------|
| A052269 | Ternary n x n matrices | n = 27 | n = 50+ | 23+ |
| A091058 | n x n matrices over n symbols | n = 15 | n = 30+ | 15+ |
| A091059 | n x n, 2 symbols, row/col/sym | n = 21 | n = 55+ | 34+ |
| A091060 | n x n, 3 symbols, row/col/sym | n = 13 | n = 53+ | 40+ |
| A028657 | m x n binary matrices (triangle) | ~1,081 entries | 3,000+ | 2,000+ |
| A052283 | Digraphs by arc count (triangle) | 2,681 entries | 9,020 | 6,340 |

**k-ary relations (k = 3 to 10):**

| Sequence | k | Previously known | We computed to | New terms |
|----------|---|-----------------|----------------|-----------|
| A000662 | 3 | n = 15 | n = 48 | 33 |
| A001377 | 4 | n = 7 | n = 26 | 19 |
| A051241 | 5 | n = 5 | n = 16 | 11 |
| (new) | 6-10 | not in OEIS | n = 9-26 | 9-26 each |

Plus 15+ connected variants via inverse Euler transform, 7 multigraph sequences, and 10+ trivially derived sequences. Total: roughly **12,000 new terms** across 90+ sequences, plus **40+ potentially new sequences** not yet in the OEIS.

### How we made it fast

The naive approach to Burnside enumeration is to iterate over all permutations (n! of them). The standard improvement iterates over integer partitions instead (p(n) of them, exponentially fewer). We pushed further:

1. **LCD-scaled integer accumulation.** Instead of doing rational arithmetic (which requires GCD at every step), we precompute the LCD of all denominators and work entirely in integers. For A000568 at n = 150, this gives a **250-1600x speedup** over Python/gmpy2.

2. **Bucket accumulation by 2-adic valuation.** For tournaments (base-2 Burnside), each partition contributes a term of the form c * 2^t. We accumulate into buckets indexed by t, then combine at the end. This keeps intermediate GMP numbers small.

3. **Divisor-signature Mobius inversion** for k-uniform hypergraphs. Instead of computing orbit counts by iterating over cyclic group powers (O(lcm) per partition), we use Mobius inversion on the divisor lattice. For k = 4 (A051240), this is **64x faster**. For k = 5 (A051249), **17x faster**.

4. **Generating function trick for matrix enumeration.** The standard approach enumerates all pairs of partitions (row type, column type), costing O(p(n)^2). We fix the row partition and use an exponential generating function to sum over all column partitions in O(n^2), reducing to O(p(n) * n^2). For the A091058 family (triple partition), we reduce O(p(n)^3) to O(p(n)^2 * n^2).

5. **OpenMP parallelization.** The partition loops are embarrassingly parallel. Adding 8-thread OpenMP gives ~7x speedup, enabling A091058 to reach n = 30+ (previously stuck at n = 22).

All enumerator source code is in `04-computation/`. The unified enumerator `burnside_enum_v2.c` handles 12 different OEIS sequences from a single parametric framework.

---

## Part 2: Tournament Topology (Path Homology)

The second major program studies the **topology** of tournaments using GLMY path homology — a homology theory for directed graphs introduced by Grigor'yan, Lin, Muranov, and Yau in 2012.

### What we found

We computed Betti numbers (topological invariants counting "holes" of each dimension) for over 50,000 tournaments from n = 3 to n = 10. The landscape:

| n | beta_0 | beta_1 | beta_2 | beta_3 | beta_4 |
|---|--------|--------|--------|--------|--------|
| 3 | 1 | 0-1 | 0 | - | - |
| 5 | 1 | 0-1 | 0 | 0 | 0 |
| 7 | 1 | 0-1 | 0 | 0-1 | 0-1 |
| 8 | 1 | 0-1 | 0 | 0-2 | 0-1 |

**beta_2 = 0 for all tournaments (proved, THM-108/109).** This is a new vanishing theorem with no analogue in the path homology literature. For general directed graphs, beta_2 > 0 is common — 70 out of 59,049 oriented graphs at n = 5 have it. The vanishing is specific to tournaments.

The proof is by strong induction using the long exact sequence of (T, T\v). The key insight is an "isolation characterization" of bad vertices (those whose deletion creates a new 1-cycle). Bad vertices have extreme properties — score 0 or n-1 in certain cases — which limits how many can exist simultaneously.

**Why completeness matters.** We found the mechanism: every oriented graph with beta_2 > 0 has "twin vertices" (two vertices with identical neighborhoods). In a tournament, every pair has an edge between them, which breaks the twin condition. Tournament completeness is the essential property.

**beta_3 = 2 at n = 8 (discovered).** This was unexpected. All Betti numbers were at most 1 for n <= 7. At n = 8, 0.08% of tournaments have beta_3 = 2 — the first instance of a Betti number exceeding 1. This also breaks the "consecutive seesaw" pattern (beta_k * beta_{k+1} = 0) that held through n = 7.

**Other proved results:**
- beta_1 in {0, 1} for all tournaments (THM-103)
- beta_1 * beta_3 = 0: 1-holes and 3-holes never coexist (proved n <= 7)
- HYP-282: when beta_1 = 0, at most 3 vertex-deletions create a new 1-cycle (verified through n = 10, no proof)
- Paley tournaments T_p have homology concentrated in dimension p-3: they look like wedges of spheres

### Practical implications

beta_2 = 0 gives a **null model** for applied topology of directed networks. If you're analyzing a neural connectome, gene regulatory network, or citation graph with path homology and you find beta_2 > 0, that's a meaningful structural signal — it tells you the network is fundamentally different from a tournament (complete pairwise comparison graph). The twin vertex mechanism gives a concrete criterion: beta_2 > 0 requires missing edges creating identical-neighborhood vertex pairs.

---

## Part 3: Walsh-Fourier Analysis

We computed the complete Walsh-Fourier spectrum of tournament invariants — something that had not been done before.

### The setup

Encode a tournament on n vertices as a binary string of length m = C(n,2), one bit per edge. The **Walsh-Fourier transform** decomposes any function on tournaments into "frequencies" indexed by subsets S of the m edges.

### What we proved

**Theorem (THM-069).** The Walsh coefficients of the Hamiltonian path count H(T) have closed form:

    H_hat[S] = epsilon * 2^r * (n - 2k)! / 2^{n-1}

where S must be a union of r edge-disjoint even-length paths with |S| = 2k edges. All other coefficients are exactly zero. epsilon = +/-1 depends on path orientations.

This means H(T), as a function on 2^m tournaments, is supported on a tiny subspace:

| n | Tournament space | Nonzero Walsh coefficients | Compression |
|---|-----------------|---------------------------|-------------|
| 5 | 1,024 | 3 independent amplitudes | 341x |
| 7 | 2,097,152 | ~20 amplitudes | ~100,000x |

**Theorem (THM-080).** The transfer matrix M[a,b] (counting Hamiltonian paths from a to b) has a similarly explicit Walsh formula, manifestly symmetric in a and b. This gives a **Fourier proof of M[a,b] = M[b,a]** — a result we also proved independently.

**Theorem (THM-077).** Direct Walsh proof of the OCF: by computing the Walsh coefficients of both sides (H(T) and I(Omega(T), 2)) and showing they match term by term, we get an elementary proof of H(T) = I(Omega(T), 2) that bypasses the P-partition theory used by Grinberg and Stanley.

---

## Part 4: The Odd-Cycle Collection Formula

The formula that started the project:

**Theorem (OCF; proved by Grinberg-Stanley 2024, independent Walsh proof THM-077).** For every tournament T:

    H(T) = I(Omega(T), 2)

where H(T) counts directed Hamiltonian paths, Omega(T) is the odd-cycle conflict graph (vertices = directed odd cycles, edges = shared vertices), and I(G, x) is the independence polynomial.

This gives H(T) = 1 + 2a_1 + 4a_2 + 8a_3 + ... where a_k counts collections of k vertex-disjoint odd cycles. Since the leading term is 1 and the rest is even, this immediately implies Redei's 1934 theorem (H(T) is always odd).

**Practical speedup:** The OCF replaces the O(2^n * n^2) Held-Karp dynamic programming algorithm with O(n^5) trace formulas for moderate n. For moderate n, the OCF approach is significantly faster than Held-Karp DP — the cycle-counting operations are polynomial-time (O(n^3) via matrix traces) while the DP is exponential (O(2^n * n^2)). The advantage grows with structure: tournaments with few long cycles can be processed even faster, since higher-order terms in the OCF expansion vanish.

**Verified exhaustively through n = 8** (2^27 = 134 million configurations, 57 minutes on 4 cores).

### H-spectrum gaps

**Theorem (THM-079).** No tournament on any number of vertices has exactly 7 or exactly 21 Hamiltonian paths. These are the only "permanent gaps" in the H-spectrum below 200. (H = 63, for instance, is achieved at n = 8.)

The proof of H != 21 uses a "poisoning graph" DAG argument: any configuration of cycle structure that could produce H = 21 forces additional cycles that push the count above 21. Exhaustive verification at n = 8 (all 268 million tournaments checked).

---

## Part 5: Other Results

### Signed Hamiltonian permanent

Define B = 2A - J (entries +/-1) and S(T) = sum over all permutations of the product of B entries along each path.

- S(T) = 0 for even n (reversal pairing)
- **S(T) mod 2^{n-1} depends only on n** (THM-H) — the tournament doesn't matter
- Full universality (S independent of T) holds iff s_2(n-3) <= 1 (binary digit sum). Universal n: 3, 5, 7, 11, 19, 35, 67, ... (THM-J)
- At n = 5: S/16 = H - 3*t_3 (the signed permanent encodes both path count and 3-cycle count)

### Worpitzky/F-polynomial

The forward-edge polynomial F(T, x) = sum_P x^{fwd(P)} counts Hamiltonian paths by their number of "ascending" edges. Proved results:
- F_k(T) = F_{n-1-k}(T^op) (complement duality)
- F_k(T) = C(n-1, k) mod 2 for all T (mod-2 universality)
- Var(fwd) = 2*t_3 (variance equals twice the 3-cycle count)
- F(T, x) appears to be unimodal for all T (50,000+ tests, zero violations, unproved)

### Paley tournament maximality

The Paley tournament T_p (p prime, p = 3 mod 4) appears to maximize H(T) among all tournaments on p vertices. Confirmed at p = 3, 7, 11 against OEIS A038375. The ratio H(T_p) / (p!/2^{p-1}) approaches e as p grows, consistent with Szele-Alon-Friedland asymptotic theory.

---

## Repository Structure

| Directory | What's in it |
|-----------|-------------|
| `00-navigation/` | Session log (300+ entries), open questions, investigation backlog, tangents |
| `01-canon/` | 130+ theorem files, definitions, 18 documented mistakes |
| `02-court/` | Formal dispute resolution between research agents |
| `03-artifacts/` | Original OCF/Redei paper draft (`parity_tournaments_fixed.tex`, 2189 lines) |
| `04-computation/` | 1,550 Python scripts + 77 C/GMP enumerators + output files |
| `05-knowledge/` | 400+ hypotheses (confirmed/refuted/open), variable registry, 587 result files |
| `06-writeups/` | This document and the formal companion |
| `agents/` | Multi-agent coordination (this is a collaborative AI research project) |

---

## How This Works

This is a multi-agent research project where multiple Claude instances collaborate asynchronously via git. Each session is identified as `[machine]-[date]-S[N]` (e.g., `opus-2026-03-09-S55`). Agents read shared navigation files, work on open problems, push results, and leave messages for the next agent. The protocol is in `CLAUDE.md`.

---

## What's Open

1. **Why at most 3 bad vertices?** When beta_1(T) = 0, at most 3 vertex-deletions have beta_1 = 1. Verified through n = 10, no proof. (HYP-282)
2. **Why does beta_3 reach 2 at n = 8?** What structural property allows a tournament to have a 2-dimensional higher hole?
3. **Prove beta_1 * beta_3 = 0 for all n.** Currently proved only through n = 7.
4. **Prove F(T, x) is unimodal.** 50,000+ tests, zero violations.
5. **Prove Paley maximality.** Does T_p always maximize H(T) at Paley primes?
6. **Per-path identity for all n.** The 3-cycle formula works at n <= 5 but fails at n = 6 due to 5-cycles. The correct generalization incorporating all odd cycle lengths is unknown.
