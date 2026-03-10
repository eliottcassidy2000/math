# Tournament Mathematics: OEIS Extensions, Path Homology, and Walsh-Fourier Analysis

This repository contains an ongoing computational and theoretical research project on **tournaments** (complete directed graphs). It has produced new terms for **90+ OEIS sequences**, proved that tournaments have no 2-dimensional topological holes (beta_2 = 0), and computed the complete Walsh-Fourier spectrum of Hamiltonian path counts.

---

## OEIS Contributions

### Summary

| Category | Sequences extended | Approximate new terms |
|----------|-------------------|----------------------|
| Graphs, digraphs, oriented graphs | 15 | 1,500+ |
| Self-complementary variants | 8 | 700+ |
| k-uniform hypergraphs (k=3..10) | 9 | 250+ |
| k-ary relations (k=3..10) | 8 | 130+ |
| k-ary n x n matrices (k=2..10) | 20+ | 800+ |
| Triangle sequences (by edge/arc count) | 2 | 8,500+ |
| Connected variants (inverse Euler) | 15+ | 600+ |
| Multigraph sequences | 7 | 230+ |
| Trivially derived (sums, halves) | 10+ | 600+ |
| **Total** | **90+** | **~12,000+** |

Additionally, **40+ potentially new sequences** have been computed (connected k-uniform hypergraphs, k-ary relations for k >= 6, connected multigraph variants) that do not yet appear in the OEIS.

### Headline Extensions

| Sequence | Description | OEIS had | We computed to | New terms |
|----------|-------------|----------|----------------|-----------|
| **A000568** | Tournaments on n nodes | n=77 | n=200+ | **123+** |
| **A000273** | Directed graphs | n=65 | n=101 | **36** |
| **A001174** | Oriented graphs | n=50 | n=200 | **150** |
| **A000595** | Binary relations | n=51 | n=200 | **149** |
| **A000171** | Self-complementary graphs | n=100 | n=439+ | **339+** |
| **A002785** | Self-comp oriented graphs | n=100 | n=300 | **200** |
| **A051240** | 4-uniform hypergraphs | n=19 | n=77+ | **58+** |
| **A052283** | Digraphs by arc count (triangle) | 2,681 entries | 9,020 entries | **6,340** |
| **A028657** | m x n binary matrices (triangle) | ~1,081 entries | 3,000+ entries | **2,000+** |
| **A091060** | n x n matrices, 3 symbols, row/col/sym | n=13 | n=53+ | **40+** |

### Algorithmic Innovations

The computation of these extensions required several new algorithmic techniques:

- **LCD-scaled integer accumulation** for A000568: 250-1600x speedup over gmpy2 DP, enabling n=200+ (3000+ digit values). Uses integer-only arithmetic with LCD scaling to avoid all rational number overhead.
- **Divisor-signature Mobius optimization** for k-uniform hypergraphs: 64-130x speedup for A051240/A051249 by replacing Burnside orbit-counting with closed-form cycle index evaluation.
- **Generating function approach** for matrix sequences: Reduces O(p(n)^2) pair-partition enumeration to O(p(n) * n^2) via single-partition exponential GF recurrence.
- **Triple partition GF trick** for A091058 family: Reduces O(p(n)^3) to O(p(n)^2 * n^2).
- **OpenMP parallelization**: ~7-8x speedup on 8 cores for the most expensive sequences.
- **Unified C+GMP framework**: Single enumerator handles 12 graph/digraph sequences via parametric edge-counting formulas.

All enumerator source code is in `04-computation/`. See `05-knowledge/results/burnside_enum_extensions.md` for the complete catalog.

### Enumerator Files

| File | Sequences | Technique |
|------|-----------|-----------|
| `burnside_enum_v2.c` | A000568, A000273, A000595, A000088, A000666, A001174, A002785, A000171, A003086, A005639, A002499, A002854 | Compressed partition + GMP |
| `a051240_gmp_enum.c` | A051240 | Divisor-signature Mobius (64x faster) |
| `k_uniform_fast_enum.c` | Any k-uniform hypergraph | General Burnside orbit counting |
| `a008406_gmp.c` | A008406 | Pair orbit generating function |
| `a052283_gmp.c` | A052283 | Directed pair orbit GF |
| `kary_matrix_fast_omp.c` | A052269-A052272, A246112-A246116 | GF + OpenMP |
| `a242095_omp.c` | A091058-A091062, A246122-A246126, A242095 triangle | Triple GF + OpenMP |
| `binary_matrix_mn_fast_gmp.c` | A028657 | m x n pair orbit GF |
| `k_ary_relations_gmp.c` | A000662, A001377, A051241, k >= 6 | General k-ary Burnside |
| `euler_transform.py` | 15+ connected variants | Inverse Euler transform |

---

## Mathematical Results

### The Odd-Cycle Collection Formula (OCF)

**Theorem (Grinberg-Stanley, 2024; THM-077, independent Walsh proof).** For every tournament T:

    H(T) = I(Omega(T), 2)

where H(T) counts directed Hamiltonian paths, Omega(T) is the odd-cycle conflict graph, and I(G, x) is the independence polynomial. This gives H(T) = 1 + 2a_1 + 4a_2 + ... where a_k counts collections of k vertex-disjoint odd cycles.

**Practical impact:** Replaces O(2^n * n^2) Held-Karp DP with O(n^5) trace formulas for moderate n.

**What's new here:** We give an independent elementary proof (THM-077) via Walsh-Fourier analysis, bypassing the P-partition theory used by Grinberg-Stanley. We also compute the complete Walsh spectrum (below).

### Walsh-Fourier Spectrum of Tournament Invariants

**Theorem (THM-069).** The Walsh-Fourier transform of H(T) has closed form:

    H_hat[S] = epsilon * 2^r * (n - 2k)! / 2^{n-1}

where S is a union of r edge-disjoint even-length paths with |S| = 2k edges, and epsilon = +/-1.

**Theorem (THM-080).** The transfer matrix M[a,b] (counting Hamiltonian paths from a to b) has a similarly explicit Walsh spectrum, manifestly symmetric in a and b. This gives a Fourier proof of M[a,b] = M[b,a].

**What's new:** No prior work computes the full Walsh spectrum of tournament path counts. The extreme sparsity (3 independent amplitudes at n=5 on a 1024-dimensional space) means tournament invariants can be learned from very few samples.

### GLMY Path Homology of Tournaments

**Theorem (THM-108/109).** For every tournament T, beta_2(T) = 0 in GLMY path homology.

This is a new vanishing result with no prior analogue. For general directed graphs, beta_2 > 0 is common. The proof uses induction via the long exact sequence of (T, T\v), with an isolation characterization of "bad" vertices.

**Additional proved results:**
- beta_1(T) in {0, 1} for all tournaments (THM-103)
- beta_1 * beta_3 = 0: mutual exclusivity of 1-holes and 3-holes (proved n <= 7)
- Rank formula: rank(d_2) = C(n,2) - n + 1 - beta_1(T)

**New discoveries:** beta_3 first appears at n=6 (1% of tournaments), reaches 2 at n=8. beta_4 onset at n=7, reaching 6 for the Paley tournament T_7. beta_5 first appears at n=8. The "defect wave" pattern: beta_1 prevalence drops (29.7% -> 14.6% -> 5.8% -> 1%) while beta_3 rises (0% -> 1% -> 7.2% -> 21%) as n grows.

| n | beta_0 | beta_1 | beta_2 | beta_3 | beta_4 | beta_5 |
|---|--------|--------|--------|--------|--------|--------|
| 3 | 1 | 0-1 | 0 | - | - | - |
| 4 | 1 | 0-1 | 0 | 0 | - | - |
| 5 | 1 | 0-1 | 0 | 0 | 0 | - |
| 6 | 1 | 0-1 | 0 | 0-1 | 0 | 0 |
| 7 | 1 | 0-1 | 0 | 0-1 | 0-6 | 0 |
| 8 | 1 | 0-1 | 0 | 0-2 | 0-5 | 0-1 |

### Signed Hamiltonian Permanent

**Theorem (THM-H).** S(T) mod 2^{n-1} depends only on n, not on the tournament.

**Theorem (THM-J).** Full universality (S(T) independent of T) holds iff s_2(n-3) <= 1, where s_2 is the binary digit sum. Universal n: 3, 5, 7, 11, 19, 35, 67, ...

### H-Spectrum Gaps

**Theorem (THM-079).** H(T) = 7 and H(T) = 21 are impossible for any tournament on any number of vertices. These are the only permanent gaps in [1, 200] through n = 9.

### Paley Tournament Maximality

**Conjecture (verified through n = 11).** Among all tournaments on p vertices (p prime, p = 3 mod 4), the Paley tournament maximizes H(T). Matches OEIS A038375 at p = 3, 7, 11.

---

## Repository Structure

| Directory | Purpose |
|-----------|---------|
| `00-navigation/` | Session log, open questions, investigation backlog, tangents |
| `01-canon/` | Proved theorems, definitions, documented mistakes |
| `02-court/` | Formal dispute resolution between research agents |
| `03-artifacts/` | LaTeX paper draft, code artifacts |
| `04-computation/` | All Python and C scripts (enumerators, homology, Walsh) |
| `05-knowledge/` | Hypotheses index, variable registry, computational results |
| `06-writeups/` | Formal and casual research summaries |
| `agents/` | Multi-agent coordination (inbox/outbox, session state) |

### Key Files

- `03-artifacts/drafts/tournaments_comprehensive.tex` — comprehensive paper (enumeration + Walsh + homology)
- `03-artifacts/drafts/parity_tournaments_fixed.tex` — original OCF/Rédei paper draft (2,189 lines)
- `04-computation/burnside_enum_v2.c` — unified graph/digraph enumerator (12 OEIS sequences)
- `04-computation/path_homology_v2.py` — GLMY path homology computation
- `05-knowledge/results/burnside_enum_extensions.md` — complete OEIS extension catalog
- `05-knowledge/hypotheses/INDEX.md` — 400+ hypotheses (confirmed/refuted/open)
- `06-writeups/formal-writeup.md` — detailed research summary with all theorem statements

---

## How This Repository Works

This is a multi-agent research project where multiple Claude instances collaborate asynchronously via git. Each session is identified as `[machine]-[date]-S[N]`. Agents read shared navigation files, work on open problems, and push results. The `CLAUDE.md` file contains the full protocol.

---

## References

- D. Grinberg, R.P. Stanley. *Counting Hamiltonian paths in tournaments.* arXiv:2412.10572 (2024).
- A. Grigor'yan, Y. Lin, Y. Muranov, S.-T. Yau. *Homologies of path complexes and digraphs.* arXiv:1207.2834 (2012).
- L. Redei. *Ein kombinatorischer Satz.* Acta Litt. Sci. Szeged 7 (1934), 39-43.
- J.W. Moon. *Topics on Tournaments.* Holt, Rinehart and Winston (1968).
- R. Forcade. *Parity of paths and circuits in tournaments.* Discrete Math. 6 (1973), 115-118.
- R. Tang, S.-T. Yau. *Homology of tournaments and path homology.* arXiv:2602.04140 (2026).
