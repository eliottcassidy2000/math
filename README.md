# Tournament Mathematics: Theory, Algorithms, and Tools

This repository is an ongoing research project on **tournaments** (complete directed graphs). It spans pure mathematics, algorithm engineering, and practical tools — all connected by a single algebraic structure: the formal group **F(x,y) = (x+y)/(1+xy)**.

**Highlights:**
- **90+ OEIS sequences extended** with 12,000+ new terms
- **114 proved theorems** including beta_2 = 0 for all tournaments (GLMY path homology)
- **Production tools** for pairwise ranking at 300K+ observations/sec
- **The formal group logarithm** (arctanh) as a universal linearizer connecting number theory, spectral theory, and practical computation

---

## Practical Tools

### FormalRank — One-pass pairwise ranking
**`04-computation/formalrank.py`**

Ranks items from pairwise comparisons using the formal group logarithm. No iteration, no hyperparameters, no convergence criterion. Evidence is perfectly additive: `score += arctanh(coupling)`.

```python
from formalrank import FormalRank
ranker = FormalRank()
ranker.add("Claude", "GPT-4", wins=156, losses=134)
ranker.add("Claude", "Gemini", wins=178, losses=112)
print(ranker.summary())
```

**Performance:** 560K observations/sec. **Applications:** LLM evaluation, A/B testing, sports, elections.

### BoostRanker — Three-channel ranking via the boost trichotomy
**`04-computation/boost_ranker_s91k.py`**

Every comparison generates a Cayley boost Q = (1+x)/(1-x) that decomposes into three independent signals:

| Channel | What it measures | Use |
|---------|-----------------|-----|
| **INERT** (mod 2) | Who wins? | Binary ranking |
| **RAMIFIED** (mod 3) | By how much? | Confidence estimation |
| **SPLIT** (mod 7) | Was this expected? | Upset/anomaly detection |

Includes automatic cycle detection (finds interviewer disagreements) and upset flagging.

### StreamingComparator — Real-time ranking engine
**`04-computation/practical_lattice_s91j.py`**

Streaming variant optimized for dashboards and live scoring. O(1) amortized updates, auto-stopping rules, win probability estimation, and "how many more comparisons needed?" queries.

**Performance:** 307K observations/sec at 1M scale.

### SpectralRanker — Exact algebraic ranking
**`04-computation/spectral_ranker.py`**

Uses the spectral gap theorem (gap = 4/C(n,2), exact) to compute ranking confidence with no floating-point error. Identifies high-leverage comparisons — which single result, if reversed, changes the ranking most.

### Modular Rank Library
**`04-computation/mod_rank_library.py`**

Small-prime Gaussian elimination for integer matrix rank. 8x memory reduction (uint8 storage), two-prime certification for mathematical correctness. Powers the path homology computation for tournaments up to n=19.

### Tournament Toolkit
**`04-computation/tournament_toolkit.py`**

Four modules: TilingGrid (triangular encoding), FlipGraph (arc transitions), RankAmbiguity (H-count based ranking), TournamentTDA (topological features for ML).

---

## Mathematical Theory

### The Formal Group F(x,y) = (x+y)/(1+xy)

This is the thread connecting everything. It is simultaneously:
- The **relativistic velocity addition** law
- The **Poincare disk** group law (hyperbolic geometry)
- The **tanh addition** formula
- The **Bayesian evidence aggregation** law (log-odds addition)
- A **formal group of height 1** at every odd prime and **height infinity** at p=2

Its logarithm **arctanh(x) = x + x^3/3 + x^5/5 + ...** linearizes all of the above. The coefficients 1/(2k+1) are reciprocals of odd numbers — the same odd numbers that appear as Hamiltonian path counts of tournaments.

### Key Results

**The Odd-Cycle Collection Formula (OCF).** For every tournament T: H(T) = I(Omega(T), 2), where I is the independence polynomial and Omega is the odd-cycle conflict graph. This gives H(T) = 1 + 2a_1 + 4a_2 + ... (always odd, by Redei's theorem). The formal group is supersingular at p=2, which is WHY H(T) is always odd.

**Walsh-Fourier Spectrum (THM-071).** The Walsh transform of H(T) has closed form with only 3 independent amplitudes at n=5 on a 1024-dimensional space. Extreme sparsity.

**Path Homology (THM-108/109).** beta_2(T) = 0 for every tournament. beta_1 in {0,1}. beta_1 * beta_3 = 0 (mutual exclusivity of holes).

**Forbidden H-values (THM-029/079).** H = 7 and H = 21 are impossible for any tournament. These are C(7,1) and C(7,5) — binomial coefficients appearing in the 7-torsion polynomial of the formal group.

**Spectral Gap (verified n=3..6).** The flip chain on isomorphism classes has gap = 4/C(n,2) and eigenvalue denominator = odd part of C(n,2).

**The Rapidity Lattice.** The Cayley transform Q(lambda) = (1+lambda)/(1-lambda) of eigenvalues factors over {2, 3, 7} (the Hurwitz primes). The rapidities arctanh(eigenvalues) form a rank-3 lattice Z*ln(2) + Z*ln(3) + Z*ln(7), guaranteed by Baker's theorem.

**Heat Kernel Algebraicity.** K(ln(q)) = Tr(q^{-T}) is a Laurent polynomial in q^{1/D}, not a transcendental function. Discreteness is an illusion of the transcendental representation.

**The Adelic Tournament Space.** Tournament eigenvalues live on a truncated adelic space A_T(n) = R x prod_{p|D} Z/p^e Z, approaching the full adeles as n grows. The flip chain is a Hecke operator on this space.

### The Boost Trichotomy

42 = 2 * 3 * 7 = INERT * RAMIFIED * SPLIT. Three independent axes:

| Axis | Prime | Eisenstein type | Computation type | Tournament role |
|------|-------|----------------|-----------------|-----------------|
| Parity | 2 | INERT | Lossless | Symmetric/antisymmetric |
| Curvature | 3 | RAMIFIED | Lossy | Eigenvalue structure |
| Position | 7 | SPLIT | Crystallization | Forbidden values |

Both 3 and 7 have class number 1 in their imaginary quadratic fields, ensuring unique factorization in the rapidity lattice.

---

## OEIS Contributions

| Sequence | Description | OEIS had | We computed to | New terms |
|----------|-------------|----------|----------------|-----------|
| **A000568** | Tournaments on n nodes | n=77 | n=200+ | **123+** |
| **A000273** | Directed graphs | n=65 | n=101 | **36** |
| **A000171** | Self-complementary graphs | n=100 | n=439+ | **339+** |
| **A052283** | Digraphs by arc count (triangle) | 2,681 entries | 9,020 entries | **6,340** |
| **A028657** | m x n binary matrices (triangle) | ~1,081 entries | 3,000+ entries | **2,000+** |

**Total: 90+ sequences, 12,000+ new terms, 40+ potentially new sequences.**

### Algorithmic Innovations

- **LCD-scaled integer accumulation** for A000568: 250-1600x speedup, enabling n=200+
- **Divisor-signature Mobius** for k-uniform hypergraphs: 64-130x speedup
- **Generating function approach** for matrix sequences: O(p(n)^2 * n^2) from O(p(n)^3)
- **Denominator killing primes**: when p = C(n,2)+1 is prime, eigenvalue fractions become trivial mod p (works at n=4,5,8,9,12,13,17,...)
- **Formal group acceleration**: [n](x) = tanh(n*arctanh(x)) gives 3589x speedup over iterative group law

---

## Repository Structure

| Directory | Purpose |
|-----------|---------|
| `00-navigation/` | Session log, open questions, investigation backlog |
| `01-canon/` | 114 proved theorems, definitions, documented mistakes |
| `02-court/` | Formal dispute resolution between research agents |
| `03-artifacts/` | Papers, engineering specs, OEIS submissions |
| `04-computation/` | All scripts: enumerators, homology, ranking tools, spectral analysis |
| `05-knowledge/` | 400+ hypotheses, variable registry, computational results |
| `06-writeups/` | Research summaries |
| `07-reflections/` | Cross-domain patterns (adelic geometry, boost trichotomy, three types) |
| `agents/` | Multi-agent coordination system |

### Key Files

**Tools:**
- `04-computation/formalrank.py` — one-pass pairwise ranking (production)
- `04-computation/boost_ranker_s91k.py` — three-channel ranking via trichotomy
- `04-computation/practical_lattice_s91j.py` — streaming comparator (307K obs/sec)
- `04-computation/spectral_ranker.py` — exact algebraic ranking
- `04-computation/mod_rank_library.py` — modular rank with 8x memory reduction
- `04-computation/tournament_toolkit.py` — tournament analysis toolkit (4 modules)

**Enumerators:**
- `04-computation/burnside_enum_v2.c` — unified graph/digraph enumerator (12 OEIS sequences)
- `04-computation/k_uniform_fast_enum.c` — general k-uniform hypergraph counter
- `04-computation/fast_a000568_v2_s90di.py` — fast tournament counting via partition sum

**Theory:**
- `03-artifacts/drafts/tournaments_comprehensive.tex` — comprehensive paper
- `03-artifacts/drafts/engineering-synthesis-2026-03-10-S53.md` — engineering roadmap (12 domains)
- `07-reflections/adelic-tournament-geometry.md` — the adelic structure of eigenvalues
- `07-reflections/three-types.md` — lossless / lossy / crystallization trichotomy

---

## How This Repository Works

This is a multi-agent research project where multiple Claude instances collaborate asynchronously via git. Each session is identified as `[machine]-[date]-S[N]`. The `CLAUDE.md` file contains the full protocol including startup sequence, session logging, and end-of-session close-out.

---

## References

- D. Grinberg, R.P. Stanley. *Counting Hamiltonian paths in tournaments.* arXiv:2412.10572 (2024).
- A. Grigor'yan, Y. Lin, Y. Muranov, S.-T. Yau. *Homologies of path complexes and digraphs.* arXiv:1207.2834 (2012).
- L. Redei. *Ein kombinatorischer Satz.* Acta Litt. Sci. Szeged 7 (1934), 39-43.
- J.W. Moon. *Topics on Tournaments.* Holt, Rinehart and Winston (1968).
- R. Forcade. *Parity of paths and circuits in tournaments.* Discrete Math. 6 (1973), 115-118.
- R. Tang, S.-T. Yau. *Homology of tournaments and path homology.* arXiv:2602.04140 (2026).
- A. Baker. *Linear forms in the logarithms of algebraic numbers.* Mathematika 13 (1966), 204-216.
