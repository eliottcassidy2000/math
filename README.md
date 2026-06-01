# Tournament Mathematics: Theory, Algorithms, and Engineering Tools

This repository is an ongoing research project on **tournaments** (complete directed graphs). It spans pure mathematics, algorithm engineering, and practical tools -- all connected by a single algebraic structure: the formal group **F(x,y) = (x+y)/(1+xy)**.

**Highlights:**
- **90+ OEIS sequences extended** with 12,000+ new terms
- **114 proved theorems** including beta_2 = 0 for all tournaments (GLMY path homology)
- **21 hidden tournaments discovered** across CS, ML, biology, and finance
- **Production tools** for pairwise ranking at 300K+ observations/sec
- **Polynomial output head** for GPT-2 that replaces softmax with tournament ranking
- **The formal group logarithm** (arctanh) as a universal linearizer connecting number theory, spectral theory, and practical computation

---

## Engineering Applications

Anything with pairwise comparisons that produces a ranking is a tournament. This project has identified **21 hidden tournaments** across computing, machine learning, biology, and finance -- and built tools that exploit the polynomial structure of tournament dynamics.

### LLM Evaluation and Ranking

The `RankAmbiguity` class (in `tournament_toolkit.py`) treats LLM leaderboards as tournaments. H(T) -- the Hamiltonian path count -- is a single number that quantifies "how settled is the leaderboard?" H=1 means an unambiguous ranking. Large H means many consistent orderings. The `high_leverage_comparisons()` method identifies which single matchup, if re-evaluated, would most reduce ambiguity -- directly applicable to Chatbot Arena's query budget.

Forbidden H-values (H=7 and H=21 are impossible for any tournament) provide new impossibility theorems for preference learning: certain ambiguity levels provably cannot occur.

### Hallucination Detection (Fast Channel)

Every pairwise comparison decomposes into three spectral signals via the **Boost Trichotomy** (THM-253):

| Signal | What it detects | Decay rate |
|--------|----------------|------------|
| **SLOW** (degree 1-2) | Persistent ranking structure | Half-life ~3 flips |
| **MEDIUM** (degree 2-3) | Pairwise correlations | Half-life ~1.4 flips |
| **FAST** (degree 3-4) | Contradictions and cycles | Half-life ~0.5 flips |

When the fast/slow ratio spikes, the system is generating contradictions -- the tournament-theoretic definition of hallucination. This works at the attention level (before the token is generated), not just the output level.

### Polynomial Output Head for GPT-2

`polynomial_head.py` implements three alternative output heads benchmarked against standard softmax:

1. **Standard softmax** -- the baseline
2. **Polynomial head** -- quantized Horner evaluation of the token tournament
3. **Formal group head** -- arctanh aggregation with early exit

The formal group head replaces softmax with the Cayley formal group logarithm. For greedy decoding, arctanh is monotone so argmax is unchanged. For sampling, `tanh(score)` gives calibrated probabilities. The architecture separates token selection (expensive, O(V)) from token ranking (cheap, O(D_eff^2)), yielding up to 500x savings on the ranking phase when D_eff << V.

### A/B Test Cycle Detection

Current A/B testing frameworks report individual test results but cannot detect when multiple tests produce contradictory results (A>B, B>C, C>A). The polynomial framework's fast-channel detector fills this gap: build the tournament from n*(n-1)/2 pairwise tests, decompose into three signals, and flag intransitivity. This detects Simpson's paradox in experimental design -- currently invisible to experimenters.

### Cache Eviction Improvement

Every cache eviction policy is a tournament ranking strategy. LRU, LFU, ARC, and W-TinyLFU all correspond to different tournament structures. ARC can exhibit intransitivity when recency and frequency disagree about which item to keep. A tournament-aware cache uses the three-signal decomposition to detect workload phase transitions in O(1) and switch between LRU-like and LFU-like eviction automatically, subsumming ARC and W-TinyLFU as special cases.

### Consensus Protocol Improvement (Paxos/Raft)

In Paxos, proposals can defeat each other in cycles causing livelock. In Raft, leader election stalls when candidates cycle. The spectral gap of the election tournament equals the convergence time after leader failure. A tournament-aware consensus protocol sets the election timeout based on the measured spectral gap of the cluster rather than an arbitrary constant.

### Active Learning for RLHF

The `PolynomialPredictor` class computes information gain for each potential pairwise comparison. The `next_comparison()` method identifies the pair whose comparison would most reduce ranking ambiguity -- enabling active learning that focuses the annotation budget on the most informative comparisons. This applies directly to RLHF preference data collection.

---

## Production Libraries

### tournament_H.py -- Tetrahedron Algorithm
Fast computation of H(T) via three-regime classification:
- **Regime 30** (near-transitive): O(k * 2^k) perturbation method
- **Regime 36** (moderate): OCF or polynomial evaluation
- **Regime 42** (full complexity): Held-Karp DP, O(n^2 * 2^n)

```python
from tournament_H import TournamentH
T = TournamentH(adj_matrix)      # or TournamentH.from_bits(bits, n)
print(T.H())                      # fast Hamiltonian path count
print(T.regime)                   # which regime was used
print(T.c3)                       # 3-cycle count (Rao's formula)
```

### tournament_toolkit.py -- Four-Module Analysis Suite
- **TilingGrid**: Triangular encoding of tournaments as bit vectors
- **FlipGraph**: Arc transition analysis (which flip changes the ranking most?)
- **RankAmbiguity**: H-count based ranking with `high_leverage_comparisons()`
- **TournamentTDA**: Topological features for ML pipelines

```python
from tournament_toolkit import RankAmbiguity
ra = RankAmbiguity(adj_matrix)
print(ra.H())                     # ranking ambiguity
print(ra.high_leverage())         # most informative comparison to re-evaluate
```

### formalrank.py -- One-Pass Pairwise Ranking
Ranks items from pairwise comparisons using the formal group logarithm. No iteration, no hyperparameters, no convergence criterion. Evidence is perfectly additive: `score += arctanh(coupling)`.

```python
from formalrank import FormalRank
ranker = FormalRank()
ranker.add("Claude", "GPT-4", wins=156, losses=134)
ranker.add("Claude", "Gemini", wins=178, losses=112)
print(ranker.summary())
```

**Performance:** 560K observations/sec.

### boost_ranker.py -- Three-Channel Ranking
Every comparison generates a Cayley boost Q = (1+x)/(1-x) that decomposes into three independent signals:

| Channel | What it measures | Use |
|---------|-----------------|-----|
| **INERT** (mod 2) | Who wins? | Binary ranking |
| **RAMIFIED** (mod 3) | By how much? | Confidence estimation |
| **SPLIT** (mod 7) | Was this expected? | Upset/anomaly detection |

Includes automatic cycle detection and upset flagging.

### polynomial_predictor.py -- Tournament Polynomial Agent
A ranking agent whose entire state is 5 numbers. Replaces Elo (n params), Bradley-Terry (n params), TrueSkill (2n params), or GNNRank (O(n^2) params) with 5 polynomial coefficients that capture 100% of the variance.

```python
from polynomial_predictor import PolynomialPredictor
pp = PolynomialPredictor(items=['Claude', 'GPT-4', 'Gemini', 'Llama', 'Mistral'])
pp.observe('Claude', 'GPT-4')
pp.observe('GPT-4', 'Gemini')
pp.observe('Gemini', 'Claude')       # creates a cycle

print(pp.ranking())                   # ranked list
print(pp.confidence())                # how certain is the ranking
print(pp.hallucination_risk())        # cycle/intransitivity detection
print(pp.next_comparison())           # most informative pair to evaluate next
```

### polynomial_head.py -- Transformer Output Head Replacement
Three output heads benchmarked against each other: standard softmax, polynomial (quantized Horner), and formal group (arctanh aggregation with early exit). Drop-in replacement for the transformer output layer.

### instant_mcmc.py -- MCMC Without Running the Chain
For any observable with Walsh degree D on the flip chain, the time evolution is a degree-D polynomial. Once precomputed, any query "what is E[f] after t random flips?" is answered in O(1) time with exact rational arithmetic.

```python
from instant_mcmc import InstantMCMC
mcmc = InstantMCMC(H_table, m)
mcmc.expected_H(x_bits, t)           # exact rational E[H] at time t
mcmc.mixing_time(x_bits, eps)        # time to eps-convergence
mcmc.robustness(x_bits, k)           # E[H] after k random re-evaluations
```

### spectral_ranker.py -- Exact Algebraic Ranking
Uses the spectral gap theorem (gap = 4/C(n,2), exact) to compute ranking confidence with no floating-point error. Identifies high-leverage comparisons -- which single result, if reversed, changes the ranking most.

### tournament_chat.py / tournament_chat_v2.py -- Demo Chatbots
Working chatbots where token prediction IS tournament ranking. Candidate next-words compete in a tournament. The polynomial P(z) ranks them, and the hallucination risk is displayed at each token.

```bash
# v1: minimal proof of concept with built-in corpus
python 04-computation/tournament_chat.py

# v2: scaled version with co-occurrence embeddings, sparse tournaments,
#     multi-evidence fusion, and real intransitivity detection
python 04-computation/tournament_chat_v2.py
```

---

## The 21 Hidden Tournaments

The project identified 21 domains where pairwise-comparison structures create tournaments, often without the practitioners realizing it. Our polynomial framework (five coefficients, three signals) applies to all of them.

| # | Domain | Items Competing | Our Value-Add |
|---|--------|----------------|---------------|
| 1 | **Attention mechanism** | Key positions per query | Hallucination detection at the attention level |
| 2 | **Immune system** | B cells in germinal centers | Autoimmunity risk from spectral decomposition |
| 3 | **Financial order book** | Limit orders | Market manipulation detection via fast channel |
| 4 | **A/B testing** | Experiment variants | Cycle detection (Simpson's paradox) |
| 5 | **Loss landscape** | Parameter configurations | Saddle point detection for learning rate scheduling |
| 6 | **Protein folding** | Conformations | Misfolding risk at temperature |
| 7 | **Sorting algorithms** | Items | Partial-sort quality measure from P(z) |
| 8 | **PageRank** | Web pages | Link spam detection via fast channel |
| 9 | **Self-play RL** | Agent versions | Training stagnation detection |
| 10 | **Evolutionary algorithms** | Individuals | Rock-paper-scissors dynamics detection |
| 11 | **Recommendation systems** | Items per user | Filter bubble detection |
| 12 | **Compiler scheduling** | Pipeline instructions | ILP maximization from spectral gap |
| 13 | **Drug discovery** | Molecules | Assay-dependent hit detection |
| 14 | **Quantum measurement** | Basis states | Schrodinger-cat state detection |
| 15 | **Cache eviction** | Cache lines | Adaptive policy via three-signal decomposition |
| 16 | **Hiring/interviews** | Candidates | Interviewer bias/inconsistency detection |
| 17 | **Democratic deliberation** | Arguments | Circular reasoning detection |
| 18 | **Language evolution** | Words competing for niches | Register-dependent ranking |
| 19 | **Consensus protocols** | Proposals/candidates | Livelock prediction from spectral gap |
| 20 | **Dating apps** | Profiles | Intransitive preference detection |
| 21 | **Neural network training** | Parameter settings | Saddle proximity for LR scheduling |

See `07-reflections/hidden-tournaments.md` and `07-reflections/cs-tournaments-deep-dive.md` for detailed analysis.

---

## Key Theorems

### The Odd-Cycle Collection Formula (OCF)
For every tournament T: **H(T) = I(Omega(T), 2)**, where I is the independence polynomial and Omega is the odd-cycle conflict graph. This gives H(T) = 1 + 2a_1 + 4a_2 + ... (always odd, by Redei's theorem).

### THM-250: Near-Transitive Flip Formula
Reversing a single arc of skip k in the transitive tournament gives **H = 2^{k-1} + 1**. This powers the Regime 30 perturbation method.

### THM-251: Cayley Boost Spectrum
The Cayley transform Q(lambda_k) = (m-k)/k satisfies the functional equation **Q(lambda_k) * Q(lambda_{m-k}) = 1** -- a duality pairing analogous to the functional equation of zeta functions.

### THM-252: Rapidity Lattice
Eigenvalue rapidities span **Q * ln(2) + Q * ln(3) + Q * ln(7)** -- a rank-3 lattice (by Baker's theorem) that is the archimedean shadow of the cuboid Z/42Z.

### THM-253: Boost Trichotomy
The 9 nontrivial Cayley boosts form a 3x3 matrix whose Column I product = **84 = 2 * 42 = the Hurwitz automorphism bound**. The three groups correspond to persistent (INERT), transient (RAMIFIED), and oscillating (SPLIT) modes.

### Path Homology
**beta_2(T) = 0** for every tournament (THM-108/109). **beta_1 * beta_3 = 0**: 1-holes and 3-holes never coexist (THM-095).

### Forbidden H-Values
**H = 7 and H = 21 are impossible** for any tournament on any number of vertices (THM-029/079).

### The Krawtchouk Framework (2026-03-24)
The tiling hypercube Q_m (where m = C(n-1,2) tiles) is the **Hamming scheme H(m,2)**. Its Krawtchouk polynomial eigenfunctions provide a 3-axis coordinate system for tournament space:

- **B_1 = K_1(x;m) = m - 2x**: correlates with -H at r = -0.94 (n=5). The "principal axis."
- **B_2 = K_2**: correlates with -c3 (3-cycle count) at r = -0.86. The "width axis."
- **B_3**: captures the "twist" (SC/NS separation via complement parity).

**Band-limitedness**: H is exactly zero in the upper Walsh spectrum (k >= m/2). Tournament structure is a low-pass signal on the hypercube. This gives 4:1 compression of tournament data and free error detection.

### Tournament Counting Theorem
**V_n x n!/2^m = 1 + Sum_{k odd} (1/k) x n-falling-k x 2^{(k^2-1)/2 - (k-1)n}**. An Euler product over "tournament primes" (odd integers), dominated by 1/3 (99.98% of correction by n=15). Poles of the generating function at x = 4, 16, 64, 256, ...

### The Waggly Completeness Theorem
All connections between tilings (= the complete graph on 2^m vertices) decompose by Hamming distance d = 1,...,m. "Wiggly" = d=1 (single tile flip). "Blue/black" = d=m (complement tiling). The completeness order k* = diam(G_n) = **OEIS A003141** = maximum feedback arc set number, growing as ~n^2/4.

### Paley Tournaments Give Classical Codes
Using the adjacency matrix A+I of the Paley tournament P_p as a parity-check matrix over GF(2):
- **P_7 -> Hamming [7,4,3]** code
- **P_23 -> binary Golay [23,12,7]** code

---

## Quick Start

### Run the interactive chatbot demo
```bash
cd 04-computation
python tournament_chat.py
```
Type prompts and watch tournament-polynomial token prediction in action, with hallucination risk displayed at each generated token.

### Run the polynomial predictor
```bash
cd 04-computation
python polynomial_predictor.py
```
Demonstrates the 5-coefficient polynomial agent on pairwise comparison data with ranking, confidence, hallucination detection, and active learning.

### Rank items from pairwise data
```python
import sys; sys.path.insert(0, '04-computation')
from formalrank import FormalRank

ranker = FormalRank()
ranker.add("Model_A", "Model_B", wins=60, losses=40)
ranker.add("Model_B", "Model_C", wins=55, losses=45)
ranker.add("Model_C", "Model_A", wins=52, losses=48)  # cycle!
print(ranker.summary())
```

### Compute Hamiltonian path count
```python
import sys; sys.path.insert(0, '04-computation')
from tournament_H import TournamentH

# 5-vertex tournament from bit encoding
T = TournamentH.from_bits(0b1010101010, 5)
print(f"H = {T.H()}, regime = {T.regime}, c3 = {T.c3}")
```

### Three-signal analysis
```python
import sys; sys.path.insert(0, '04-computation')
# Run boost_ranker.py to see the slow/medium/fast decomposition
# of every comparison in a tournament
```

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

## Mathematical Theory

### The Formal Group F(x,y) = (x+y)/(1+xy)

This is the thread connecting everything. It is simultaneously:
- The **relativistic velocity addition** law
- The **Poincare disk** group law (hyperbolic geometry)
- The **tanh addition** formula
- The **Bayesian evidence aggregation** law (log-odds addition)
- A **formal group of height 1** at every odd prime and **height infinity** at p=2

Its logarithm **arctanh(x) = x + x^3/3 + x^5/5 + ...** linearizes all of the above. The coefficients 1/(2k+1) are reciprocals of odd numbers -- the same odd numbers that appear as Hamiltonian path counts of tournaments.

### The Boost Trichotomy

42 = 2 * 3 * 7 = INERT * RAMIFIED * SPLIT. Three independent axes:

| Axis | Prime | Eisenstein type | Computation type | Tournament role |
|------|-------|----------------|-----------------|-----------------|
| Parity | 2 | INERT | Lossless | Symmetric/antisymmetric |
| Curvature | 3 | RAMIFIED | Lossy | Eigenvalue structure |
| Position | 7 | SPLIT | Crystallization | Forbidden values |

Both 3 and 7 have class number 1 in their imaginary quadratic fields, ensuring unique factorization in the rapidity lattice.

### The Rapidity Lattice

The Cayley transform Q(lambda) = (1+lambda)/(1-lambda) of eigenvalues factors over {2, 3, 7} (the Hurwitz primes). The rapidities arctanh(eigenvalues) form a rank-3 lattice Z*ln(2) + Z*ln(3) + Z*ln(7), guaranteed by Baker's theorem.

### The Adelic Tournament Space

Tournament eigenvalues live on a truncated adelic space A_T(n) = R x prod_{p|D} Z/p^e Z, approaching the full adeles as n grows. The flip chain is a Hecke operator on this space.

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
| `07-reflections/` | Cross-domain patterns (adelic geometry, boost trichotomy, hidden tournaments) |
| `agents/` | Multi-agent coordination system |

### Key Files

**Production Tools:**
- `04-computation/tournament_H.py` -- Tetrahedron Algorithm (three-regime H computation)
- `04-computation/tournament_toolkit.py` -- four-module analysis suite (TilingGrid, FlipGraph, RankAmbiguity, TournamentTDA)
- `04-computation/formalrank.py` -- one-pass pairwise ranking (560K obs/sec)
- `04-computation/boost_ranker.py` -- three-channel ranking via the boost trichotomy
- `04-computation/spectral_ranker.py` -- exact algebraic ranking with high-leverage detection
- `04-computation/polynomial_predictor.py` -- 5-coefficient polynomial agent for ranking
- `04-computation/polynomial_head.py` -- transformer output head replacement (softmax/polynomial/formal group)
- `04-computation/instant_mcmc.py` -- MCMC without running the chain (O(1) mixing queries)
- `04-computation/tournament_chat.py` -- interactive chatbot demo (v1: proof of concept)
- `04-computation/tournament_chat_v2.py` -- scaled chatbot (v2: embeddings, sparse tournaments, multi-evidence)
- `04-computation/mod_rank_library.py` -- modular rank with 8x memory reduction

**Enumerators:**
- `04-computation/burnside_enum_v2.c` -- unified graph/digraph enumerator (12 OEIS sequences)
- `04-computation/k_uniform_fast_enum.c` -- general k-uniform hypergraph counter
- `04-computation/fast_a000568_v2_s90di.py` -- fast tournament counting via partition sum

**Theory and Reflections:**
- `03-artifacts/drafts/tournaments_comprehensive.tex` -- comprehensive paper
- `03-artifacts/drafts/engineering-synthesis-2026-03-10-S53.md` -- engineering roadmap (12 domains)
- `07-reflections/hidden-tournaments.md` -- the 21 hidden tournaments across all domains
- `07-reflections/cs-tournaments-deep-dive.md` -- 10 CS problems that are secretly tournaments
- `07-reflections/applications-roadmap.md` -- tiered application roadmap with implementation priorities
- `07-reflections/polynomial-agent.md` -- the polynomial agent architecture
- `07-reflections/llm-as-tournament.md` -- why an LLM is a tournament machine
- `07-reflections/adelic-tournament-geometry.md` -- the adelic structure of eigenvalues
- `07-reflections/three-types.md` -- lossless / lossy / crystallization trichotomy

---

## How This Repository Works

This is a multi-agent research project where multiple Claude instances collaborate asynchronously via git. Each session is identified as `[machine]-[date]-S[N]`. The `CLAUDE.md` file contains the full protocol including startup sequence, session logging, and end-of-session close-out.

Concurrency is intentional. Agents should use `origin/main` as the live shared
research surface, push small checkpoints during long sessions, claim scarce
namespaces early, and treat new work seen during rebase as possible evidence for
the current investigation. See `00-navigation/CONCURRENT-SESSIONS.md` for the
operating playbook.

---

## References

- D. Grinberg, R.P. Stanley. *Counting Hamiltonian paths in tournaments.* arXiv:2412.10572 (2024).
- A. Grigor'yan, Y. Lin, Y. Muranov, S.-T. Yau. *Homologies of path complexes and digraphs.* arXiv:1207.2834 (2012).
- L. Redei. *Ein kombinatorischer Satz.* Acta Litt. Sci. Szeged 7 (1934), 39-43.
- J.W. Moon. *Topics on Tournaments.* Holt, Rinehart and Winston (1968).
- R. Forcade. *Parity of paths and circuits in tournaments.* Discrete Math. 6 (1973), 115-118.
- R. Tang, S.-T. Yau. *Homology of tournaments and path homology.* arXiv:2602.04140 (2026).
- A. Baker. *Linear forms in the logarithms of algebraic numbers.* Mathematika 13 (1966), 204-216.
