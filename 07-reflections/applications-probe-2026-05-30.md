# Applications Probe: What This Work Might Unlock

**Session:** codex-2026-05-30-applications-probe  
**Status:** speculative search note, not canon  
**Prompt:** search for possible applications of the current tournament/OCF/residue
toolkit; be creative, include mathematical problems, practical uses, and speedups.

## The Reusable Machinery

The repo now has several pieces that are unusually portable:

```text
H(T)                         exact count of consistent total orderings
OCF                          H(T) = I(Omega(T), 2)
Omega(T)                     conflict graph of odd directed cycles
alpha-vector                 disjoint odd-cycle packing profile
Delta H / deletion residue   most informative comparison or support
formal group arctanh         one-pass additive evidence aggregation
Krawtchouk/Walsh structure    low-frequency compression on tournament bits
path-homology constraints     beta_2 = 0 and beta_1/beta_3 separation
residue calculus             choose supports -> project/forget -> inspect survivor
```

The shared application pattern is:

```text
pairwise decisions create a tournament
cycles are not noise by default
H measures residual ordering freedom
Omega identifies the independent packets of conflict
residue features identify which comparison, vertex, or projection is carrying the obstruction
```

## Best Near-Term Application Bets

### 1. Active Pairwise Reranking For RAG And LLM Evaluation

External work on pairwise ranking prompting is now explicitly framing reranking
as active learning from noisy, order-sensitive, sometimes intransitive pairwise
comparisons. That is almost exactly the local `RankAmbiguity` problem.

Local upgrade:

- Use `H(T)` as the stopping statistic: stop querying when the top-k order has
  low residual Hamiltonian-path mass, not merely when a sorting algorithm has
  spent its budget.
- Use `Delta H` or expected `I(Omega,2)` drop as the acquisition function:
  query the pair that destroys the most independent odd-cycle mass.
- Use forbidden or near-forbidden scalar residues as a data-quality alarm:
  if a weighted/partial comparison system realizes hard-tournament-forbidden
  ambiguity values, the dataset is not behaving like a complete binary
  tournament.

Practical product shape:

```text
input: candidate docs/models/items and pairwise oracle
state: weighted tournament plus confidence intervals
query rule: maximize expected drop in H or top-k H
output: top-k ranking, ambiguity mass, named contradiction cycles
```

Why this might be a real speedup: active PRP wants better quality per LLM call.
The current active learners focus on winner recovery or noisy rank aggregation.
OCF can target the exact obstruction to a stable ranking: independent odd-cycle
packets. If most real reranking sets are near-transitive, the perturbation
regime should avoid most all-pairs comparisons.

First experiment: wrap a BEIR-style reranking benchmark with a tiny `H_topk`
tracker for `k <= 8`. Compare "query by smallest margin" against "query by
expected H drop" under the same call budget.

### 2. Lossless A/B Testing And Experiment Triage

Most experimentation systems collapse pairwise outcomes into scalar lift or
posterior probability. That destroys cycles such as:

```text
A beats B on retention
B beats C on conversion
C beats A on revenue quality
```

Local upgrade:

- Treat variants as vertices and metric-conditioned wins as tournament edges.
- Report `H`, `alpha_1`, `alpha_2`, and named cycles instead of a single winner.
- Use the high-leverage edge as the next experiment to rerun, segment, or power.

Wild but plausible theorem target: in multi-metric A/B systems, Simpson-style
paradoxes should appear as rank-2 Omega residues before they appear in global
lift aggregates. The THM-025 near-kill pattern is the model: a tiny survivor
can carry the real obstruction.

### 3. Clinical Or Policy Network Meta-Analysis

Network meta-analysis compares treatments through direct and indirect pairwise
evidence. Inconsistency loops are already a core issue. OCF gives a sharper
language:

```text
cycle = irreducible tradeoff or evidence inconsistency
independent cycles = multiple non-overlapping inconsistency packets
H = number of viable treatment orderings
Delta H = next comparison/trial with maximum ordering value
```

This should be framed as an audit tool, not an automatic clinical decision
tool. The value is that it names exactly which evidence loop keeps guidelines
ambiguous.

First experiment: take a public network meta-analysis graph, threshold each
pairwise comparison into a tournament by effect direction, then compute
cycle-packet features. Compare them with standard inconsistency diagnostics.

### 4. Attention And Decoding Diagnostics

The repo's Cartan bridge suggests splitting each attention or logit interaction
matrix into:

```text
antisymmetric part = tournament / competition
symmetric part     = metric / self-knowledge / compatibility
commutator         = entanglement between competition and metric confidence
```

Hypothesis: hallucination is not just low top-token margin. It is a mismatch
between a decisive antisymmetric tournament and a confused symmetric metric, or
the reverse: high symmetric confidence with cyclic antisymmetric preferences.

Practical speedups:

- adaptive depth: early-exit only when both the logit gap and Cartan
  consistency are high;
- speculative decoding filter: reject branches whose top-k token tournament
  has rising fast-channel cycle mass;
- calibration: use cycle mass as a zero-parameter confidence feature before
  adding learned probes.

First experiment: for a small local transformer, log per-token `(gap,
intransitivity, Cartan commutator norm)` on factual QA versus creative
continuation. The target is separation, not a trained classifier.

### 5. Tournament Fingerprints For Comparison Databases

Pairwise-comparison datasets are expensive to store, index, deduplicate, and
compare. The Krawtchouk/Walsh story says tournament structure is low-frequency
on the tiling cube.

Application:

- store a compact fingerprint: score sequence, `H`, `c3`, SCC defect, residue
  profile, low Krawtchouk coordinates;
- use it as a fast pre-screen before canonical isomorphism or expensive
  similarity search;
- treat high-frequency Krawtchouk energy as an error-detection checksum for
  corrupted or partial comparison data.

This is not only compression. It is a diagnostic index: a sports season, model
leaderboard, product preference panel, or market microstructure snapshot can be
queried by cycle profile instead of raw edge table.

First experiment: benchmark the existing `tournament_fingerprint.py` against
canonicalization for `n=7,8` representatives, then add residue-rank features
from `tournament_tda.py`.

## Best Pure-Math Targets

### A. Hamiltonian-Path Maximizers Beyond Szele-Alon

Alon's theorem puts the maximum number of Hamiltonian paths within a polynomial
factor of the random mean. The repo has more detailed structure:

- `H(T)=I(Omega,2)` separates total odd-cycle mass from disjoint packing mass.
- Paley wins the small prime cases `p=3,7,11`.
- Interval/cyclic constructions plausibly overtake Paley among larger
  circulants because higher `alpha_k` terms dominate.
- Krawtchouk coordinates give a compressed search basis for extremizers.

Possible result: prove an explicit Paley-to-interval phase transition among
circulant tournaments by expressing the first few `alpha_k` terms through
additive energy or Fejer-kernel quantities.

This is the cleanest path from repo computations to a classical theorem-shaped
paper.

### B. Which Odd Integers Are H-Values?

The forbidden values `7` and `21` are no longer isolated curiosities. The H=63
unlock shows a complete-Omega single-core mechanism:

```text
H = 1 + 2*r_core(signature)
```

and the image gaps `r_core != 3,10` look like the exact complete-core
explanation for `H != 7,21`.

Possible result: a finite-state or recurrence proof classifying the image of
single-core complete-Omega signatures. More ambitious: prove that every
permanent H-forbidden value is caused by a finite set of impossible residue
shapes in the OCF alpha-vector.

First experiment: convert weighted signatures into an automaton modulo powers
of 2 and 3. The immediate target is still `3` and `10`, not all odd values.

### C. Paley/Circulant Path-Homology Formula

Tang-Yau's circulant path-homology paper uses shift symmetry and Fourier
decomposition to reduce GLMY boundary-rank computations to finite symbol
matrices. This directly matches the repo's Paley/circulant machinery.

Possible result: closed formulas for low Betti numbers of Paley tournaments,
or a general even-Betti vanishing range for tournaments. The universal
`beta_2=0` theorem is the base case; circulants are the first family where a
closed formula feels plausible.

First experiment: implement the symbol-matrix recipe for Paley connection sets
`p=7,11,19,23`, then compare with existing `circulant_homology.py` data.

### D. Real-Root Failure Classification

Universal real-rootedness is false at `n=9`, but the known failure is
structured: a near-kill deletion residue with small rank-2 survivor.

Possible result: classify first failures by deletion-residue rank and prove a
weaker ultra-log-concavity statement for tournament conflict graphs even after
real-rootedness dies.

First experiment: exhaust or heavily sample n=9 non-real-root examples and log
`max deletion loss`, `rank_res`, alpha residue, score sequence, and SCC defect.

## Wilder Application Hypotheses

### 6. Compiler And Build Scheduling

Items are instructions, tasks, or build targets. Edge `i -> j` means "i should
precede j" according to a local heuristic. Real schedulers have cycles because
latency, cache locality, register pressure, and parallelism disagree.

Use `H` as scheduling freedom, Omega cycles as irreducible heuristic conflicts,
and `Delta H` as the next benchmark or constraint to resolve. A near-transitive
schedule should fall into the perturbative regime, giving a shortcut over
generic feedback-arc search.

### 7. Cache Eviction And Memory Tiering

Cache candidates compete by recency, frequency, write cost, object size, and
future-likelihood signals. A scalar policy collapses these signals. A
tournament policy preserves cycles:

```text
A beats B by recency
B beats C by frequency
C beats A by size or write penalty
```

Cycle mass becomes a phase-change detector. When fast-channel cycle mass rises,
the workload is changing and the cache should switch policy or increase
exploration.

### 8. Markets And Wash-Trading Rings

Order flow naturally creates directed pairwise relations: who sells to whom,
whose quotes dominate whose, which venue leads another. Independent directed
cycles are plausible signatures of rings, latency arbitrage loops, or synthetic
self-crossing. The advantage over generic graph cycle detection is the
ranking-pressure interpretation: `H` says how much ordering ambiguity the cycle
structure creates.

### 9. Consensus And Leader Election

Raft/Paxos-style systems usually hide pairwise preference cycles inside
timeouts and randomized elections. Build a tournament among candidate leaders
from freshness, latency, term, and quorum evidence. If `H` is high or Omega has
independent cycles, the system is in a livelock-prone regime. Timeouts should
respond to measured cycle temperature, not a fixed constant.

### 10. Protein Or Molecule Conformer Selection

Pairwise "A is lower energy than B under condition c" comparisons often depend
on solvent, temperature, scoring model, or assay. If those contexts disagree,
the result is a tournament with real cycles, not just noise. OCF could separate
one metastable tradeoff from several independent metastable packets.

This is speculative, but it has the right shape: small top-k candidate sets,
expensive pairwise evidence, and genuine context-dependent intransitivity.

### 11. Legal, Debate, And Negotiation Strategy

Arguments often beat each other pairwise by audience or criterion:

```text
A beats B on precedent
B beats C on equity
C beats A on administrability
```

`H` measures how many narrative orderings remain viable. Omega cycles identify
the irreducible tradeoffs. `Delta H` suggests which factual stipulation or
concession would collapse the most ambiguity.

### 12. Curriculum And Test Design

Problems, skills, or lessons can form a prerequisite tournament from pairwise
"should come before" judgments. Cycles reveal conceptual mutual dependencies.
Independent cycles are separable curriculum bottlenecks. A high-leverage edge
is the prerequisite relation worth testing with students.

## Speedup Bets

### Speedup 1: H-Drop Active Learning Beats Margin Sampling

Margin sampling chooses uncertain pairwise edges. H-drop sampling chooses edges
that collapse many Hamiltonian paths. These are not the same: a moderately
uncertain edge inside many independent cycles can be more valuable than the
least certain edge in an isolated triangle.

Benchmark target: PRP reranking, Chatbot Arena subsets, synthetic Bradley-Terry
data with injected cycles.

### Speedup 2: Near-Transitive Perturbation Is The Default Real-World Case

Most practical comparison systems have a rough latent ranking plus a small
number of upsets. That is Regime 30. If the upset count is `k`, the local
method is closer to `O(k*2^k)` than Held-Karp `O(n^2*2^n)`.

Benchmark target: sports seasons, leaderboard snapshots, A/B variants, and
RAG top-k reranking sets.

### Speedup 3: Residue Features Pre-Filter Expensive Algebra

Before computing full Omega, path homology, or canonical classes, compute:

```text
score sequence
c3
SCC defect
max deletion loss
small residue alpha/rank
low Krawtchouk coordinates
```

These should separate most mundane cases from exact kills, near kills,
Paley-flat families, and real-root failure candidates.

### Speedup 4: Instant MCMC For Robustness Audits

If an observable has bounded Walsh degree, its expected value after random edge
flips is a finite polynomial in time. That means robustness to label noise can
be answered without running the Markov chain.

Application: "How stable is this ranking if 3 percent of pairwise judgments are
wrong?" Precompute once, answer many perturbation queries in O(1).

## What I Would Actually Do Next

1. **Build an OCF active reranker prototype.** Use `H_topk` for `k <= 8`,
   query by expected H-drop, compare against margin and random on a public PRP
   reranking task.

2. **Run an n=9 real-root failure census.** This is the best small pure-math
   experiment because the residue-rank hypothesis is already concrete.

3. **Implement Paley symbol matrices.** This connects directly to new external
   path-homology language and could produce a clean closed-family theorem.

4. **Make `tournament_tda.py` output an application-grade JSON report.** Fields:
   ranking ambiguity, named 3-cycles, SCC defect, residue rank, high-leverage
   comparisons, and a "near-transitive/regime" tag.

5. **Apply the report to one real dataset.** The strongest demos are RAG
   reranking, sports, or A/B tests because the input format is just pairwise
   comparisons.

## Source Trail

Local sources read:

- `README.md`
- `00-navigation/CONCEPT-MAP.md`
- `00-navigation/OPEN-QUESTIONS.md`
- `07-reflections/applications-roadmap.md`
- `07-reflections/applications-vision.md`
- `07-reflections/breakthrough-candidate-map.md`
- `07-reflections/residue-calculus-theses.md`
- `07-reflections/tournament-compression-and-beyond.md`
- `07-reflections/fractal-tournament-codec.md`
- `04-computation/tournament_toolkit/README.md`
- `04-computation/tournament_tda.py`

External primary/source pages checked:

- Grinberg and Stanley, Redei-Berge symmetric function page:
  <https://scholarship.miami.edu/esploro/outputs/preprint/The-Redei--Berge-symmetric-function-of-a/991032796911302976>
- Mitrovic, noncommuting Redei-Berge function:
  <https://arxiv.org/abs/2504.20968>
- Tang and Yau, path homology of circulant digraphs:
  <https://arxiv.org/abs/2602.04140>
- Active Learners as Efficient PRP Rerankers:
  <https://arxiv.org/abs/2605.14236>
- PARWiS active pairwise comparisons:
  <https://arxiv.org/abs/2603.01171>
- Alon, maximum number of Hamiltonian paths in tournaments:
  <https://web.math.princeton.edu/~nalon/PDFS/hamilton.pdf>

