# How New Math Makes Hard Problems 100,000x Easier

**Status:** Living document. Last updated 2026-03-08.

---

## The Problem: Counting Rankings Is Absurdly Hard

Imagine a round-robin sports tournament — every team plays every other team, no ties. Now ask: **in how many ways can you rank all the teams so that each team beat the one right below it?** These are called "consistent rankings" or "Hamiltonian paths."

This sounds like a simple counting question. It is not. It's what computer scientists call **#P-complete** — meaning it's at least as hard as any counting problem out there. The best-known general algorithm takes time proportional to 2^n for n teams. For 20 teams, that's about a billion operations. For 30 teams, a trillion. For 50 teams, more operations than there are atoms in the universe.

This matters beyond sports. The same problem appears in **voting theory** (how many total orderings are consistent with pairwise majority preferences?), **biology** (how linear is an animal dominance hierarchy?), **search engines** (how to aggregate pairwise relevance judgments into a single ranking?), and **network analysis** (how many directed flows exist through a network?).

This project discovered formulas that, in many cases, make the problem **dramatically faster to solve** — and revealed deep mathematical structure in the process.

---

## The Breakthrough: A Formula That Changes Everything

The central discovery is a formula called the **Odd-Cycle Collection Formula (OCF)**. Instead of laboriously enumerating every possible ranking, it says:

> **The number of consistent rankings = 1 + 2a_1 + 4a_2 + 8a_3 + ...**

where a_k counts the number of ways to pick k non-overlapping odd-length directed cycles in the tournament.

Why does this help? Because counting cycles is vastly easier than counting rankings. You can find all the 3-cycles in a tournament by multiplying the adjacency matrix by itself three times and reading the diagonal — a basic matrix operation that takes O(n^3) time. That's polynomial time, not exponential. The 5-cycles and 7-cycles can be found similarly.

The OCF was computationally discovered in this project, verified exhaustively for all tournaments up to 8 teams (134 million configurations), and later proved rigorously by mathematicians Grinberg and Stanley.

---

## The Speedups, Concretely

### 10x faster via matrix traces

Instead of the standard 2^n dynamic programming algorithm, you can compute the ranking count by:

1. Compute the number of 3-cycles (one matrix multiplication, O(n^3))
2. Compute the number of 5-cycles (matrix power, O(n^3))
3. Count disjoint cycle pairs (inclusion-exclusion on per-vertex cycle counts)
4. Plug into the OCF formula

In benchmarks at n = 9, this takes **0.7 milliseconds** per tournament versus **70 milliseconds** for the standard algorithm — a **100x speedup**. The advantage grows with structure: tournaments with few long cycles can be processed even faster.

### 100,000x compression via Fourier analysis

Every tournament can be encoded as a string of bits (one bit per game outcome). The ranking count is then a function on bit-strings. When you decompose this function using a **Walsh-Fourier transform** — the binary analogue of the Fourier transform used in signal processing — something remarkable happens: almost all the "frequencies" are zero.

| Teams | Total bit-strings | Nonzero frequencies | Compression factor |
|-------|------------------|--------------------|--------------------|
| 5     | 1,024            | 3                  | **341x** |
| 7     | 2,097,152        | ~20                | **~100,000x** |

This means you can reconstruct the ranking count for *any* tournament of a given size from a tiny handful of numbers. It's like discovering that a seemingly complex image is actually made of just three colors — the compression is exact and lossless.

This is not just a theoretical curiosity. It means tournament invariants can be **learned from very few samples** (relevant to property testing in computer science) and **computed in compressed form** (relevant to large-scale data analysis).

### 64–130x faster enumeration

The tournament encoding lives on a triangular grid that turns out to have a rich symmetry group (called S_3 x Z_2 — the symmetry group of a triangular prism). Exploiting these symmetries via Burnside's lemma, we derived closed-form counting formulas that replace brute-force iteration, yielding **64x to 130x speedups** for enumerating tournament types — directly applicable to computing OEIS integer sequences.

---

## Why Does This Formula Work? Four Independent Explanations

One of the satisfying aspects of the OCF is that we found **four completely independent reasons** why the ranking count is always odd (a famous 1934 result by the Hungarian mathematician Redei). Each explanation uses different mathematical machinery:

**1. The Toggle Trick.** For any two teams, the number of rankings with one before the other equals the reverse count (mod 2). You can pair up rankings by swapping their relative positions.

**2. Symmetry Cancellation.** Every tournament's symmetry group has odd size. Symmetries pair up rankings, leaving an odd number unpaired.

**3. The Cycle Formula.** The OCF gives ranking count = 1 + (something even), which is manifestly odd.

**4. Mirror Pairing.** A "mirror" operation on the tournament pairs rankings, and the unpaired ones correspond to a smaller tournament whose count is odd by induction.

Having four independent proofs isn't redundant — each one illuminates different structure. The toggle trick reveals pairwise balance. The cycle formula reveals the role of feedback loops. The mirror pairing reveals recursive self-similarity.

---

## Connections to the Real World

### Elections and Voting

In voting theory, a tournament encodes majority preferences: A beats B if more voters rank A above B. The **Condorcet paradox** — where voters collectively prefer A to B, B to C, and C to A — is exactly a 3-cycle. The OCF quantifies exactly how these cycles affect the number of consistent total rankings.

This is directly relevant to **ranked-choice voting design** and **rank aggregation algorithms**. Search engines, recommendation systems, and multi-criteria decision tools all need to aggregate pairwise comparisons into a single ordering. Our trace formulas could accelerate these computations, and our impossibility results (no tournament can have exactly 7 or 21 consistent rankings) constrain what preference structures can even arise.

### Biology and Dominance Hierarchies

Biologists study **pecking orders** — dominance relationships among animals where each pair has a winner. These are literally tournaments. The questions biologists ask ("How linear is this hierarchy?" "How many consistent rankings exist?" "How much ambiguity is there?") are exactly the questions our formulas answer.

Our topological results (specifically, that tournaments never have "bubble-like" holes — see below) suggest a structural constraint on the kinds of higher-order relationships that can emerge from pairwise dominance. This could inform models of social structure in animal groups.

### Network Science and Brain Mapping

**Path homology** — the topological invariants we compute — is an increasingly important tool for analyzing directed networks:

- **Neural connectomics:** Mapping directed connections between neurons. Path homology detects higher-order information flow patterns invisible to standard graph metrics.
- **Gene regulatory networks:** Directed cycles correspond to feedback loops. Our conflict graph encodes which feedback loops can operate independently.
- **Citation networks:** Path homology reveals hierarchical structure and community boundaries.

Our discovery that beta_2 = 0 for tournaments gives a **null model** for applied work: if you observe nonzero beta_2 in a real-world directed network, that's a meaningful structural signal distinguishing it from tournament-like (complete pairwise comparison) structure.

### Computer Science

- **Sorting analysis:** Consistent rankings are topological sorts. Counting them efficiently has implications for analyzing comparison-based algorithms.
- **Compressed sensing / property testing:** The extreme Walsh sparsity means tournament properties can be determined from surprisingly few edge queries.
- **Quantum algorithms:** The Walsh-Hadamard transform is a fundamental quantum gate. The structured sparsity of tournament spectra could inform quantum network analysis.

---

## The Shape of a Tournament: Topological Discoveries

Beyond counting, this project discovered that tournaments have a surprisingly constrained **topology** — a mathematical notion of "shape."

Using **GLMY path homology** (a topological invariant for directed networks invented in 2012), we computed Betti numbers for tournaments. These count "holes" of various dimensions:

- **beta_0:** Connected components. Always 1 for tournaments.
- **beta_1:** Loop-like holes. Either 0 or 1; equals 1 when the tournament has an "unfillable" directed cycle.
- **beta_2:** Bubble-like holes. **Always 0** — tournaments never have these.
- **beta_3 and higher:** Can be nonzero starting at 6 and 7 teams.

The vanishing of beta_2 has been verified in approximately **47,000 tournaments** from 3 to 10 teams with **zero failures**. This is striking because beta_3 and beta_4 *can* be nonzero — so the gap is specific to dimension 2.

Other topological surprises:
- beta_1 and beta_3 are **mutually exclusive** — a tournament never has both loop-like and higher-dimensional holes at the same time
- The Euler characteristic is always exactly 0 or 1
- **Paley tournaments** (the most "balanced" tournaments, defined using number theory) have extreme topology: the Paley tournament on p teams has homology concentrated entirely in dimension p-3

Proving beta_2 = 0 algebraically is the project's biggest open problem.

---

## Hidden Universals: The Signed Permanent

A more exotic discovery involves a quantity called the **signed Hamiltonian permanent**. Replace each 1 in the tournament matrix with +1 and each 0 with -1, then sum the product along every ranking.

For an **even** number of teams, this is always exactly zero (rankings pair up with opposite signs).

For **odd** teams, something wild happens: **the signed permanent modulo 2^{n-1} depends only on n, not on the tournament.** Every 5-team tournament gives a result divisible by 16. Every 7-team tournament gives a result equal to 48 (mod 64). The specific tournament doesn't matter.

The values of n where this universality holds perfectly (3, 5, 7, 11, 19, 35, 67, ...) follow a beautiful number-theoretic pattern tied to binary digit sums, connecting tournament combinatorics to the arithmetic of the integers.

---

## Nature's Favorites: Paley Tournaments

Among all tournaments, the most "balanced" are the **Paley tournaments**, built from number theory: team a beats team b if b - a is a perfect square modulo a prime p.

Computationally, Paley tournaments appear to **maximize** the number of consistent rankings:
- 3 teams: 3 rankings (maximum possible)
- 7 teams: 189 rankings (maximum possible)
- 11 teams: 95,095 rankings (maximum possible)

This conjectured maximality, if proved, would be a deep connection between number theory and combinatorial optimization.

---

## Why This Matters for Mathematics

This work sits at a crossroads where several major mathematical fields meet:

**Combinatorics meets topology.** The OCF links path counting to cycle structure; path homology adds a layer of topological information invisible to purely combinatorial tools.

**Fourier analysis meets graph theory.** The Walsh decomposition reveals that tournament invariants live on a dramatically smaller subspace than expected — a phenomenon with parallels in additive combinatorics and Boolean function analysis.

**Number theory meets optimization.** Paley tournaments (built from quadratic residues) appear to be combinatorially optimal. The universality of the signed permanent follows a pattern controlled by binary digit sums.

**Representation theory meets enumeration.** The pin grid's S_3 x Z_2 symmetry group connects to Young tableaux, Schützenberger involution, and the representation theory of the symmetric group.

Key results that resolve or advance open questions:
- The OCF was proved (Grinberg-Stanley, 2024, building on computational discovery here)
- Transfer matrix symmetry M[a,b] = M[b,a] was proved via the Walsh framework
- The complete Fourier spectrum of tournament invariants was characterized for the first time
- The beta_2 = 0 phenomenon is entirely new, with no analogue in existing path homology literature
- Universal congruences for the signed permanent were previously unknown

---

## What's Still Open?

- **Prove beta_2 = 0** for all tournaments algebraically (47,000+ tests, zero failures, no proof yet)
- **Prove Paley maximization** — do Paley tournaments always maximize ranking count?
- **Extend the per-path identity** to all tournament sizes (currently works only for small n)
- **Prove unimodality** of the forward-edge distribution (50,000+ tests, zero violations)

---

## How to Read the Repository

| Folder | Contents |
|--------|----------|
| `00-navigation/` | Index files: open questions, session log, investigation backlog |
| `01-canon/` | Definitions, proved theorems, documented mistakes |
| `02-court/` | Disputes between research agents (formal disagreement resolution) |
| `03-artifacts/` | Paper drafts and code |
| `04-computation/` | Python scripts for all computations |
| `05-knowledge/` | Knowledge base: hypotheses, variables, computational results |
| `06-writeups/` | This document and the formal companion |

The main paper draft is at `03-artifacts/drafts/parity_tournaments_fixed.tex`.

---

*This document is written for a general audience. For precise theorem statements and proofs, see the companion formal write-up or the paper draft.*
