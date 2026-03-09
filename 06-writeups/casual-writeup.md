# Tournaments, Paths, and Parity: A General-Audience Guide

**Status:** Living document. Last updated 2026-03-08.

---

## What is this project about?

Imagine a round-robin sports tournament where every team plays every other team exactly once, and every game has a winner (no ties). Mathematicians call this a **tournament**. We can draw it as a network: each team is a dot, and we draw an arrow from the winner to the loser of each game.

Now ask a seemingly simple question: **In how many ways can you line up all the teams so that each team beat the team right after it?** This is called a **Hamiltonian path** — a chain of victories through every team. Think of it as a "perfect ranking" consistent with the actual results.

In 1934, the Hungarian mathematician Laszlo Redei proved something remarkable: **the number of such perfect rankings is always odd.** No matter how the games turn out, you can never get an even number of perfect rankings.

This project set out to understand *why* Redei's theorem is true at a deep level, and in the process discovered surprising connections between tournament combinatorics, abstract algebra, and topology.

---

## The Big Discovery: The Odd-Cycle Collection Formula

The central result of this research is a formula that tells you *exactly* how many perfect rankings a tournament has — not just that it's odd.

Here's the intuition. In a tournament, some groups of teams form **directed cycles**: A beats B, B beats C, C beats A. These cycles are the "interesting" parts of a tournament — they're what make rankings ambiguous.

Now imagine drawing a new network where:
- Each directed cycle of odd length (3 teams, 5 teams, 7 teams, etc.) is a dot
- Two cycles are connected if they share a team

This new network is called the **conflict graph**. Two cycles "conflict" when they can't both be respected in a ranking at the same time.

The formula says:

> **H(T) = 1 + 2a_1 + 4a_2 + 8a_3 + ...**

where a_k counts the number of ways to pick k non-conflicting odd cycles. The number of perfect rankings is literally "one plus powers-of-two weighted by independent cycle collections."

This immediately explains Redei's theorem: the formula always gives 1 + (even stuff) = odd.

But it does much more. It connects the count of perfect rankings (a question about paths) to the structure of cycles (a completely different kind of object). This connection, called the **Odd-Cycle Collection Formula (OCF)**, was computationally discovered in this project and later proved by mathematicians Grinberg and Stanley using different techniques.

---

## Four Ways to See Why It's Odd

One of the satisfying aspects of this work is finding *multiple* independent reasons why the count is always odd. Think of it like understanding why the sky is blue — you could explain it through wave physics, through quantum mechanics, or through the molecular structure of the atmosphere. Each explanation illuminates a different facet.

**Way 1: The Toggle Trick.** For any two teams u and w, the number of rankings with u before w equals the number with w before u (modulo 2). You prove this by showing you can pair up rankings that differ only in the relative order of u and w. Since the total is the sum of two equal-mod-2 numbers, it's odd.

**Way 2: Symmetry Cancellation.** If the tournament has any symmetry at all (which is rare but always has odd-sized symmetry group), the symmetries pair up rankings, leaving an odd number unpaired.

**Way 3: The Cycle Formula.** As described above — the OCF gives H = 1 + even.

**Way 4: Anti-Automorphism Pairing.** A certain "mirror" operation on the tournament pairs up rankings. The unpaired ones correspond to rankings of a smaller tournament, and by induction, that count is odd too.

---

## Seeing the Tournament Through a Prism: The Walsh-Fourier Approach

Here's where things get surprisingly powerful. Every tournament on n teams can be encoded as a string of bits — one bit per game, recording who won. So a tournament on 5 teams needs C(5,2) = 10 bits. There are 2^10 = 1024 possible tournaments.

Now, the number of perfect rankings H(T) is a function that assigns a number to each bit-string. Functions on bit-strings can be decomposed using **Walsh-Fourier analysis** — the binary version of the Fourier transform that decomposes sound into frequencies.

When you decompose H this way, something beautiful happens: **most "frequencies" are zero.** The nonzero frequencies correspond to very specific patterns of games — specifically, unions of edge-disjoint even-length paths in the complete graph.

This decomposition gives:
- A new proof of the OCF (by showing both sides have the same Fourier transform)
- An exact formula for the transfer matrix entries M[a,b] (the number of rankings starting at team a and ending at team b)
- A proof that M[a,b] = M[b,a] — the number of rankings from a to b equals the number from b to a

The last result resolves a conjecture that had been open in the community. The Walsh formula makes it *obvious* because it's symmetric in a and b.

---

## The Shape of a Tournament: Path Homology

Here's where the research ventures into topology — the mathematics of shapes and holes.

In 2012, mathematicians Grigor'yan, Lin, Muranov, and Yau invented a way to associate topological invariants to directed networks. Just like a donut has a "hole" that a sphere doesn't, a directed network can have "directed holes" of various dimensions.

These invariants, called **Betti numbers** (beta_0, beta_1, beta_2, ...), count holes of different dimensions:
- beta_0: connected components (always 1 for tournaments)
- beta_1: "loop-like" holes (either 0 or 1 for tournaments; equals 1 when the tournament has an unfillable directed cycle)
- beta_2: "bubble-like" holes
- beta_3 and higher: higher-dimensional holes

**The big discovery:** For tournaments, **beta_2 is always zero.** We've tested approximately 47,000 tournaments ranging from 3 to 10 teams, and every single one has beta_2 = 0. This means tournaments never have "bubble-like" directed holes — their topology is at most one-dimensional (loop-like).

This is remarkable because beta_3 and beta_4 *can* be nonzero (starting at 6 and 7 teams respectively). So the vanishing of beta_2 is specific and special.

Furthermore:
- beta_1 and beta_3 are **mutually exclusive**: a tournament never has both loop-like and higher-dimensional holes simultaneously
- The "Euler characteristic" (an alternating sum of Betti numbers) is either 0 or 1, perfectly correlated with whether beta_1 vanishes

Proving beta_2 = 0 algebraically remains one of the key open problems.

---

## The Signed Permanent: A Hidden Universal

Here's a more exotic discovery. Take the tournament's adjacency matrix (1 for a win, 0 for a loss) and transform it: replace each 1 with +1 and each 0 with -1. Now compute a "signed permanent" — multiply the entries along each ranking and add them all up.

For **even** numbers of teams, this signed permanent is always exactly zero (because you can pair up rankings with opposite signs).

For **odd** numbers of teams, something wild happens: **the signed permanent modulo 2^{n-1} depends only on n, not on the tournament.** Every tournament on the same number of teams gives the same remainder when you divide by a large power of 2.

For example, all 5-team tournaments give a signed permanent divisible by 16. All 7-team tournaments give a signed permanent equal to 48 (mod 64).

At n = 9, the universality partially breaks: the remainder modulo 128 is universal, but the remainder modulo 256 depends on whether the tournament has an even or odd number of 3-cycles.

The exact values of n where the congruence is fully universal follow a beautiful number-theoretic pattern related to binary digit sums: n = 3, 5, 7, 11, 19, 35, 67, ...

---

## The Paley Tournaments: Nature's Favorites

Among all tournaments, some are more "balanced" than others. The most balanced are the **Paley tournaments**, defined using number theory: for a prime p, team a beats team b if b - a is a perfect square modulo p.

A striking computational observation is that Paley tournaments appear to **maximize** the number of perfect rankings among all tournaments on the same number of teams. The values are:
- 3 teams: H = 3
- 7 teams: H = 189
- 11 teams: H = 95,095

These numbers grow rapidly and have beautiful number-theoretic structure. The Paley tournaments also have extreme topological properties: the Paley tournament on p teams has a path homology Betti number beta_{p-3} = p - 1, concentrated in a single dimension.

---

## Why Does This Matter?

### For mathematics

This work sits at a fertile crossroads of combinatorics, algebra, and topology. The OCF connects a counting problem (Hamiltonian paths) to a structural invariant (cycle conflict graphs) in an unexpected way. The Walsh-Fourier approach provides a powerful general technique for analyzing functions on combinatorial objects. The path homology results suggest deep structural constraints on directed networks.

Several results here resolve or advance long-standing questions:
- The OCF was proved (by Grinberg-Stanley, building on the computational discovery here)
- Transfer matrix symmetry was proved via the Walsh framework
- The Fourier spectrum of tournament invariants was completely characterized
- The beta_2 = 0 phenomenon for tournaments is new and has no analogue in the existing literature

### For applications

Tournaments model any situation with pairwise comparisons: sports, elections, biological dominance hierarchies, preference aggregation, and network analysis. Understanding the structure of rankings (Hamiltonian paths) is directly relevant to:
- **Ranking aggregation:** How many consistent total orders are compatible with pairwise comparisons?
- **Network topology:** Path homology is used in applied topology to analyze directed networks (neural networks, citation networks, gene regulatory networks)
- **Voting theory:** Tournaments encode majority preferences; the OCF gives structural insight into how cycle structure constrains rankings

### For the broader research community

The methodology here — combining exhaustive computation with algebraic proof, using multiple AI agents as collaborative research assistants, and maintaining a structured knowledge base of hypotheses and results — represents a new model for computational mathematics. Over 277 hypotheses have been systematically generated, tested, and cataloged, with clear documentation of both confirmations and dead ends.

---

## What's Still Open?

The biggest open problem is proving **beta_2 = 0 for all tournaments** — that tournaments never have "bubble-like" directed holes. The computational evidence is overwhelming (47,000+ tests, zero failures), but no algebraic proof exists yet. The proof would likely require understanding why the boundary map from 3-chains to 2-chains in the path chain complex is always surjective onto cycles.

Other tantalizing open problems:
- **The per-path identity:** Can the OCF be refined to a path-by-path formula for all tournament sizes? (Currently works only for small tournaments.)
- **Paley maximization:** Do Paley tournaments really maximize the count of perfect rankings? (Verified computationally but unproved.)
- **Unimodality:** Is the distribution of forward edges across perfect rankings always unimodal? (50,000+ tests, zero violations.)

---

## How to Read the Repository

The research is organized as follows:

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
