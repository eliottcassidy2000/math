---
id: HYP-9462
title: "Two Thomassen-merge conjectures. (A) The trivial group and {7, 21}: the Hamiltonian-path counts of asymmetric tournaments (Aut = 1) are exactly the odd numbers other than 3, 7 and 21. (B) Hypohamiltonian side: in a uniquely Hamiltonian GRAPH no two consecutive vertices of the Hamiltonian cycle are both non-critical (G - v Hamiltonian), the opposite of the single run of THM-4868 (C) for tournaments. Pre-registered tests: (A) at H in [225, 1125], (B) exhaustive at n = 8 and random at n = 12-15"
status: PRE-REGISTERED (kind-pasteur-2026-11-04-S10); confirmatory runs not yet made at the time of this commit
session: kind-pasteur-2026-11-04-S10
date: 2026-11-04
related: [THM-4868, THM-1370, THM-4810, THM-4812, THM-4747, THM-4743, THM-4851, THM-4505, THM-4720, THM-4721, HYP-9340, HYP-9420, HYP-9459, HYP-9461]
note: 05-knowledge/results/thomassen_merge_hypohamiltonian_groups_20261104.md (to be written after the runs)
script: 04-computation/experiments/asym_spectrum_closure_20261104.py, 04-computation/experiments/uh_graph_noncritical_20261104.c
---

# HYP-9462 - the asymmetric path spectrum and the non-critical vertices of uniquely Hamiltonian graphs

## (A) The trivial group and {7, 21}

H(T) is the number of Hamiltonian paths, and Aut(T) the automorphism group, which always has odd order.

**Known in the repository.**
- THM-1370: H is never 7 or 21.
- THM-4810 / HYP-9420: every other odd number up to 6000 is attained, and Chan-Pak completeness is OPEN.
- THM-4747: |Aut(T)| divides H(T), since Aut acts freely on Hamiltonian paths.
- THM-4812: H / |Aut| is never 7 or 21.

**Conjecture (A).** The set {H(T) : Aut(T) = 1} is exactly the set of odd numbers other than 3, 7 and 21.

**Proved.** 3 is excluded. H(T) = 3 forces a strong component with H = 3, which is the 3-cycle, so Aut(T) contains Z_3.

**Exploratory evidence** (FINITE-EXACT, `thomassen_merge_20261104_explore.out`, made before this commit):
- On H <= 224 the asymmetric holes are exactly {3, 7, 21}. This range is complete: strong tournaments of order >= 11
  have H >= 225.
- Among strong classes with H <= 224, only 3 needs a nontrivial automorphism group. The |Aut| values seen are
  1, 3, 5, 7, 9, 21, and the largest belongs to the Paley tournament P_7 (H = 189). |Aut| | H holds throughout.
- Uniquely Hamiltonian tournaments are all asymmetric (THM-4868 (E)). Their path counts at orders 4-14 (complete
  below 1125) generate, under products, every odd number up to 1125 except 24 values: 3, 7, 13, 21, 41, 43, 67, 77,
  83, 163, 173, 231, 233, 247, 251, 277, 299, 421, 449, 503, 533, 569, 637, 821.

**Pre-registered prediction (A1).** Take two inputs:
- every strong class of order <= 10 with H <= 1125 and |Aut| = 1 (`asym_spectrum_autspec_20261104 autspec 10 1125`);
- the uniquely Hamiltonian path counts of orders 4-14.

The product closure of their H values covers every odd number in [225, 1125]. This is a sufficient construction: a
miss would refute A1, not (A) itself, since strong asymmetric classes of orders 11-13 also have H <= 1125.

## (B) Non-critical vertices of uniquely Hamiltonian graphs

A graph is uniquely Hamiltonian (UH) if it has exactly one Hamiltonian cycle C. A vertex v is non-critical if G - v is
Hamiltonian.
- For UH tournaments, the non-critical vertices form ONE run of C (THM-4868 (C)).
- For hypohamiltonian graphs, every vertex is non-critical and G has no Hamiltonian cycle (Thomassen's planar examples
  on 105 vertices, 1976; the smallest known has 40).

**Conjecture (B).** In every uniquely Hamiltonian graph, no two consecutive vertices of C are both non-critical.

**Proved partial case.** Write C = c_0 c_1 ... c_(n-1), and suppose c_0 and c_1 are both non-critical, with Z_i a
Hamiltonian cycle of G - c_i.
- If c_i is adjacent to two Z_i-consecutive vertices, inserting it gives a Hamiltonian cycle of G, which must be C.
  So either Z_0 = C - c_0 + c_(n-1)c_1 and Z_1 = C - c_1 + c_0c_2 (the shortcut forms), or N(c_i) contains no two
  Z_i-consecutive vertices.
- If both Z_i are shortcuts, then c_(n-1) -> c_1 -> c_0 -> c_2 -> ... -> c_(n-1) is a second Hamiltonian cycle,
  which is a contradiction. The remaining case is open.

**Exploratory evidence** (before this commit):
- Exhaustive over every graph containing the cycle 0..n-1 at n = 4-7: 644 UH graphs, no consecutive non-critical
  pair.
- Random UH graphs: 924 (Python) and 1,800 (C) at n = 6-11, none.
- Two non-critical vertices can be adjacent through a chord (first at n = 9). So the non-critical set is independent
  along C but not in G.
- In the square-sum UH graphs Q_32 and Q_33 (THM-4851), the non-critical sets have sizes 6 and 9, are pairwise
  non-adjacent in the graph, and contain no two cycle-consecutive vertices.

**Pre-registered predictions.**
- **(B1)** `uh_graph_noncritical_20261104 8 8 0 1 0.5` checks exhaustively every chord set on the 8-cycle (2^20 sets).
  It finds no consecutive non-critical pair.
- **(B2)** `uh_graph_noncritical_20261104 12 15 750 20261104 0.12 0.18 0.24 0.30` draws 3,000 random UH graphs,
  750 at each n = 12-15. It finds no consecutive non-critical pair. The rejection loop is capped at 2 * 10^7 trials
  and reports its yield.

Every outcome will be reported, with the counterexample if there is one.

## Results

(Filled after the confirmatory runs.)
