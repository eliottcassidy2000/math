# Paley doubles and the odd spectrum of even graphs

Session mac-mini-2026-11-09-paleyeven, 2026-11-09. Canon: THM-4995, THM-4996, HYP-9575. Proofs: [PROOF.md](PROOF.md).
Audit: [AUDIT_P1.md](AUDIT_P1.md), PASS WITH CORRECTIONS (1 MAJOR, 13 MINOR; all applied; MISTAKE-746), with an
addendum for the order-16 example. Code and outputs: `04-computation/experiments/paleyeven_20261109/`.
Pre-registration (Q1-Q4, committed before the runs): `PREREGISTRATION.md` there.

The owner's directive had two parts:
- build escapers of large module defect from Paley halves (HYP-9565 (A));
- merge in the observation that tournaments and automorphism-even graphs are equinumerous. The question is how an odd
  spectrum missing 7 and 21 shows up on the even-graph side, read through Chan-Pak (Conjecture 6.2, "h takes all odd
  values except 7 and 21", and their notions of almost complete and concise counting functions). Use arXiv 2610.31289
  and 2610.28149 as inspiration.

## 1. Paley halves: doubles of doubly regular tournaments (THM-4995)

**Priority first.** The "module defect" md, which is THM-4985's "decomposability defect", is the classical
**simplicity index** s(T): the fewest arc reversals that create a nontrivial module.
- Müller and Pelant (1974) proved s(T) <= (n - 1)/2, with equality exactly for doubly regular tournaments.
- See also Boussaïri-Lakhlifi-Talbaoui (2022). Audit P1 caught this missing priority (MISTAKE-746).

The **double** D(B) of a tournament B is defined as follows.
- Two copies A and B' of B.
- Every arc goes from A to B', except a perfect matching of back arcs b_x -> a_x.
- For B = QR_p this is the Paley double T_p, of order 2p.

**Theorem (PROVED; new).** For every doubly regular B of order q >= 7, md(D(B)) = lambda(D(B)) = delta0(D(B)) =
(q + 1)/2. D(B) is strong and prime, Aut(B) acts on it diagonally, and D(QR_p) is self-converse.
- The proof is a thirteen-case analysis over how a vertex set meets the two halves.
- Its input is Müller-Pelant's md(B) = (q - 1)/2, which is re-proved here from the Paley/Hadamard identity
  S S^T = qI - J by a Parseval and Cauchy-Schwarz bound.

So the Paley halves supply the family HYP-9565 (A) asked for: simplicity index (q + 1)/2, unbounded.

| object | order | md = lambda | H | H(T - v) | status |
|---|---|---|---|---|---|
| D(B8), B8 near-regular | 16 | 4 | 144,227,835 | 8 values, 16,750,389-21,469,761 | **certified escaper** (all 9 values witnessed) |
| T_7 = D(QR_7) | 14 | 4 | 4,730,397 | 686,457 | not certifiable: above the best known counts at orders 13 and 12; undecided |
| D(C_9{1,2,3,5}) | 18 | 5 | 6,450,103,269 | 773,858,997 | **certified escaper** (1.19 E[H(17)], 1.21 E[H(16)]) |
| T_11 = D(QR_11) | 22 | 6 | 17,667,825,987,845 | 1,822,800,617,515 | escaper under (S_0.4): 0.363 E[H(21)], 0.393 E[H(20)] |
| chain QR_7 x 3 | 21 | 4 | 110,624,292,765 | 9.48e9 / 1.64e10 | far below E[H(20)] |

Further notes on the table:
- **Unconditional records.** The previous record simplicity index of an escaper was 3 (THM-4985, order 12). This
  session certifies index 4 at order 16 and index 5 at order 18.
- **Brégman.** Brégman's permanent bound, applied to Hamiltonian paths as in Alon (1990), shows H(D(B)) and every
  H(D(B) - v) lie below the means E[H(2q - 1)] and E[H(2q - 2)] for every q >= 55. The explicit bound is
  R(q) <= q log(3/4) + 2.85 log q + 4.8.
- **Conditional escapers.** Suppose the odd values up to c m!/2^(m - 1) are all attained at order m (hypothesis
  (S_c); the exact censuses show c = 2.06 at m = 11). Then D(B) is an escaper for all large q, so HYP-9565 (A) holds
  conditionally.
- **Negative findings.**
  - The back arcs multiply H heavily: H(T_7)/H(QR_7)^2 = 132 and H(T_11)/H(QR_11)^2 = 1954. So the heuristic
    "H about H(A)H(B) m^{O(r)}" of HYP-9565 was far too optimistic at small orders.
  - Denser back-arc patterns lower md: 2-regular patterns give md 2, and the Fano incidence gives md 0.
  - Order-15 attempts reach md 4, but each has one deletion above the best known order-13 count.

## 2. Even graphs and an odd spectrum: dominating sets (THM-4996)

Tournaments and even graphs are equinumerous.
- An even graph is one where no automorphism reverses an odd number of edges. The notion is von Brömssen's.
- Royle-Praeger-Glasby-Freedman-Devillers proved the equinumerosity in 2023 (THM-4531).
- There is no natural bijection (THM-4722), so H cannot be carried over.
- Session232 studied automorphism orbits of Hamiltonian paths (THM-4987-4990).

This session looks instead for an **intrinsic odd count on graphs**, the graph analogue of Rédei's theorem.
Brouwer, Csorba and Schrijver (2009) proved that every graph has an odd number d(G) of dominating sets. Like H, d is
multiplicative, over disjoint unions instead of ordinal sums.

**Results.**
- **Completeness over all graphs (PROVED, elementary).** d takes every odd value on all graphs, concisely: on
  ceil(log2(k + 1)) vertices, via threshold graphs. The cone step is the x = 1 case of Alikhani-Peng's join formula.
- **Exact census below 512 (FINITE-EXACT; reproduced independently by audit P1).** An exhaustive census covers every
  connected graph with d <= 511: 474,101 classes, all connected graphs on <= 9 vertices and the relevant ones on 10
  and 11.
  - On even graphs, d takes **every odd value below 512 except 3, 7 and 13**.
  - The only graphs with d in {3, 7, 13} are K_2, K_3 and K_4 - e (with isolated vertices). Each has adjacent twins,
    whose transposition reverses exactly one edge.
  - 21 = d(K_1 v (P_3 u K_1)) is even.
  - The complementary class of odd graphs misses exactly {1, 5, 37}.
- **Parity code validation.** The Euler transform of the census's connected-even counts reproduces A000568 through
  order 9: a count-by-count re-check of the equinumerosity, which von Brömssen had verified to n = 10.
- **Constructive evidence.** Unions, cones and joins with proved parity reach all but about 0.6% of the odd values up
  to 2^22. The rest cluster near 2^m - 2^k.
- **The OCF statistic, for comparison.** I(G, 2) restricted to even graphs misses 5, 7, 13, 15, 21, 23 and 31 below
  32. It contains 7 and 21 but many more, so it is not the right analogue.

**Reading, in Chan-Pak's terms (ANALOGY + CONJECTURE, HYP-9575 (B)).**

| tournaments | even graphs |
|---|---|
| H(T) odd (Rédei) | d(G) odd (Brouwer-Csorba-Schrijver) |
| H multiplicative over ordinal sums; strong components | d multiplicative over disjoint unions; components |
| Chan-Pak Conj. 6.2: exceptions {7, 21} | HYP-9575 (B): exceptions {3, 7, 13} (exact below 512) |
| 7 and 21 excluded by the floors of strong components (THM-1370) | 3, 7, 13 excluded because their only realisers have adjacent twins |
| equinumerous with even graphs (no natural bijection) | equinumerous with tournaments; the odd graphs miss {1, 5, 37} |

The number 7 is excluded on both sides. The mechanisms differ: on the graph side the exclusion is a pure
automorphism-parity statement. 21 survives on the graph side because it has a realiser whose only twins are
non-adjacent.

## 3. Inspiration from the two papers (ANALOGY only)

Both papers are in p-adic geometry. Session232 already recorded source notes (`session232_automorphisms_20261109/`).
Two ideas were useful here.
- **Zhang (arXiv 2610.28149): a genuine chain of r Nygaard charts glued by Frobenius, contracting to the naive
  stack.** Here:
  - the *naive* object is the ordinal sum A (+) B, where H is multiplicative (the Künneth-type statement, which holds
    only for this naive gluing);
  - the *genuine* object glues the charts along a matched set of back arcs;
  - the extra Hamiltonian paths through the gluing play the role of the exceptional fibre: the factors 132 and 1954
    above;
  - the Paley chain QR_7 x 3 is the r = 3 chain.
- **Shao (arXiv 2610.31289): an exhaustion by quasi-compact opens has an almost constant cohomology tower.** Here the
  even-graph spectrum is exhausted by vertex count:
  - the exceptional set below 2^k stabilises at {3, 7, 13} from k = 4 on (exact to 2^9);
  - the hole tower of H stabilises at {7, 21} (klein's censuses).

  Neither is a theorem transfer.

## 4. Certificates, outcomes, process notes

- **Order 18.**
  - 773,858,997 in S(16): `wit_d9_del_m16.out`.
  - 6,450,103,269 in S(17): `wit_d9_H_m17.out`.
- **Order 16.**
  - 144,227,835 in S(15): `wit_d8_H_m15.out`.
  - The eight deletion values in S(14): `wit_d8_del_m14_part1-3.out`.
- **Recounts.** Every witness is re-counted by wit.c's DP, by tprop.c and by audit P1's two algorithms.
- **Pre-registered outcomes** (`PREREGISTRATION.md`):
  - **Q1**: md = lambda = delta0 = 4 CONFIRMED. Escaper status UNDECIDED: not certifiable; it fails if THM-329's
    a(12) is exact. This is corrected from "REFUTED" after audit P1.
  - **Q2** CONFIRMED for md = lambda = 6. In-session certification of T_11 was not achieved.
  - **Q3** CONFIRMED and PROVED, for all doubly regular B with q >= 7. Its input lemma is Müller-Pelant's.
  - **Q4** REFUTED: the exceptional set is {3, 7, 13}, not {3, 7, 21}.
- **Process.**
  - A first order-17 witness run (seed 1709) ended, probably with a witness, but its output was lost: the file had
    been committed, and a `git pull --rebase` replaced it with a new inode. The rerun (seed 4242) found a witness.
  - One order-14 search initially received its eight targets as a single zsh word and searched only the first. It was
    rerun with explicit arguments.
  - Both are recorded as gotchas.
