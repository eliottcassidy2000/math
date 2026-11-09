# Paley doubles and the odd spectrum of even graphs

Session mac-mini-2026-11-09-paleyeven, 2026-11-09. Canon: THM-4995, THM-4996, HYP-9575. Proofs: [PROOF.md](PROOF.md).
Code and outputs: `04-computation/experiments/paleyeven_20261109/`. Pre-registration (Q1-Q4, committed before the
runs): `PREREGISTRATION.md` there.

The owner's directive had two parts:
- build escapers of large module defect from Paley halves (HYP-9565 (A));
- merge in the observation that tournaments and automorphism-even graphs are equinumerous. The question is how an odd
  spectrum missing 7 and 21 shows up on the even-graph side, read through Chan-Pak (Conjecture 6.2, "h takes all odd
  values except 7 and 21", and their notions of almost complete and concise counting functions). Use arXiv 2610.31289
  and 2610.28149 as inspiration.

## 1. Paley halves: doubles of doubly regular tournaments (THM-4995)

The **double** D(B) of a tournament B is defined as follows.
- Two copies A and B' of B.
- Every arc goes from A to B', except a perfect matching of back arcs b_x -> a_x.
- For B = QR_p this is the Paley double T_p, of order 2p.

**Theorem (PROVED).** For every doubly regular B of order q >= 7, md(D(B)) = lambda(D(B)) = delta0(D(B)) = (q + 1)/2.
D(B) is strong and prime, Aut(B) acts on it diagonally, and D(QR_p) is self-converse.
- The proof is a thirteen-case analysis over how a vertex set meets the two halves.
- Its input is the lemma md(B) = (q - 1)/2. That lemma follows from the Paley/Hadamard identity S S^T = qI - J by a
  Parseval and Cauchy-Schwarz bound, for q >= 11 (q = 7 checked).

So the Paley halves supply the family that HYP-9565 (A) asked for: decomposability defect (q + 1)/2, unbounded.

| object | order | md = lambda | H | H(T - v) | position |
|---|---|---|---|---|---|
| T_7 = D(QR_7) | 14 | 4 | 4,730,397 | 686,457 | above THM-329's best order-13 count 3,719,831: a first creator unless max H(13) >= 4,730,397 |
| D(C_9{1,2,3,5}) | 18 | 5 | 6,450,103,269 | 773,858,997 | 1.19 E[H(17)], 1.21 E[H(16)]: both certified by witnesses, so a certified escaper |
| T_11 = D(QR_11) | 22 | 6 | 17,667,825,987,845 | 1,822,800,617,515 | 0.363 E[H(21)], 0.393 E[H(20)]: escaper under (S_0.4) |
| chain QR_7 x 3 | 21 | 4 | 110,624,292,765 | 9.48e9 / 1.64e10 | far below E[H(20)] |

Further notes on the table:
- **Brégman.** Brégman's permanent bound shows H(D(B)) and every H(D(B) - v) lie below the means E[H(2q - 1)] and
  E[H(2q - 2)] for every q >= 55 (log-ratio about q log(3/4)).
- **Conditional escapers.** If the odd values up to c m!/2^(m - 1) are all attained at order m (hypothesis (S_c);
  the exact censuses show c = 2.06 at m = 11), then D(B) is an escaper for all large q. So HYP-9565 (A) holds
  conditionally.
- **Unconditional record.** D(C_9{1,2,3,5}) is a certified escaper of order 18 with defect 5 (previous record 3,
  THM-4985).
- **Negative findings.**
  - The back arcs multiply H heavily: H(T_7)/H(QR_7)^2 = 132 and H(T_11)/H(QR_11)^2 = 1954. So the heuristic
    "H about H(A)H(B) m^{O(r)}" of HYP-9565 was far too optimistic at small orders.
  - Denser back-arc patterns lower the defect: 2-regular patterns give md 2, and the Fano incidence gives md 0
    (`paley_mdefect.out`).

## 2. Even graphs and an odd spectrum: Brouwer's dominating-set count (THM-4996)

Tournaments and Royle-even graphs are equinumerous. A graph is Royle-even if no automorphism reverses an odd number of
edges (RPGFD 2023; THM-4531). There is no natural bijection (THM-4722), so H cannot be carried over. Session232 studied
automorphism orbits of Hamiltonian paths (THM-4987-4990).

This session looks instead for an **intrinsic odd count on graphs**, the graph analogue of Rédei's theorem. Brouwer
(2009) proved that every graph has an odd number d(G) of dominating sets. Like H, d is multiplicative, over disjoint
unions instead of ordinal sums.

**Results.**
- **Completeness over all graphs (PROVED, elementary).** d takes every odd value on all graphs, concisely: on
  ceil(log2(k + 1)) vertices, via threshold graphs.
- **Exact census below 512 (FINITE-EXACT).** An exhaustive census covers every connected graph with d <= 511:
  474,101 classes, all connected graphs on <= 9 vertices and the relevant ones on 10 and 11. It is exhaustive because
  adding a joined vertex strictly increases d.
  - On Royle-even graphs, d takes **every odd value below 512 except 3, 7 and 13**.
  - The only graphs with d in {3, 7, 13} are K_2, K_3 and K_4 - e (with isolated vertices). Each has adjacent twins,
    whose transposition reverses exactly one edge.
  - 21 = d(K_1 v (P_3 u K_1)) is even.
- **Parity code validation.** The census reproduces A001349 and A000088. The Euler transform of its connected-even
  counts is A000568 through order 9 (1, 1, 2, 4, 12, 56, 456, 6,880, 191,536): the equinumerosity theorem, reproduced
  class by class.
- **Constructive evidence.** Unions, cones and joins with proved parity reach all but about 0.6% of the odd values up
  to 2^22. The rest lie next to powers of 2, where cones cannot reach because 3, 7 and 13 are not even values.
- **The OCF statistic, for comparison.** I(G, 2) restricted to even graphs misses 5, 7, 13, 15, 21, 23 and 31 below
  32. It contains 7 and 21 but many more, so it is not the right analogue.

**Reading, in Chan-Pak's terms (ANALOGY + CONJECTURE, HYP-9575 (B)).**

| tournaments | Royle-even graphs |
|---|---|
| H(T) odd (Rédei) | d(G) odd (Brouwer) |
| H multiplicative over ordinal sums; strong components | d multiplicative over disjoint unions; components |
| Chan-Pak Conj. 6.2: exceptions {7, 21} | HYP-9575 (B): exceptions {3, 7, 13} (exact below 512) |
| 7 and 21 excluded by the floors of strong components (THM-1370) | 3, 7, 13 excluded because their only realisers have adjacent twins |
| equinumerous with even graphs (no natural bijection) | equinumerous with tournaments |

The number 7 is excluded on both sides. The mechanisms differ: on the graph side the exclusion is a pure
automorphism-parity statement. 21 survives on the graph side because it has a realiser whose only twins are non-adjacent.

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

## 4. Certificates and status

- **Deletion certificate.** 773,858,997 in S(16): `wit_d9_del_m16.out`, a 16-tournament re-counted by direct DP.
- **H certificate.** 6,450,103,269 in S(17): `wit_d9_H_m17.out`, a 17-tournament re-counted by two programs.
- **Pre-registered outcomes** (`PREREGISTRATION.md`):
  - **Q1** REFUTED: T_7 has md 4, but its H exceeds the best known order-13 count.
  - **Q2** CONFIRMED for md = lambda = 6; certification of T_11 not achieved (conditional only).
  - **Q3** CONFIRMED and PROVED, for all doubly regular B with q >= 7.
  - **Q4** REFUTED: the exceptional set is {3, 7, 13}, not {3, 7, 21}.
