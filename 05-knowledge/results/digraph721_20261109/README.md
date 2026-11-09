# digraph721: even-graph Chan-Pak by offset moves

Session mac-mini-2026-11-09-digraph721 (2026-11-09).

**Canon.**
- THM-5055: a second, independent proof of THM-5040, the sharp even-graph Chan-Pak theorem (HYP-9575 (B)).
- THM-5056: the odd-graph Chan-Pak theorem, exceptions exactly {1, 5, 37}, each value on ceil(log2(k + 1)) vertices.
- HYP-9625: the odd-graph conjecture with its reduction to primes; PROVED by THM-5056 the same day.

[PROOF.md](PROOF.md). Code: `04-computation/experiments/digraph721_20261109/`.

## Owner directive

The directive (2026-11-09) asked for three things:
- keep going and prove the even-graph Chan-Pak conjecture;
- study Hamiltonian path counts of digraphs as square binary matrices, and the conjecture that totals 7 or 21 force a
  loop (a diagonal 1) or a pair of mirrored zeros;
- relate an n x n matrix to an (n+1)-tournament: the diagonal plus two tilings.

The digraph/matrix part was shared out:
- **Session234** (first on that prompt) refuted the matrix conjecture. A semicomplete 4-vertex matrix with no loop
  and no mirrored zero pair has H = 7, and a 5-vertex one has H = 21. It also built the all-order semicomplete
  classification and core grammar.
- **opus-2026-11-09-S89** took the (n+1)-tiling dictionary (THM-5041) and the ternary side (THM-5042).
- **kind-pasteur-2026-11-09-S18** proved the Grinberg-Stanley weighted cycle formula for digraphs and base-3 port
  chains (THM-5070, THM-5071).

This session's own early check of the tiling dictionary (scratch; tournaments with n <= 5) found no functional relation
between H(D) and H of the two associated (n+1)-tournaments. It took the even-graph lane.

## Results

1. **Even-graph Chan-Pak, sharp (THM-5055; the same statement as THM-5040).** For every n, the dominating-set counts
   of even graphs on n vertices are exactly the odd numbers <= M(n), except 3, 7, 13 and X(n) = 2^n - n - 2 (odd
   n >= 5). Here M(n) = 2^n - n for odd n and 2^n - n - 1 for even n. So every odd number other than 3, 7 and 13 is
   d of an even graph, on ceil(log2(k + 1)) vertices or one more.
2. **The engine: a hub move.** For any graph G, vertex u with |N[u]| = r, and l >= 1, let PH(G, u, l) be G plus a hub
   on every vertex but u, plus l leaves on the hub.
   - d(PH) = d(G) + (2^r - 1) 2^(|G| + l - r).
   - PH is even whenever G is even and the hub has strictly maximum degree.
   - Two cases preserve the offset c = d - 2^(n-1):
     - R (exclude a pendant, one leaf), taking the order from n to n + 2;
     - PH(3, 2) (exclude a degree-2 vertex, two leaves), taking the order from n to n + 3.
   - Starting from a "good" seed (two pendants with non-universal neighbours, and a degree-2 vertex), a single offset
     is therefore realised at every larger order.
3. **Seed lemma.** For every odd j, the value 2^(n-1) - j has a good realisation for all n >= L(j) + 8.
   - For L(j) >= 13 the proof is by cases: a cone over (X u E_2) u K_1, an E_2-join over a padded graph, or reduction
     to a smaller offset by padding.
   - For L(j) <= 12 it is a finite computation.
4. **How the families are generated.** The odd numbers are written by their leading bits:
   - a cone over a padded graph writes 10;
   - an E_2-join writes 11, with a borrow of 2;
   - the hub moves carry the values that neither can reach at any order. These are the values just below 2^(n-1)
     (which need prime graphs) and the shifted exceptions 2^(n-1) + {3, 7, 13} and 3 * 2^(n-2) + {1, 5, 11}.
5. **Why shift 3.** R alone is not enough. 11, 23 and 47 have no even realisation with a pendant at any order (they are
   prime, and none of their connected even realisers has a pendant). Through R this blocks the offset -81 at every
   even order. So the four-move closure is exact to n = 26 but misses 2^(n-1) - 81 for even n >= 84. PH(3, 2) changes
   the parity of the order and removes the block.
6. **Odd graphs (THM-5056; HYP-9625 PROVED).** The d-values of odd graphs are exactly the odd numbers other than 1, 5
   and 37, each attained on ceil(log2(k + 1)) vertices, so oddness costs no vertices.
   - K_2 v G writes every leading 11.
   - Cones over odd graphs write 10.
   - Offsets 2^(L-1) + {1, 5, 37} come from the THM-5055 moves on graphs carrying a persistent pair of adjacent
     twins, which keeps every move odd. The four seeds are K_4 - e, K_2 u C_4, K_2 v (P_3 u E_2) and K_2 u T_55.
   - The base below 512 is the census.
   - Example: 1061 = 2^10 + 37, which is prime and needs a prime graph, is the hub move PH(K_2 u T_55, w, 2) on 11
     vertices.

## Verification

- **Lemma checks.** The lemmas are proved in PROOF.md and checked by brute force: 706 single moves and 217 chains, with
  evenness from enumerating Aut(G).
- **Closure.** The base closure was computed twice (C, and an independent numpy re-implementation). It gives
  E(n) = F(n) for n <= 26, the fixed offsets from order 8, and the seed base for L(j) <= 12.
- **Witnesses.** 6,147 explicit witness graphs (up to 21 vertices) were checked by brute force for d, for evenness
  (via a generating set of Aut(G)) and for flags, with 0 failures.
- **Case split.** The case split of the seed lemma was checked against the computed tables for 13 <= L(j) <= 15.
- **Agreement with THM-5040.** THM-5040's independent closure gives the same E(n) to n = 25.
- **Odd graphs (THM-5056).** Explicit realisers following the proof were checked by brute force for d and for oddness
  (via Aut(G)): every odd k < 8192 (4,093 graphs), 300 random k < 2^20, and 57 offset values up to 20 vertices. All
  pass.

## Pre-registration (PREREGISTRATION.md)

| | prediction | outcome |
|---|---|---|
| R1 (p = 0.50) | explicit even constructions realise every odd k in [15, 2^24] | CONFIRMED, and in fact for all k (the closure to n = 26 covers [1, 2^25]) |
| R2 (p = 0.25) | a complete proof of HYP-9575 (B) this session | CONFIRMED (THM-5055; THM-5040 found one concurrently) |
| R3 (p = 0.70) | partial hubs are essential to R1 | CONFIRMED: without hub moves the closure first fails at 501 = X(9) and misses 2,096,928 values of F(n) for n <= 26 |

## Credits

- von Brömssen (even graphs); RPGFD 2023 (equinumerosity with tournaments).
- Brouwer-Csorba-Schrijver (d odd); Alikhani-Peng (join formula); Alikhani, Graphs Combin. 29 (2013) (completeness
  over all graphs).
- Chan-Pak Conjecture 6.2 (the tournament analogue, {7, 21}, still open).
- THM-4996 (census, parity calculus, the exceptional sets below 512).
- THM-5040 (opus-2026-11-09-S89: the sharp statement, the exclusions, and a parallel proof).
- Session234 and kind-pasteur S18 (digraph lane).
