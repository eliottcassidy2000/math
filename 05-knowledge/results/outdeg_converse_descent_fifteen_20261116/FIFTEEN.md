# Fifteen and twenty-five: moves on five labels, supersingular primes, the first 25 square-sum nodes — typed

Owner: *"Think the 15 distinct moves on 5 labels and the 15 super singular primes and monster and baby monster and first
25 nodes of the square-sum graph cut up in chunks of 5 and arranged multiple meaningful ways you discover. Think the
fermat and ogg primes. Also extensively immerse yourself in the many prior sessions on the first 25 primes, and draw
analogy to the 25 first square-sum graphs."*

## 0. What the repo already settled (inherited verdicts, unchanged)

- THM-4860 (`machines_automata_fractran_tournaments_dictionary_20261101`): the prime square-sum graph (primes `p ~ q`
  iff `p + q` is a square) is bipartite by `p mod 4`; the 25 primes below 100 by Ogg layer (Fi22 and M∖B perfect
  matchings, pentads edgeless) are **NUMEROLOGY**; the first 50,396 primes admit a consecutive-square-sum arrangement and
  50,396 is the least such count.
- THM-4861 (`fermat_monster_hamiltonian_counts_prime_classes_20261102`): Fermat primes as Hamiltonian path counts with
  least orders 3, 4, 6, 8, 11; Monster layers vs residue classes **COINCIDENCE**.
- THM-4505 (square-sum degree law at 15, endpoints 8, 9), THM-4840 (`7/15/5 = PG(2,2)`, the duads of six, `1 + 4`;
  Hamming(15) = equal-parity-degree subgraphs of `K_6`), LEM-003 / THM-2477 (`H(RT_5) = 15`, two Hamiltonian cycles).
- Memory and `25-prime groups vs square-sum`: NUMEROLOGY, "don't redo". Nothing below re-opens these; the new facts are
  small computations, each typed.

## 1. Fifteen moves on five labels ([receipts/fifteen.out](receipts/fifteen.out) (1))

The 15 moves are the double transpositions `(ab)(cd)` of `S_5` = the 15 involutions of `A_5 ≅ PSL(2,5)` = the 15 edges of
the Petersen graph (disjoint duads) = the 15 two-fold axes of the icosahedron. The regular tournament `RT_5` (circulant
`{1, 2}`) has exactly 15 Hamiltonian paths and `|Aut(RT_5)| = 5`. Under this `Z_5` the paths fall into **three orbits of
five**, labelled by the step patterns `(1,1,1,1)`, `(2,2,2,2)`, `(2,1,1,2)`; the double transpositions fall into three
orbits of five under the same `Z_5` (conjugation by the 5-cycle). Any orbit-matching bijection is `Z_5`-equivariant;
nothing finer exists (the full `S_5` acts on the involutions, only `Z_5` on the paths). P-7 confirmed. **Type: ANALOGY**
(a counting coincidence sharing a `Z_5`).

## 2. The supersingular fifteen, the Monster, the Baby Monster, Fermat ([receipts/fifteen.out](receipts/fifteen.out) (3))

Ogg's 15 supersingular primes `2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 41, 47, 59, 71` (the primes dividing `|M|`) all
lie below 100, so the 25 primes below 100 split 15 + 10, the ten being `37, 43, 53, 61, 67, 73, 79, 83, 89, 97`; the Baby
Monster's primes are the fifteen minus `29, 41, 59, 71`; the Fermat primes `3, 5, 17` are supersingular, `257` and
`65537` are not. In THM-4860's prime square-sum graph on the 25 primes (24 edges, degrees ≤ 4, no Hamiltonian path —
P-8 confirmed — three isolated vertices `37, 43, 67`): the fifteen supersingular primes induce **a tree** (14 edges,
connected: `2–7, 2–23, 2–47, 3–13, 5–11, 5–31, 5–59, 7–29, 13–23, 17–19, 17–47, 23–41, 29–71, 41–59`), the ten others
induce the single edge `61 + 83 = 144` (P-14) and contain all three isolated vertices, and 9 edges cross. The
supersingular primes are the small primes, where the square-sum graph is densest, and `67`'s isolation is `67 + q` never a
square for `q < 100` (`67 + 54`, `67 + 77` are not prime sums); no mechanism relates genus-zero `X_0(p)⁺` to square sums.
**Type: NUMEROLOGY**, consistent with THM-4860 and THM-4861.

## 3. The first 25 square-sum nodes in chunks of five ([receipts/squaresum25.out](receipts/squaresum25.out), [receipts/fifteen.out](receipts/fifteen.out) (2), [receipts/hampath_count.out](receipts/hampath_count.out))

The square-sum graph on `1..25` (`i ~ j` iff `i + j` is a square) has 32 edges and degrees
`4,3,4,3,3,3,3,2,2,2,3,3,3,3,3,2,2,1,2,2,2,2,2,3,2`.

- **Antidiagonal matchings.** The edges are the six matchings `i + j = s`: sizes `1, 4, 7, 12, 7, 1` for
  `s = 4, 9, 16, 25, 36, 49`; the middle one, `s = 25 = 5²`, is a perfect matching of `1..24`.
- **Rows of five** (`1–5 | 6–10 | 11–15 | 16–20 | 21–25`): every row has an internal edge (`1+3, 7+9, 11+14, 16+20,
  24+25`; P-5 refuted), the row-pair edge matrix is `[[4,3,5,1,4],[3,4,1,4,0],[5,1,4,0,5],[1,4,0,4,0],[4,0,5,0,2]]`, and
  the row pairs with no edge are `(2,5), (3,4), (4,5)`.
- **Residue classes mod 5.** Internal edges only in classes `0` (`5+20, 10+15`), `2` (`2+7`), `3` (`3+13, 13+23`); the
  class-pair matrix is `[[4,6,0,0,2],[6,0,0,3,5],[0,0,2,5,6],[0,3,5,4,0],[2,5,6,0,0]]`.
- **The unique Hamiltonian path of `1..15`**, `8,1,15,10,6,3,13,12,4,5,11,14,2,7,9`, cut into chunks of five:
  `(8,1,15,10,6) (3,13,12,4,5) (11,14,2,7,9)`, chunk sums `40, 37, 43`.
- **Hamiltonian path counts** (`hampath_count.c`, DFS from the forced degree-1 endpoint): `n = 15, 16, 17`: one path each
  (endpoints `8–9`, `8–16`, `16–17`); `n = 18, …, 22`: none; `n = 23`: 3; `n = 24`: none; `n = 25`: **10**, all starting
  at `18` (degree 1) and ending at `2, 3, 4, 8, 9, 10, 11, 13, 22, 23`. P-9 confirmed.

## 4. The analogy "first 25 primes ↔ first 25 square-sum nodes"

Same construction (sum a square), different arithmetic: the integer graph on `1..25` is connected, not bipartite (the
5-cycle `1, 3, 6, 10, 15`), 32 edges, 10 Hamiltonian paths; the prime graph on the 25 primes below 100 is bipartite
(`p mod 4`), has 24 edges, four components (three isolated), no Hamiltonian path — and the first prime count with a
Hamiltonian path is 50,396 (THM-4860) against 15 for the integers. Chunks of five align nothing: the primes' rows of five
have internal edges only in rows 1–2, the integers' in every row. **Type: ANALOGY** (the shared construction) /
**NUMEROLOGY** (any pairing of the two 25-sets). The repo's prior sessions on the first 25 primes (THM-4860 §Ogg, THM-4861
§Monster) reached the same verdict; the only proved structure on the prime side is the mod-4 bipartition and the 50,396
threshold, on the integer side the degree law at 15 (THM-4505).
