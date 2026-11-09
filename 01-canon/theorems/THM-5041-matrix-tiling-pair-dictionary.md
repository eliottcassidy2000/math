---
id: THM-5041
title: "An n x n binary matrix is a diagonal (n loop bits) plus an ordered pair (U, L) of order-(n+1) tilings. Tournaments are exactly the pairs (t, complement t), i.e. the metagraph lines with an orientation. Semicomplete digraphs are the covering pairs, oriented graphs the disjoint pairs, undirected graphs the equal pairs; transpose swaps U and L. The Hamiltonian path total H ignores the diagonal and is tile-wise additive: H(i->j) + H(j->i) = H(i<->j) + H(i, j missing) at every pair. So the two defect states of a tile carry the same total as the metagraph arc-flip pair through that tile. The owner's conjecture (totals 7 and 21 need a loop or a mirrored zero pair) is restated and corrected in these coordinates"
status: >
  PROVED, elementary. This is a dictionary that consolidates KNOWN pieces; audit2 PASS WITH CORRECTIONS (applied).
  - **KNOWN coordinates.** Session89 (05-knowledge/results/session89_consecutive_tournament_pairs_20261028/encoding/
    ENCODING_AND_ACTIONS.md) already has:
    - section 1: the coordinates x_ij = M_ij, y_ij = M_ji, d_i = M_ii, and symmetric matrices as equal pairs;
    - section 6: the transpose law;
    - section 4: the discrepancy graph, whose zero set is the defect set here.
  - **KNOWN additivity.** Tile-wise additivity follows at once from Session234 DG7: the increment for adding a -> b
    never involves b -> a.
  - **Other KNOWN sources.** THM-830 section 1 (tilings are labelled tournaments; lines are converse pairs). Session234
    DG1 and DG2.
  - **Contributed here:**
    - the tiling reading through THM-830 with its convention fixed (the tiling of A is L(A));
    - the class table;
    - the reversal formula;
    - the explicit metagraph-flip corollary and the signed expansion;
    - the restatement of the owner's conjecture.
session: opus-2026-11-09-S89
date: 2026-11-09
related: [THM-830, THM-4731, THM-5071, THM-5040, HYP-9610]
---

# THM-5041 - Matrices as a diagonal plus a pair of tilings

## 1. Dictionary

Let A be an n x n 0/1 matrix, read as a digraph with possible loops. By THM-830 section 1, a tiling of order n + 1 is
the same as a 0/1 array indexed by the unordered pairs {i < j} of [n] (the path gaps), and it is a labelled tournament
on [n]. Define:
- diag(A) = (A_11, ..., A_nn);
- U(A)_ij = A_ij for i < j;
- L(A)_ij = A_ji for i < j.

**(1a) Bijection.** The map A -> (diag A, U(A), L(A)) is a bijection from n x n matrices to
F_2^n x Tile(n+1) x Tile(n+1). The bit count is n^2 = n + 2 C(n,2). (Session89 splits n^2 differently, as
C(n,2) + C(n+1,2), with the diagonal inside the second tiling as root arcs. Its coordinates x = U, y = L and d = diag
are the ones used here (Session89 encoding, section 1).)

**(1b) Pair states.** At a pair p = {i < j}, the state (U_p, L_p) is one of:
- (1, 0): i -> j;
- (0, 1): j -> i;
- (1, 1): reciprocal;
- (0, 0): missing (the owner's "mirrored zero pair").

**(1c) Classes.** Every statement here ignores loops:

| Digraph class | Condition on (U, L) |
| --- | --- |
| Tournament | L = complement of U: the ordered line (t, t^c) of the metagraph |
| Semicomplete | covering pair, U or L = 1 at every tile |
| Oriented | disjoint pair: U_p L_p = 0 (never both 1) at every tile |
| Undirected (symmetric) | U = L |

Each loopless digraph is a tiling U together with a **defect set** Delta = {p : U_p = L_p}, where the tile type is
U_p. So there are 4^C(n,2) loopless digraphs, 3^C(n,2) semicomplete ones (the ternary triangle) and 2^C(n,2)
tournaments.

**(1d) Symmetries.**
- Transpose: (U, L) -> (L, U).
- Complement J - A: complements all three coordinates.
- Reversal relabelling i -> n + 1 - i: (U, L) -> (s L, s U), where (sX)_ij = X_(n+1-j, n+1-i). This is THM-830's
  staircase reflection sigma, (a, b) -> (n+1-b, n+1-a), transported to gap pairs.
- **Convention.** Under THM-830 (1.1), (G_t)_(i,j) = t_(i+1,j). So the tiling of the tournament A is L(A), and U(A) is
  the tiling of A^op. Audit2 checked (1.3) and the reversal formula over all tilings for n = 4, 5, 6.

On tournaments these give THM-830's dictionary: lines are converse pairs, and blue lines are the
reversal-anti-automorphic tournaments.

## 2. Path totals in these coordinates

**(2a)** H(A) depends neither on diag(A) (Session234 DG1) nor on the order of the pair: H(U, L) = H(L, U), since a
path reversed is a path of the transpose.
- On lines (tournaments), H is odd (Rédei).
- On equal pairs (U, U), H = 2 hp(G_U) for n >= 2, where hp counts undirected Hamiltonian paths; this is even.
- In general, H mod 2 counts packings of defect cycles (THM-5071 (3)). A defect cycle is a cycle of D through a
  (1,1) tile, or a cycle of the loopless complement through a (0,0) tile.

**(2b) Tile-wise additivity.** For every digraph and every pair p, with all other entries fixed,

    H(p = (1,0)) + H(p = (0,1)) = H(p = (1,1)) + H(p = (0,0)).

*Proof.* This also follows at once from Session234 DG7. Directly: H = sum over vertex orderings of the product of A
over consecutive pairs. This is multilinear in the arc
indicators, and no term contains both A_ij and A_ji: a Hamiltonian path cannot use both arcs of one pair. So
H = alpha + beta U_p + gamma L_p, and the identity follows. QED

**(2c) Metagraph-flip corollary.** Let T be a tournament (an ordered line) and p a tile. Let T^(+p) be T with the arc at p
reversed (the single-tile flip). Then

    H(T with p reciprocal) + H(T with p missing) = H(T) + H(T^(+p)).

- Both one-defect totals have the same parity.
- Their difference is the number of Hamiltonian paths of T or T^(+p) that use the pair p.
- Iterating over the reciprocal set R of a semicomplete D gives the signed expansion
  H(D) = sum over A subset of R of (-1)^|A| sum over orientations s of R \ A of H(T_s - A). Here T_s - A is the tournament
  resolution with the pairs of A deleted. Compare Session234 DG6, which sums over resolutions with weights 2^(b - k(P)).

## 3. The owner's 7/21 conjecture in these coordinates

The conjecture was that a total of 7 or 21 forces a 1 on the diagonal or a mirrored zero pair.
- **Loops.** The diagonal never affects H (2a), so a loop cannot be what forces 7 or 21.
- **Missing pairs.** A mirrored zero pair is the state (0, 0) at some tile.
- **The correct statement** (Session234, DG2), in these coordinates:
  - for H = 7, the pair (U, L) has a (0, 0) tile or at least two (1, 1) tiles; on covering pairs, exactly two;
  - for H = 21, it has a (0, 0) tile or at least one (1, 1) tile; on covering pairs, one to three;
  - in particular it is never a line.
- **Refutation.** Session234 refuted the original statement with semicomplete examples, that is, covering pairs with
  only (1, 1) defects.

**The ternary triangle.**
- The covering pairs (each tile in {(1,0), (0,1), (1,1)}) are exactly the owner's "binary hypotenuse, ternary
  triangle" picture once loops are added: 2^n 3^C(n,2) matrices.
- Their H-spectrum is every positive integer (THM-5071 (4a)).
- The binary triangle (lines) misses 7 and 21 among the odd numbers (THM-1370), and conjecturally nothing else
  (Chan-Pak).
- A third pair of this kind is the dominating-set spectrum of graphs: even graphs, which are equinumerous with
  tournaments, miss {3, 7, 13} (THM-5040), while odd graphs miss {1, 5, 37} (THM-5056 first, THM-5042). See HYP-9610 for the typed
  dictionary.
