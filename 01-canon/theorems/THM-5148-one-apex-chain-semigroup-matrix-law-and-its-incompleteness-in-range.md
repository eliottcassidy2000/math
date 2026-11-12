---
id: THM-5148
title: "The one-apex chain semigroup: a matrix law for joins with an inserted transitive chain, and its incompleteness in range. (A) For rooted tournaments (X, u), (Y, v) and the ordinal sum T = X => C => Y with a transitive chain C of i >= 1 vertices in the middle and the arc between u and v reversed: H(T) = H(X) H(Y) + 2^i W+_X(u) W-_Y(v) and W+_T(v) = H(X) W+_Y(v) + 2^i W+_X(u) tau(Y - v), where W+ (W-) is the sum of H(X - P) over the directed paths P starting (ending) at the root and tau the same sum over all paths. (B) For one-apex blocks with apex words of values p_0, ..., p_k the count is the top-left entry of M(p_k, i_k)...M(p_1, i_1) M(p_0, 1), M(p, i) = [[g(p), 2^i (R(p)+1)], [p+1, 2^(i+L(p))]]. (C) This semigroup is not complete in range: its full closure misses 1214 odd numbers below 2^14, and products of at most three matrices miss 5.26 per cent of the odd numbers of [2^17, 2^18), with mean multiplicity 82.9."
status: >
  (A) PROVED from the KNOWN increment formula of the Hamiltonian path count (THM-4876; odd-cycle formula). The join
  law without an inserted chain is KNOWN (THM-4866 (B), Lorentzian form in rooted triples; Session263, Theorem RJ2);
  (A) is its case 'one side made neutral by a chain', with the root moved to the new apex. Not found in one sweep
  of the repo (sub-agent digest, 2026-11-11) in this form; no novelty is claimed beyond the formulation.
  (B) PROVED from (A) and the KNOWN one-apex data (THM-4875, THM-4885: g, rooted triple (g(v), e(v), 1 + v));
  brute force 1470 of 1470 chains of depth 2 and 3 (dynamic programming over subsets).
  (C) FINITE-EXACT, pre-registered (commit 7479352f36): predictions C1, C2, C5 CONFIRMED; C3 (fewer than 0.5 per
  cent missed at depth 3) and C4 (full closure below 2^14 misses only 7 and 21) FAILED.
  AUDIT: pending at this commit (PRE-AUDIT DRAFT).
  The Chan-Pak conjecture remains OPEN; this theorem is a negative result about one candidate family. Completeness
  in range is carried by the joins without an inserted chain and by the bank of small rooted tournaments
  (THM-4885 (D): every odd number up to 2^35 except 7 and 21).
session: kind-pasteur-2026-11-11-S25
source: 05-knowledge/results/kps_s25_euler_field_bridge_law_20261111/README.md
scripts:
  - 05-knowledge/results/kps_s25_euler_field_bridge_law_20261111/scripts/chain_semigroup.py
related:
  - 01-canon/theorems/THM-4866-single-cycle-floor-bottom-band-uniquely-hamiltonian-one-arc-lorentzian-joins-and-crevasses.md
  - 01-canon/theorems/THM-4885-one-apex-join-representation-mass-pi2-over-12.md
  - 01-canon/theorems/THM-5141-one-apex-counts-three-quarters-exponent-reversal-graph-energy.md
---

Notation. `H` = number of Hamiltonian paths; for a tournament `X` and a vertex `u`:
`W+_X(u) = sum over directed paths P starting at u of H(X - P)`, `W-_X(u)` the same for paths ending at `u`,
`tau(Z) = sum over all directed paths P of Z, the empty one included, of H(Z - P)`; `H` of the empty tournament is
1. `g(0) = 1`, `g(2v+1) = g(v)`, `g(2v) = g(v) + v + (v mod 2)`; `L(p)` = bit length of `p`; `R(p)` = reversal of
the complement of `p` on `L(p)` bits. Proofs and data: part B of the source note.

## (A) Join with an inserted chain

`T = X => C => Y`, `C` transitive with `i >= 1` vertices, all arcs from earlier to later blocks except `v -> u`:

    H(T) = H(X) H(Y) + 2^i W+_X(u) W-_Y(v),        W+_T(v) = H(X) W+_Y(v) + 2^i W+_X(u) tau(Y - v).

So the pair `(H, W+)` at the moving root is transported by the matrix
`[[H(Y), 2^i W-_Y(v)], [W+_Y(v), 2^i tau(Y - v)]]`.

## (B) One-apex blocks

A one-apex tournament rooted at its apex, apex word of value `p` (letter 1 = the apex beats the spine vertex, first
spine vertex = top bit): `H = g(p)`, `W+ = p + 1`, `W- = R(p) + 1`, `tau = 2^(L(p))`. Hence for a chain of one-apex
blocks `Y_0, ..., Y_k` with chains of `i_1, ..., i_k >= 1` vertices between them and the arcs between consecutive
apexes reversed,

    H = top-left entry of M(p_k, i_k) ... M(p_1, i_1) M(p_0, 1),     M(p, i) = [[ g(p), 2^i (R(p)+1) ], [ p+1, 2^(i+L(p)) ]].

## (C) Coverage (exact, pre-registered)

- Depth at most 2: 11.06 per cent of the odd numbers of `[2^17, 2^18)` are missed (7250 of 65536); 16,430 below
  `2^18`.
- Depth at most 3: 5.26 per cent (3446 of 65536); 9190 below `2^18`; 10,866,860 products, mean multiplicity 82.9
  per odd number.
- Full closure (all depths): 414 odd numbers missed below `2^12`, 1214 below `2^14` (28,155 live states). First
  misses: 7, 13, 21, 29, 35, 41, 43, 49, 55, 59, 67, 71, 77, 83, 89, 101.

A pure join gives `N = A_1 A_2 + 2^i W_1 W_2` with `W >= A` on both sides, so `N >= 3 A_1 A_2`: the values are
numerous and unevenly distributed.
