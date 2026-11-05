---
id: THM-4879
title: "Two-port transitive seam gadgets and the seam parity. In SEAM's three-coordinate transfer v = (H, S, R), the transitive gadget TT_k (vertex 0 the source) with entry u and exit z has a closed integer matrix. Its H-row is (1, -[u = k-1], 2^(k-1-u)). The R-row is (2^z, -mu, c), where mu = 2^u if z = k-1 > u, mu = 2^z if u = k-1 > z, and mu = 0 otherwise. The R-to-R coefficient is c = 2^(k-1-u+z) for z <= u and c = 2^(k-z+u) 3^(z-u-1) for z > u, because each vertex strictly between an earlier entry and a later exit has three admissible slots. Separately, every rooted tournament has R ≡ S (mod 2). So on the plane S = 0 every state has H odd and R even. Every gadget whose S-row vanishes has an even entry R_z(G), and its 2x2 restriction has determinant ≡ M_(z,u)(G) (mod 2)"
status: >
  PROVED (kind-pasteur-2026-11-04-S11; elementary).
  (A) restates the general row formulas of the seam transfer (KNOWN in substance: SEAM.md eqs. (1)-(3) build the
  matrix; opus-S78 THM-4875 (B) gives closed entries in the converse convention). (B) and (C) are proved here by
  direct counting. (D) is a one-line consequence of Redei's theorem and the singleton seam.
  CHECKED: (B) against SEAM's own matrix() for all 204 cases k <= 8; (D) on all 202,013 labelled rooted tournaments
  of order <= 6 and 300 random ones of orders 7-10; (E) on all 1,995 S-row-zero gadget matrices of order <= 7 (per-order
  files), where M_(z,u)(G) and R_z(G) were also always even.
  Audit pending.
session: kind-pasteur-2026-11-04-S11
date: 2026-11-04
related: [THM-4875, THM-4865, HYP-9471, THM-4876]
script: 04-computation/experiments/seam_gadget_density_20261104/ (transitive_two_port.py, parity_check.py)
---

# THM-4879 - two-port transitive seam gadgets and the seam parity

## Setting (SEAM convention)

For a tournament X with a marked exit a:
- **The state** v(X, a) = (H, S, R), where H = H(X), S = S_a(X) is the number of Hamiltonian paths starting at a, and
  R = R_a(X) = sum over U ∋ a of E(X[U], a) H(X \ U).
- **Conventions.** E counts paths ending at a, H(empty) = 1, and singleton paths are allowed.
- **The seam.** For a marked gadget (G, u, z) (u the entry, z the exit, u = z allowed), X ⊙ G keeps all internal arcs
  and orients every cross arc from G to X, except the seam arc a -> u. Its exit is z.
- **The transfer** (SEAM.md (1)) is v(X ⊙ G, z) = M_G v(X, a).

## Statement

**(A) Row formulas** (KNOWN in substance; restated in this convention). For every marked gadget (G, u, z):

    H' = h(G) H  -  E_u(G) S  +  L_u(G) R,
    S' = s_z(G) H  -  P_(z->u)(G) S  +  L_(z,u)(G) R,
    R' = R_z(G) H  -  M_(z,u)(G) S  +  (K1 + K2)(G) R,

where all counts are of ordered vertex-disjoint path tuples covering V(G):
- **Path counts on G.** h, E_u and s_z count Hamiltonian paths of G (all of them, those ending at u, and those
  starting at z). P_(z->u) counts Hamiltonian paths from z to u.
- **L_u** counts pairs (A1, A2): A1 any path, possibly empty, and A2 nonempty, starting at u.
- **L_(z,u)** counts pairs (A1, A2), both nonempty, with A1 starting at z and A2 starting at u.
- **R_z(G)** is the R-statistic of (G, z).
- **M_(z,u)** counts pairs (P, Q), both nonempty, with P ending at z and Q ending at u.
- **K1** counts triples (A1, A2, Q): A2 a path from u to z, with A1 and Q any paths, possibly empty.
- **K2** counts triples (P, Q1, Q3): P ending at z, Q1 any path, possibly empty, and Q3 nonempty, starting at u.

**(B) Two-port transitive gadgets.** Let G = TT_k on 0, ..., k-1 with i -> j iff i < j. Then:

    H-row:  (1, -[u = k-1], 2^(k-1-u)),
    S-row:  k = 1: (1, -1, 0);
            z = 0, k >= 2: (1, 0, 0) if u = 0, and the H-row if u >= 1;
            z >= 1: (0, 0, 2^(k-1-z)) if u = 0, and (0, 0, 0) if u >= 1,
    R-row:  (2^z, -mu, c), where
            mu = 2^u if z = k-1 > u,   mu = 2^z if u = k-1 > z,   mu = 0 otherwise;
            c = 2^(k-1-u+z) if z <= u,   c = 2^(k-z+u) 3^(z-u-1) if z > u.

The single-port case u = z has only powers of 2. It is the transitive-gadget formula of Session 195
(`session195/rooted-interfaces/EQUAL_RESPONSE_BARRIER.md` section 2, in the four-coordinate convention there). The
factor 3^(z-u-1) appears exactly when the entry precedes the exit.

**(C) Ranks.** Note M_(z,u)(G) = 0 when u = z. Among transitive gadgets:
- the matrix is unimodular iff (u, z) = (0, k-1), with det 1;
- every other two-port transitive gadget with k >= 2 is singular.

(This was checked for k <= 8. The general proof is the determinant of the closed form, (B).)

**(D) The seam parity.** For every rooted tournament (X, a), R_a(X) ≡ S_a(X) (mod 2).

**(E) Consequences on the plane P = {S = 0}.**
- States in P have H odd and R even.
- If the S-row of M_G vanishes (G maps every state into P), then R_z(G) is even.
- The determinant of the restriction [[h(G), L_u(G)], [R_z(G), K1 + K2]] acting on (H, R) is ≡ M_(z,u)(G)
  (mod 2).

## Proof

**(A).**
- **Shape of a path.** In X ⊙ G every cross arc runs from G to X except a -> u. So a path moves from X to G only
  through a -> u, at most once. A Hamiltonian path of X ⊙ G is therefore either:
  - (i) a G-path followed by an X-path; or
  - (ii) [A1][B1][A2][B2], with A1 a G-path (possibly empty), B1 an X-path ending at a, A2 a G-path starting at u,
    and B2 an X-path (possibly empty).
- **H'.** In (i), every G -> X junction is an arc except last = u, first = a. That gives h(G)H - E_u(G)S.
  - In (ii), the junctions A1 -> B1 and A2 -> B2 are arcs, because u lies in A2 and a in B1.
  - The X-part of (ii) contributes sum over B1 of E(X[B1], a) H(X \ B1) = R, and the G-part contributes L_u(G).
- **S'.** This is the same split restricted to paths whose first vertex is z: in (i) the G-path starts at z, and in
  (ii) A1 is nonempty and starts at z.
- **R'.** Count pairs (P, Q) with P ending at z and Q covering the rest. The seam arc is used by neither, by P, or by
  Q.
  - *Neither.* P lies in G (it ends in G and cannot re-enter G). Q is a G-part followed by all of X, with the junction
    rule as in (i). This gives R_z(G)H - M_(z,u)(G)S.
  - *By P.* P = [A1][B1 -> a] u [A2 ending at z], and Q is a G-part followed by an X-part. The X-parts give R, and
    the G-parts give K1.
  - *By Q.* Q = [Q1][Q2 ending at a] u [Q3][Q4], and P lies in G. The X-parts give R, and the G-parts give K2.
  - In each case the forbidden junction (u, a) cannot occur, because u and a are already placed.

**(B).** Every vertex subset of TT_k carries exactly one path, in increasing order. Each count in (A) is therefore
a number of slot assignments.
- **Basic counts.** h = 1, E_u = [u = k-1], s_z = [z = 0], P_(z->u) = [z = 0, u = k-1], and R_z(G) = 2^z (a subset with
  maximum z).
- **L_u = 2^(k-1-u).** A2 has minimum u, and each vertex above u is free.
- **L_(z,u).** Nonzero only if min(z, u) = 0 and z ≠ u. It is 2^(k-1-u) for z = 0 < u and 2^(k-1-z) for u = 0 < z.
- **M_(z,u).** Nonzero only if z ≠ u and max(z, u) = k-1. It is 2^u when z = k-1, and 2^z when u = k-1.
- **K1** needs u <= z.
  - For u < z, A2 = {u, z} ∪ (any subset of the vertices strictly between), and every other vertex chooses A1 or Q.
    So K1 = sum over S of 2^(k-2-|S|) = 2^(k-1-z+u) 3^(z-u-1).
  - For u = z, K1 = 2^(k-1).
- **K2** needs u ≠ z. A vertex x outside {u, z} may join P if x < z, join Q3 if x > u, and join Q1 always.
  - For z < u the slot counts multiply to 2^z * 1 * 2^(k-1-u).
  - For z > u they multiply to 2^u 3^(z-u-1) 2^(k-1-z), since every vertex strictly between has all three slots.
- **c = K1 + K2.** Summing gives the stated c in all three cases (u = z, z < u, z > u).

**(D).** The singleton seam gives H(X ⊙ {z}) = H - S + R (the first row of SEAM's singleton matrix [1 -1 1; 1 -1 0;
1 0 1], and klein's T_1 = h - e + B in THM-4865). X ⊙ {z} and X are tournaments, so both H-values are odd (Redei).
Hence R - S is even.

**(E).**
- **P-states.** If S = 0 then R is even by (D), and H is odd by Redei.
- **R_z(G) even.** If the S-row of M_G vanishes, then on input (1, 1, 1) the output's S-coordinate is 0. Since S' =
  s_z(G)H - ..., this forces s_z(G) = 0, and (D) for (G, z) gives R_z(G) even.
- **Parity of the R-row.** The output always lies in P, so R' ≡ S' = 0 (mod 2) on every realisable input. The inputs
  (1, 1, 1) and (1, 0, 2) are realisable (singleton; TT_2 at its sink) and give R_z(G) ≡ 0 and
  R_z(G) - M + (K1 + K2) ≡ 0. So K1 + K2 ≡ M (mod 2).
- **Determinant.** It equals h(G)(K1 + K2) - L_u R_z(G) ≡ K1 + K2 ≡ M_(z,u)(G) (mod 2), since h(G) is odd and
  R_z(G) is even. QED.

## Checks

- `transitive_two_port.py` compares the closed form (B) with SEAM's own `matrix()` (an independent construction: it
  contracts the prefix) for every k <= 8 and every port pair: 204 cases, 0 mismatches.
- `parity_check.py` checks (D) on all 202,013 labelled rooted tournaments of order <= 6 and 300 random rooted
  tournaments of orders 7-10, with 0 violations.
- (E) was checked on the 1,995 S-row-zero matrices in the per-order files of order <= 7 (`gadget_matrices.c`). In all
  of them det ≡ M_(z,u) (mod 2), and also M_(z,u)(G) and R_z(G) are even, so every such 2x2 restriction of order <= 7
  has even determinant.
  - Whether M_(z,u)(G) is always even for S-row-zero gadgets is OPEN. It is the order-8 test P5 of HYP-9471, which
    was blind when this file was written.

## What this does and does not give

- **Explicit family.** (B) is an explicit all-order family of composable pieces, and the only one here with a
  3-power entry.
- **Not a spectrum route.** Seam words built only from transitive gadgets are thin. Over k <= 10 and all port pairs
  they attain only about 3/4 of the odd numbers <= 2,000: 261 are missed, among them 11, 19, 35 and 37
  (FINITE-EXACT, `ttbfs` in the session's note). They cannot carry the {7, 21} spectrum alone.
- **Parity.** (D) and (E) are parity facts. They forbid odd-determinant maps on P only for gadgets with
  M_(z,u)(G) odd, and no such gadget has been seen.

## Credits

- SEAM.md (2026-10-29) for the transfer, its spanning-prefix construction and the singleton matrix.
- opus-S78 THM-4875 (B) for closed entries in the converse convention.
- Session 195 for the single-port transitive formulas.
- klein THM-4865 for the tower step T_1.
