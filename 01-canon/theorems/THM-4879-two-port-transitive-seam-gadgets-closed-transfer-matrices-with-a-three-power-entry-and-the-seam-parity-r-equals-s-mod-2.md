---
id: THM-4879
title: "Two-port transitive seam gadgets and seam parities. In SEAM's three-coordinate transfer v = (H, S, R), the transitive gadget TT_k (vertex 0 the source) with entry u and exit z has a closed integer matrix. Its H-row is (1, -[u = k-1], 2^(k-1-u)). Its R-row is (2^z, -mu, c), with mu = 2^u if z = k-1 > u, mu = 2^z if u = k-1 > z, and mu = 0 otherwise. The R-to-R coefficient is c = 2^(k-1-u+z) for z <= u and c = 2^(k-z+u) 3^(z-u-1) for z > u, because each vertex strictly between an earlier entry and a later exit has three admissible slots. Every rooted tournament has R ≡ S (mod 2) (KNOWN, Session29). The two-endpoint refinement: for distinct z, u, the counts of ordered two-path covers starting at z and u, and of those ending at z and u, are both ≡ P(z->u) + P(u->z) (mod 2), and every nonempty tournament has an odd number of directed paths. Hence every gadget that maps the plane S = 0 into itself acts there by a 2x2 map of even determinant, at every order"
status: >
  PROVED (kind-pasteur-2026-11-04-S11; elementary).
  (A) restates the general row formulas of the seam transfer (KNOWN in substance: SEAM.md eqs. (1)-(3) build the
  matrix; opus-S78 THM-4875 (B) gives closed entries in the converse convention).
  (B), (C), (E) and (F) are proved here by direct counting and parity.
  (D) is KNOWN: Session29 session29-converse-congruences REPORT section 5 (2026-10-22, same convention), and the mod-2
  reduction of Session 193 multiple_triangles PROOF.md section 2 applied to the converse. The singleton-seam proof
  given here is a shorter alternative.
  (F) is new as far as the audit found. It upgrades the order-8 check of (E) (HYP-9471 P5) to all orders.
  AUDIT: independent adversarial audit PASS WITH CORRECTIONS (4 major, all in the companion note or the credit of
  (D); 12 minor; all applied; MISTAKE-717), 04-computation/experiments/seam_gadget_thickness_20261104_audit_S11/AUDIT.md.
  The audit's own code reproduced:
  (A) by permutation brute force (all labelled gadgets of order <= 5 with every prefix of order <= 3, plus random
  cases);
  (B) and (C) for k <= 9 by three routes;
  (D), (F1) and (F2) on all labelled tournaments of order <= 7 (89 million ordered pairs), plus random ones to
  order 10;
  (E) and (F3) on all 22,522 S-row-zero matrices of order <= 8.
session: kind-pasteur-2026-11-04-S11
date: 2026-11-04
related: [THM-4875, THM-4865, HYP-9471, THM-4876]
script: 04-computation/experiments/seam_gadget_density_20261104/ (transitive_two_port.py, parity_check.py)
---

# THM-4879 - two-port transitive seam gadgets and seam parities

## Setting (SEAM convention)

For a tournament X with a marked exit a:
- **The state** v(X, a) = (H, S, R), where H = H(X), S = S_a(X) is the number of Hamiltonian paths starting at a, and
  R = R_a(X) = sum over U ∋ a of E(X[U], a) H(X \ U).
- **Conventions.** E counts paths ending at a, H(empty) = 1, and singleton paths are allowed.
- **The seam.** For a marked gadget (G, u, z) (u the entry, z the exit, u = z allowed), X ⊙ G keeps all internal arcs
  and orients every cross arc from G to X, except the seam arc a -> u. Its exit is z.
- **The transfer** (SEAM.md (1)) is v(X ⊙ G, z) = M_G v(X, a).
- **Words are gadgets.** Every seam word is itself a marked gadget, so every statement below about gadgets also holds
  for words.

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
(`session195/rooted-interfaces/EQUAL_RESPONSE_BARRIER.md` section 2, in the converse four-coordinate convention
there). The factor 3^(z-u-1) appears exactly when the entry precedes the exit.

**(C) Ranks.** For k >= 2, the matrix of (TT_k, u, z) has determinant 1 if (u, z) = (0, k-1) and is singular
otherwise. (k = 1 is the singleton, det 1.) Note also M_(z,u)(G) = 0 when u = z.

**(D) The seam parity** (KNOWN: Session29, see Status). For every rooted tournament (X, a),
R_a(X) ≡ S_a(X) (mod 2).

**(E) Consequences on the plane P = {S = 0}.**
- States in P have H odd and R even.
- Suppose G maps every state into P, that is, its S-row (s_z, -P_(z->u), L_(z,u)) vanishes. Then R_z(G) is even.
- For such G (S-row zero), the determinant of the restriction [[h(G), L_u(G)], [R_z(G), K1 + K2]] acting on (H, R)
  is ≡ M_(z,u)(G) (mod 2).
  - Without the hypothesis this congruence fails, for example for 171,897 of the distinct order-8 matrices (audit).

**(F) Two-endpoint parity.** Write Pi(W) for the number of nonempty directed paths of a tournament W, of any length,
single vertices included. For a tournament G and distinct vertices z, u, define the Lmix and Mmix analogues of L and
M: Lmix counts covers (A from z, B ending at u), and Mmix counts covers (A ending at z, B from u).

- **(F1)** Pi(W) is odd for every nonempty tournament W.
- **(F2)** L_(z,u)(G) ≡ M_(z,u)(G) ≡ P_(z->u)(G) + P_(u->z)(G) (mod 2), and Lmix ≡ Mmix (mod 2).
- **(F3)** For a gadget (G, u, z), these are equivalent: it maps P into P; s_z(G) = L_(z,u)(G) = 0; its S-row is zero.
  (The last two agree because P_(z->u) <= s_z.) For every such gadget, at every order, M_(z,u)(G) is even, and so
  its 2x2 restriction to P has even determinant. No gadget, and no seam word, acts on P by a unimodular map or by a
  map of odd determinant.

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

**(C).** Read the closed form (B).
- **u, z >= 1.** The S-row is zero.
- **u = 0, 1 <= z <= k-2.** Every row has middle entry 0: the H-row since u ≠ k-1, the S-row (0, 0, 2^(k-1-z)), and
  the R-row since mu = 0. So the middle column is zero.
- **u = z = 0.** The middle column is again zero.
- **z = 0 < u.** The S-row equals the H-row.
- **(u, z) = (0, k-1).** The matrix is [[1, 0, 2^(k-1)], [0, 0, 1], [2^(k-1), -1, c]]. Expanding along the middle
  column gives det = (-(-1)) * det [[1, 2^(k-1)], [0, 1]] = 1.

These cases cover all port pairs with k >= 2. (Checked numerically for k <= 9.)

**(D)** (alternative proof). The singleton seam gives H(X ⊙ {z}) = H - S + R (the first row of SEAM's singleton
matrix [1 -1 1; 1 -1 0; 1 0 1], and klein's T_1 = h - e + B in THM-4865). X ⊙ {z} and X are tournaments, so both
H-values are odd (Redei). Hence R - S is even.
- Session29 derived the same congruence from H(X ⊙ X^op) = H² - S² + R².
- Session 193's identity gives the exact value R - S = 2V in the converse reading.

**(E).**
- **P-states.** If S = 0 then R is even by (D), and H is odd by Redei.
- **R_z(G) even.** By (A) the S-row is (s_z(G), -P_(z->u)(G), L_(z,u)(G)). If it vanishes, then s_z(G) = 0, and (D)
  for (G, z) gives R_z(G) ≡ s_z(G) = 0.
- **Parity of the R-row.** The output always lies in P, so R' ≡ S' = 0 (mod 2) on every realisable input. The
  realisable inputs (1, 1, 1) (singleton) and (1, 0, 2) (TT_2 at its sink) give, respectively,
  R_z(G) - M + (K1 + K2) ≡ 0 and R_z(G) + 2(K1 + K2) ≡ R_z(G) ≡ 0. So K1 + K2 ≡ M (mod 2).
- **Determinant.** It equals h(G)(K1 + K2) - L_u R_z(G) ≡ K1 + K2 ≡ M_(z,u)(G) (mod 2), since h(G) is odd and
  R_z(G) is even.

**(F).**
- **The reduction mod 2.** Write pi_z(T) for the number of paths of T, of any vertex set, ending at z. Every
  H(T \ U) is odd, so R_z(T) ≡ pi_z(T) (mod 2). With (D), pi_z(T) ≡ S_z(T) for every rooted tournament (T, z).
- **(F1).** Pi(W) = sum over v of pi_v(W) ≡ sum over v of S_v(W) = H(W) ≡ 1.
  - Alternatively, Pi(W) is the sum over nonempty U of H(W[U]), a sum of 2^|W| - 1 odd numbers.
- **The auxiliary tournament.** Let T = G + w, where w beats u only and every other vertex of G beats w.
  - A Hamiltonian path of T starting at z either ends with w, preceded by a vertex other than u, or has w in the
    middle, followed by u. So S_z(T) = (s_z(G) - P_(z->u)(G)) + L_(z,u)(G).
  - A path of T ending at z either avoids w, or has the form [A] w [C], where C is a path of G from u to z and A is
    a path of G - V(C), possibly empty. So pi_z(T) = pi_z(G) + sum over u->z paths C of (1 + Pi(G - V(C))).
  - By (F1), 1 + Pi(W) is even unless W is empty, so the sum is ≡ P_(u->z)(G).
- **(F2), the first congruence.** pi_z(T) ≡ S_z(T) and pi_z(G) ≡ s_z(G) give
  L_(z,u)(G) ≡ P_(z->u)(G) + P_(u->z)(G).
- **(F2), M.** Reversing both paths gives M_(z,u)(G) = L_(z,u)(G^op), and the right-hand side is symmetric under
  reversal.
- **(F2), Lmix and Mmix.** By (D) for G^op, e(W, u) ≡ sum over V ⊆ W with V ∋ u of s(V, u). Substituting this into
  Lmix = sum_U s(U, z) e(U^c, u) and into Mmix = sum_U e(U, z) s(U^c, u) turns both into the sum over disjoint
  (W ∋ z, V ∋ u) of s(W, z) s(V, u). The same substitution gives M ≡ L directly.
- **(F3).** If s_z(G) = 0 then P_(z->u)(G) = 0. So the S-row is zero iff s_z = L_(z,u) = 0, iff G maps P into P.
  - Then P_(u->z)(G) is even by (F2), and M_(z,u)(G) ≡ L_(z,u)(G) = 0.
  - (E) then gives an even determinant. If z = u, M = 0. QED.

## Checks

- **(B).** `transitive_two_port.py` compares the closed form with SEAM's own `matrix()` (an independent construction:
  it contracts the prefix) for every k <= 8 and every port pair: 204 cases, 0 mismatches. The audit used three
  routes to k <= 9 (285 cases).
- **(D).** `parity_check.py` checks all 202,013 labelled rooted tournaments of order <= 6 and 300 random rooted
  tournaments of orders 7-10, with 0 violations. The audit extended this to 14.9 million rooted tournaments.
- **(E).** Checked on the 1,995 S-row-zero matrices in the per-order files of order <= 7. In all of them
  det ≡ M_(z,u) (mod 2), and M_(z,u)(G) and R_z(G) are even. The audit checked all 22,522 distinct S-row-zero
  matrices of order <= 8.
- **(F).**
  - **Session search.** An exhaustive search of all labelled tournaments of order 2-6 and all ordered pairs z ≠ u
    (1,004,340 rows) looked for mod-2 linear relations among (s_z, e_z, s_u, e_u, P_(z->u), P_(u->z), L, Lmix,
    Mmix). The relation space has dimension 3, spanned by M ≡ L, M ≡ P_(z->u) + P_(u->z) and Lmix ≡ Mmix (the
    audit's count; the session's first summary listed only two).
  - **Audit brute force.** No violations over all labelled tournaments of order <= 7 (89,084,724 rows), plus random
    ones of orders 8-10. The proof's two counting identities hold exactly on every G of order 2-5 (21,300 rows).

## What this does and does not give

- **Explicit family.** (B) is an explicit all-order family of composable pieces, and the only one here with a
  3-power entry.
- **Not a spectrum route.** Seam words built only from transitive gadgets are thin. Over k <= 10 and all port pairs
  they attain only about 3/4 of the odd numbers <= 2,000: 261 are missed, among them 11, 19, 35 and 37
  (FINITE-EXACT, `transitive_words_bfs.py`, `transitive_words_bfs.out`). The audit showed the 261 misses are
  unchanged for every k_max from 10 to 30, and that larger k cannot add values <= 2,000. They cannot carry the
  {7, 21} spectrum alone.
- **Parity.** (E) shows that a plane map of odd determinant could only come from a gadget with M_(z,u)(G) odd, and
  (F) shows no such gadget exists. So any self-covering argument on the plane S = 0 must align lattice cosets: no
  plane map is unimodular, and S78's descent lemma for unimodular maps is moot there.

## Credits

- SEAM.md (2026-10-29) for the transfer, its spanning-prefix construction and the singleton matrix.
- opus-S78 THM-4875 (B) for closed entries in the converse convention.
- Session29 (`session29_framed_profiles_20261022/session29-converse-congruences/REPORT.md` section 5, 2026-10-22)
  for (D).
- Session 193 (`session193_triangulated_proof_20261105/multiple_triangles/PROOF.md` section 2) for the exact
  identity whose mod-2 shadow is (D).
- Session 195 for the single-port transitive formulas.
- klein THM-4865 for the tower step T_1.
- The independent audit S11 for the Lmix ≡ Mmix relation and its proof, the direct proof of M ≡ L, and the
  equivalence in (F3).
