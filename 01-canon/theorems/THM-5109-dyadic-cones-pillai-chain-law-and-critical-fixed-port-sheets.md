---
id: THM-5109
title: "Dyadic cones, the Pillai chain law and critical fixed-port sheets. (A)-(B) INHERITED: the one-arc response row as path sums and the closed transitive gadget table with its 3^(p-q-1) entry are THM-4879 (A)-(B) (this machine's S11, SEAM converse convention), restated in the direct convention. (C) A root response B = 2^j (a dyadic port) makes the near-join a pure cone, h + 2^K for every K >= j + 1; the strong tournaments of order <= 8 with two-path cover count Q a power of two have H = 5, 19, 217, 223, 599, so for example 5^a 19^b + 2^K is a Hamiltonian-path count for K >= 5a + 8b + 1. (D) A cone step from a dyadic port again ends at a dyadic port iff h + 3^(q+j-1) = 2^t with p = 2q + j, or (exit at the gadget source) 2^(p-1) h - e + 3^(p-2) 2^(j+1) is a power of two. (E) By Ellison's bound the only solutions of 1 <= 2^x - 3^y <= 2^y are (1,0), (2,1), (5,3), (8,5); dyadic cone steps are 1 -> 9 and one explicit e-relation family (example 33 -> 289), and no two dyadic cone steps chain. (F) On a post-bridge sheet the fixed-port letters S^j B tile the odd lattice with unique predecessors (mass exactly one)"
status: >
  PROVED (new): (C)-(F).
  INHERITED (restated in the direct convention): (A)-(B) are THM-4879 (A)-(B), this machine's own S11, in SEAM's converse
  convention: H-row (1, -[u = k-1], 2^(k-1-u)), R-row (2^z, -mu, c) with c = 2^(k-z+u) 3^(z-u-1) for z > u. Under the
  converse, position i becomes m - i and S becomes e, and the entries agree with (B) below. Also SEAM.md and THM-4875 (B);
  Session 195 for the single-port case. The pre-audit draft of this file typed (A)-(B) as new; MISTAKE-750.
  - Cited inputs: Ellison (1971) for |2^x - 3^y| > 2^x e^(-x/10), x >= 12 except 13, 14, 16, 19, 27 (as quoted by
    Waldschmidt, arXiv:0908.4031, section 2); and THM-4810's singleton strong fibres F_s(5) and F_s(13) for the
    exclusion of the values 5 and 13 in (E2).
  - The h = 13 part of that exclusion also uses the order-<= 8 port census, which contains the source over the unique
    order-5 strong class with H = 13.
  INHERITED and credited:
  - the one-arc join law and its first two rows (SEAM.md; THM-4875 (B));
  - near-joins h + 2^k B and B >= h (HYP-9453; THM-4865); (C) is their specialisation to B = 2^j;
  - the two-path cover count Q (Session229 proper-graft row; THM-5107, with 4 | Q);
  - the fixed-port letters, closed form and decoders (THM-5085, THM-5096); (F) recasts the decoder as an exact tiling.
  FINITE-EXACT (checks):
  - Theorem B against brute force on 4,440 cases, with 0 mismatches;
  - e = O - E on 468 ports;
  - the dyadic census of order <= 8 (69 of 148 attained odd values <= 300 have a dyadic port; 181 distinct dyadic
    triples);
  - the classification and the no-chain statement on all abstract triples with j <= 9, and on the classified ones
    with j <= 17;
  - pre-registered P4 (no two-step chain among the realised bases of order <= 8, m <= 20).
  NOT CLAIMED: Chan-Pak, binary or ternary (OPEN).
session: kind-pasteur-2026-11-09-S20
source: 05-knowledge/results/kps_s20_dyadic_cones_pillai_20261109/README.md
proof: 05-knowledge/results/kps_s20_dyadic_cones_pillai_20261109/PROOF.md
related:
  - HYP-9656
  - HYP-9453
  - THM-4865
  - THM-4875
  - THM-4879
  - HYP-9471
  - THM-4810
  - THM-5107
  - THM-5085
  - THM-5096
  - HYP-9640
  - HYP-9480
  - THM-5040
  - THM-5055
  - THM-1370 (h-spectrum-omits-7-21-all-n)
---

# THM-5109 - dyadic cones, the Pillai chain law and critical fixed-port sheets

**Notation.**
- A rooted tournament (T, v) has triple (h, e, R):
  - h = H(T);
  - e is the number of Hamiltonian paths ending at v;
  - R = B_v(T) is the root response, the sum of H(T - V(P)) over the paths P that start at v (THM-4865's B_t).
- In fixed-port notation, R = E + O and e = O - E.
- The one-arc join of (T, v) with a gadget (G, u, z) is N = T => G with u -> v. Its new root is z (SEAM.md).
- Q(X) = sum over S of H(X[S]) H(X[V - S]) is the two-path cover count, the response of a new source above X
  (Session229; THM-5107).

## (A) The response row (INHERITED: THM-4879 (A), converse convention)

    R(N, z) = B_z(G) h - J_1 e + (J_2 + K) R,

where:
- J_1 = sum over paths C of G from z avoiding u of s_u(G - V(C));
- J_2 + K = sum over the same C of E_u(G - V(C)) + sum over paths C of G from z ending at u of Q(G - V(C)).

Here B_z(G) >= 1, so the old value h always enters the new response.

## (B) Transitive gadgets (INHERITED: THM-4879 (B), converse convention)

Let G = TT_m, with entry at position p and exit at position q. Then
- h' = h + 2^(p-1) R - [p = 1] e;
- e' = [q = m] h - [p = 1][q = m] e + E_uz R, where E_uz = 2^(p-1) if q = m > p, 2^(q-1) if q < p = m, and 0
  otherwise.

For p >= 2, R' is:
- 2^(m-q) h' if q >= p;
- 2^(m-q) h + 3^(p-q-1) 2^(m-p+q) R if 2 <= q < p;
- 2^(m-1) h - 2^(m-p) e + 3^(p-2) 2^(m-p+1) R if q = 1.

For p = 1, R' = 2^(m-1)(h + R) if q = 1, and 2^(m-q) h' if q >= 2.

## (C) Dyadic cones

**The cone.** If R = 2^j, then h + 2^K is a Hamiltonian-path count for every K >= j + 1. This is the near-join with
p = K - j + 1.

**Source-port families.**
- Q is multiplicative over ordinal sums. A new source above X is a dyadic port exactly when Q(X) is a power of two.
- Among strong tournaments of order <= 8, exactly 11 classes have this property, with H in {5, 19, 217, 223, 599}
  and Q = 2^5, 2^8, 2^12, 2^12, 2^13.
- Hence 5^a 19^b 217^c 223^d 599^f + 2^K is a Hamiltonian-path count whenever K > 5a + 8b + 12c + 12d + 13f.
- The cones exist, but they cannot be chained ((E)).

## (D) The chain law

Take a cone step (entry p >= 2) from a dyadic port (h, e, 2^j). The exit g_q is a dyadic port iff q < p and one of
the following holds:
- **(i)** q >= 2, p = 2q + j and h + 3^(q+j-1) = 2^t. Then R' = 2^(m-q+t).
- **(ii)** q = 1 and 2^(p-1) h - e + 3^(p-2) 2^(j+1) = 2^s.
  - When e = 0, this holds iff p = j + 2 and h + 3^j is a power of two.

## (E) Pillai rigidity

- **(E1)** The only solutions of 1 <= 2^x - 3^y <= 2^y are (x, y) = (1, 0), (2, 1), (5, 3), (8, 5). For x >= 28
  Ellison gives 2^x - 3^y > 2^(0.855x) > 2^y; smaller x is checked exhaustively.
- **(E2)** Combined with B >= h, every dyadic cone step is one of:
  - 1 -> 9, from (1, 1, 1) with (p, q) = (4, 2), or from (1, 0, 2) with (p, q) = (3, 1);
  - a step with q = 1 and e > 0, p in {3, 5, 7}, and

        e = 2^(p-1) h - d_p 2^(j+1),    d_p = 1, 5, 13,    0 < e <= h <= 2^j.

    The realised example is 33 -> 289, with (h, e, B) = (33, 4, 64).

  The values 5 and 13, allowed by (E1), carry no suitable dyadic port. Each has a single strong fibre (THM-4810). The
  port responses at H = 5 have odd parts {1, 3, 5, 7, 9}, with the dyadic ones at j >= 5. At H = 13 the odd parts are
  {7, 9, 13, 15}.
- **(E3)** No dyadic port admits two consecutive dyadic cone steps with the new ports inside the attached gadgets.

## (F) Fixed-port sheets are critical

On the post-bridge sheet with h_u = c, the maps B: (h, R) -> (h + R - c, 2R - c) and S^jB:
(h, R) -> (h + 2^j R, 2^(j+1) R - c) have these properties:
- they send the lattice {h, R odd} bijectively onto the classes v_2(R + c) = 1 and v_2(R + c) = j + 1 respectively;
- every point therefore has exactly one predecessor;
- the determinants are 2^(j+1), and the sum over j of 2^(-(j+1)) is 1.

## Proofs and checks

Proofs are in [PROOF.md](../../05-knowledge/results/kps_s20_dyadic_cones_pillai_20261109/PROOF.md). Checks, censuses
and the pre-registered tests are in the [source packet](../../05-knowledge/results/kps_s20_dyadic_cones_pillai_20261109/README.md).

**Consequence (heuristic, see HYP-9656).**
- Pure cone inductions, the mechanism of the even-graph proofs THM-5040 and THM-5055, cannot be run in tournaments.
- Any one-arc construction that reaches every large N must use root responses whose odd parts grow with N.
