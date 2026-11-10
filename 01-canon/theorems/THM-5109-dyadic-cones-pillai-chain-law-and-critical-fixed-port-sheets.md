---
id: THM-5109
title: "Dyadic cones, the Pillai chain law and critical fixed-port sheets. (A)-(B) INHERITED (THM-4879): the one-arc response row and the closed transitive gadget table with its 3^(p-q-1) entry. (C) INHERITED mechanism (Session236 offsets section 2; THM-5085 (3)): a root response B = 2^j makes the near-join a pure cone h + 2^K for every K >= j + 1. New here: the census of the 11 strong classes of order 3-8 with two-path cover count Q a power of two (H = 5, 19 known; 217, 223, 599 new), the family 5^a 19^b 217^c 223^d 599^f + 2^K, and the window bound: since B >= h, every cone value lies in (2^K, 3 2^(K-1)], so cones never reach an odd N whose binary expansion begins 11. (D) With the new port at the exit inside the attached transitive gadget, a cone step re-lands on a dyadic port iff h + 3^(q+j-1) = 2^t with p = 2q + j, or (exit at the gadget source) 2^(p-1) h - e + 3^(p-2) 2^(j+1) is a power of two. (E) By Ellison's bound the only solutions of 1 <= 2^x - 3^y <= 2^y are (1,0), (2,1), (5,3), (8,5). The inside-gadget re-landings are 1 -> 9 and a family with e = 2^(p-1) h - d_p 2^(j+1), d_p = 1, 5, 13 (example 33 -> 289), and no two such steps chain. Re-landings at old vertices exist; their chains are OPEN. (F) On a post-bridge sheet the fixed-port letters S^j B tile the odd lattice off the line R = -c, so every point with R > 0 has a unique predecessor (mass exactly one; an elementary reformulation of known decoders)"
status: >
  PROVED (new): the window bound in (C); (D); (E).
  PROVED (elementary): (F), a reformulation of KNOWN decoders.
  NEW DATA (FINITE-EXACT): the census of strong classes of order 3-8 with Q a power of two; the product family is an
  immediate corollary of inherited ingredients.
  INHERITED, credited:
  - (A)-(B) are THM-4879 (A)-(B), this machine's own S11, in SEAM's converse convention (position i <-> m - i,
    S <-> e). The pre-audit draft typed them as new; MISTAKE-750.
  - (C)'s cone mechanism, with the same example 69, 133, 261, ..., is Session236 offsets/PROOF.md section 2 (seeds
    (5,0,32) and (9,2,16)) and this machine's THM-5085 (3) (S19 PROOF sections 2 (b) and 3), both on top of the
    near-join of HYP-9453 and THM-4865. The pre-audit draft typed it as new; MISTAKE-751.
  - The pairs (H, Q) = (5, 32) (Session227 multiroot section 5) and (19, 256) (Session228 unsaturated_module, Session229
    overlap).
  - Q is multiplicative over ordinal sums (Session227 multiroot section 5; Session229 dynamic eq. (6)).
  - The source response is the product of the blocks' Q (Session229 overlap section 1, eq. (1)).
  - B >= h (HYP-9453); R = 1 iff T is transitive and rooted at its sink (Session236).
  - (F)'s decoder content: Session197 ONE_CORE_BINARY_INVERSE section 3, THM-5096, THM-5085 (1)-(2). The "critical"
    reading has a quantitative precursor in Session206's (ln 2) X mass (HYP-9480 (C)/(D)).
  CITED INPUTS:
  - Ellison (1971): |2^x - 3^y| > 2^x e^(-x/10) for x >= 12, x not in {13, 14, 16, 19, 27} (Waldschmidt,
    arXiv:0908.4031, section 2.2.2);
  - THM-4810's singleton strong fibres F_s(5) and F_s(13).
  SCOPE: (D)-(E) concern the new dyadic port at the exit inside the attached transitive gadget.
  - Outside the gadget, one-step re-landings exist at old vertices:
    - TT_2 rooted at its sink, cone to 5, old source (5, 0, 32);
    - (1, 0, 4) -> 9 with (9, 2, 16);
    - (33, 0, 128) -> 289 with (289, 64, 512).
  - No two-step re-landing is known. The audit searched, with re-rooting allowed, up to outputs of order 15. OPEN.
  FINITE-EXACT checks:
  - Theorem B: h' and e' on 4,440 cases and R' on 3,768 (2,688 with p >= 2, 1,080 with p = 1), 0 mismatches; the
    audit's 7,956 rows, 0 mismatches;
  - e = O - E on 468 ports;
  - the dyadic census of order <= 8 (69 of 148 attained odd values <= 300; 181 dyadic triples);
  - the classification and the no-chain statement on all abstract triples with j <= 9 (author) and j <= 14 (audit,
    44.7 million states), and on the classified ones with j <= 17;
  - pre-registered P4.
  AUDIT: PASS WITH CORRECTIONS (5 major, 9 minor, all applied; MISTAKE-751),
  05-knowledge/results/kps_s20_dyadic_cones_pillai_20261109/audit/AUDIT.md.
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
  - THM-5085
  - THM-5096
  - THM-5107
  - THM-5070
  - HYP-9640
  - HYP-9480
  - THM-5040
  - THM-5055
  - THM-1370 (h-spectrum-omits-7-21-all-n)
  - 05-knowledge/results/session236_chanpak_offsets_orbits_20261109/offsets/PROOF.md
---

# THM-5109 - dyadic cones, the Pillai chain law and critical fixed-port sheets

**Notation.**
- A rooted tournament (T, v) has triple (h, e, R):
  - h = H(T);
  - e is the number of Hamiltonian paths ending at v;
  - R = B_v(T) is the root response, the sum of H(T - V(P)) over the paths P that start at v (THM-4865's B_t).
- In fixed-port notation, R = E + O and e = O - E.
- The one-arc join of (T, v) with a gadget (G, u, z) is N = T => G with u -> v. Its new root is z (SEAM.md).
- Q(X) = sum over S of H(X[S]) H(X[V - S]) is the two-path cover count (Session229; THM-5107). It is the response of a
  new source above X.

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

**The cone** (INHERITED mechanism: Session236 offsets section 2; THM-5085 (3); the near-join of HYP-9453 and
THM-4865). If R = 2^j, then h + 2^K is a Hamiltonian-path count for every K >= j + 1. This is the near-join with
p = K - j + 1.

**The window bound** (new). Since B >= h, a dyadic port has h <= 2^j <= 2^(K-1). So every cone value satisfies
2^K < h + 2^K <= 3 * 2^(K-1).
- An odd cone value therefore has a binary expansion beginning 10.
- No cone from any dyadic port of any tournament reaches an odd N in (3 * 2^(K-1), 2^(K+1)), for example 7, 13, 15,
  25-31, 49-63.

**Source-port families.**
- Inherited facts:
  - Q is multiplicative over ordinal sums (Session227, Session229);
  - a new source above X has response Q(X) (Session229), so it is a dyadic port exactly when Q(X) is a power of two.
- **New census.** Exactly 11 strong classes of order 3-8 have Q a power of two:

  | H | Q | classes |
  |---|---|---|
  | 5 | 2^5 | 1 |
  | 19 | 2^8 | 2 |
  | 217 | 2^12 | 2 |
  | 223 | 2^12 | 2 |
  | 599 | 2^13 | 4 |

  The pairs for 5 and 19 were known (Session227-229); 217, 223 and 599 are new. The trivial order-1 tournament also
  has Q = 2.
- **Product family.** 5^a 19^b 217^c 223^d 599^f + 2^K is a Hamiltonian-path count whenever
  K > 5a + 8b + 12c + 12d + 13f.

## (D) The chain law (new port at the exit inside the attached gadget)

Take a cone step (entry p >= 2) from a dyadic port (h, e, 2^j). The exit g_q is a dyadic port iff q < p and one of
the following holds:
- **(i)** q >= 2, p = 2q + j and h + 3^(q+j-1) = 2^t. Then R' = 2^(m-q+t).
- **(ii)** q = 1 and 2^(p-1) h - e + 3^(p-2) 2^(j+1) = 2^s.
  - When e = 0, this holds iff p = j + 2 and h + 3^j is a power of two.

## (E) Pillai rigidity (same scope)

- **(E1)** The only solutions of 1 <= 2^x - 3^y <= 2^y are (x, y) = (1, 0), (2, 1), (5, 3), (8, 5). For x >= 28
  Ellison gives 2^x - 3^y > 2^(0.855x) > 2^y; smaller x is checked exhaustively.
- **(E2)** Combined with B >= h, every dyadic cone step whose new port is the exit inside the gadget is one of:
  - 1 -> 9, from (1, 1, 1) with (p, q) = (4, 2), or from (1, 0, 2) with (p, q) = (3, 1);
  - a step with q = 1 and e > 0, p in {3, 5, 7}, and

        e = 2^(p-1) h - d_p 2^(j+1),    d_p = 1, 5, 13,    0 < e <= h <= 2^j.

    The realised example is 33 -> 289, with (h, e, B) = (33, 4, 64).

  The values 5 and 13, allowed by (E1), carry no suitable dyadic port:
  - each has a single strong fibre (THM-4810);
  - the port responses at H = 5 have odd parts {1, 3, 5, 7, 9}, with the dyadic ones at j >= 5;
  - at H = 13 the odd parts are {7, 9, 13, 15}.
- **(E3)** No dyadic port admits two consecutive such steps.
  - This also holds for abstract triples. The abstract first steps from 5 and 13 lead to 133 and 2061, and neither
    chains (PROOF.md).
  - Re-landings at old vertices are outside this scope (see Status); whether they can chain is OPEN.

## (F) Fixed-port sheets are critical (elementary reformulation of known decoders)

On the post-bridge sheet with h_u = c, consider the maps B: (h, R) -> (h + R - c, 2R - c) and S^jB:
(h, R) -> (h + 2^j R, 2^(j+1) R - c).
- They map the odd lattice {h, R odd} injectively onto the classes v_2(R + c) = 1 and v_2(R + c) = j + 1
  respectively.
- So their images tile the odd lattice minus the line R = -c, which has no predecessor.
- Every point with R > 0 has exactly one predecessor, and that predecessor has R > 0. Forward images can leave R > 0.
- The determinants are 2^(j+1), and the sum over j of 2^(-(j+1)) is 1.
- Predecessors of realised points need not be realised.

## Proofs and checks

Proofs are in [PROOF.md](../../05-knowledge/results/kps_s20_dyadic_cones_pillai_20261109/PROOF.md). Checks, censuses,
the pre-registered tests and the [independent audit](../../05-knowledge/results/kps_s20_dyadic_cones_pillai_20261109/audit/AUDIT.md)
are in the [source packet](../../05-knowledge/results/kps_s20_dyadic_cones_pillai_20261109/README.md).

## What is proved about cone inductions

- **Window bound.** By the window bound in (C), cones alone never reach odd numbers whose binary expansion begins 11.
  A cone induction in the style of the even-graph proofs THM-5040 and THM-5055 therefore needs other moves for half
  of every dyadic window.
- **Inside the gadget.** By (E3), cone steps whose new dyadic port is the exit inside a transitive gadget never chain.
- **Heuristic only.** Any one-arc construction reaching every large N must use root responses whose odd parts grow
  with N. This is conditional on HYP-9656 (A), which is OPEN.
