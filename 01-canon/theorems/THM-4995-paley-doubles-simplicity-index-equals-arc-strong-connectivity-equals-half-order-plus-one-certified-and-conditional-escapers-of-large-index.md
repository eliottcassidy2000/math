---
id: THM-4995
title: "Paley doubles. For a doubly regular tournament B of order q >= 7, let D(B) be two copies of B with every arc going from the first to the second, except a matched set of back arcs b_x -> a_x. md is the simplicity index (fewest arc reversals to a non-simple tournament; Müller-Pelant 1974; THM-4985's decomposability defect). (1) PROVED: md(D(B)) = lambda(D(B)) = delta0(D(B)) = (q + 1)/2; D(B) is strong and prime; D(QR_p) is self-converse. The input is Müller-Pelant's md(B) = (q - 1)/2 for doubly regular B, re-proved here from S S^T = qI - J by a Parseval bound. (2) FINITE-EXACT: T_7 = D(QR_7) has md 4 but H = 4,730,397, above the best known order-13 count, so it is not certifiable; T_11 has md 6 and H = 17,667,825,987,845, with all deletions 1,822,800,617,515. (3) PROVED-CONDITIONAL: under (S_c) (all odd values up to c m!/2^(m - 1) except 7 and 21 occur at order m), D(B) is an escaper for all large q, so escapers of unbounded simplicity index exist (HYP-9565 (A)); with c = 1 every q >= 55 works (Brégman's bound, Alon's method). (4) FINITE-EXACT: certified escapers of simplicity index 5 at order 18 (the double of C_9{1,2,3,5}) and 4 at order 16 (the double of a near-regular 8-tournament); the previous record was 3"
status: >
  PROVED + FINITE-EXACT + PROVED-CONDITIONAL. Independent audit P1 PASS WITH CORRECTIONS (1 MAJOR: the simplicity-index
  priority; 13 MINOR); all applied; MISTAKE-746. The order-16 example (4b) was added after audit P1 and checked in the
  follow-up audit (AUDIT_P1.md, addendum: (4b) CONFIRMED in full; Finding A1 on the order-15 data applied).
  - (1) PROVED (PROOF.md Part A).
    - Lemma A0 (md <= lambda <= delta0): md <= min(delta0, twin distance) is BLT 2022, Proposition 3.1. md <= lambda
      holds because a non-strong tournament is not simple.
    - Lemma A1 (md(B) = (q - 1)/2 for doubly regular B) is the equality case of Müller-Pelant 1974. The Parseval proof
      here is an alternative.
    - Theorem A2 (a thirteen-case analysis) is new.
    - FINITE-EXACT checks: md(QR_p) = lambda = (p - 1)/2 for p = 7, 11, 19, 23; md, twin distance, lambda and delta0
      of T_7 and T_11 over all 2^14 and 2^22 vertex sets (mdefect.c; reproduced by the audit's tstruct.c).
  - (2) FINITE-EXACT: subset DP (tprop.c; reproduced by the audit with a subset DP and with inclusion-exclusion).
    - T_7 is an escaper iff 4,730,397 is in S(13) and 686,457 is in S(12), and neither is known.
    - It is not one if a(12) = 531,205 (THM-329's best order-12 count, "strongly believed" maximal) or if
      max H(13) < 4,730,397.
  - (3) PROVED-CONDITIONAL on (S_c).
    - Brégman's permanent bound (lgamma, checked against exact factorials to 9e-13) covers 55 <= q <= 20,000
      (bregman_table.out).
    - For every q >= 7, both natural-log ratios satisfy R(q) <= q log(3/4) + 2.85 log q + 4.8, which is negative for
      q >= 57 (Robbins' Stirling bounds; derived and checked to q = 2·10^6 by audit P1, bregman_tail.out).
    - (S_c) is OPEN. klein's exact censuses show S(10) and S(11) solid to 2.07 and 2.06 times m!/2^(m - 1) (THM-4889).
  - (4) FINITE-EXACT. md, lambda and delta0 over all vertex sets; H by subset DP. Every membership has an explicit
    witness (wit.c), re-counted by wit.c's direct DP, by tprop.c, and by the audit's two algorithms.
  KNOWN and credited:
  - J. Müller, J. Pelant, On strongly homogeneous tournaments, Czechoslovak Math. J. 24 (1974) 378-391: the
    simplicity index s(T) <= (n - 1)/2, with equality iff T is doubly regular;
  - H. Boussaïri, A. Lakhlifi, M. Talbaoui, The simplicity index of tournaments, Contrib. Discrete Math. 17 (2022)
    172-186, arXiv:1907.11777 (Proposition 3.1: s(T) <= min(delta0, twin distance));
  - S. Kirkland (1993): lambda as the reversal index;
  - K. B. Reid, E. Brown, Doubly regular tournaments are equivalent to skew Hadamard matrices, J. Combin. Theory
    Ser. A 12 (1972) 332-338;
  - L. M. Brégman (1973): the permanent bound;
  - N. Alon, The maximum number of Hamiltonian paths in tournaments, Combinatorica 10 (1990) 319-324: Brégman's bound
    applied to Hamiltonian paths;
  - THM-4985 / HYP-9565 (mac-mini-2026-11-09-genesis): escapers and the Paley-halves heuristic; their "decomposability
    defect" is the simplicity index (credit added there after audit P1);
  - Session230 (primes §5): the bounded-defect question;
  - THM-4811 (B): lambda versus minimum semidegree;
  - THM-329 (opus-2026-05-27-S6): best known counts at orders 12 and 13;
  - HYP-9365 / THM-4757: H(QR_p) 2^(p - 1)/p! -> e;
  - wit.c and the attachment form (HYP-2879, THM-4097, THM-4114; THM-4942 (1)).
  NEW (no prior art found by audit P1):
  - the doubles D(B) and Theorem A2;
  - the conditional escapers of unbounded simplicity index;
  - the certified escapers of (4).
  NOT CLAIMED:
  - that T_7 is a first creator: it is one iff 4,730,397 is not in S(13), which holds if max H(13) < 4,730,397, and
    max H(13) is unknown;
  - that T_11 is an escaper (only conditional);
  - any statement about non-doubly-regular blocks beyond the computed examples.
session: mac-mini-2026-11-09-paleyeven
date: 2026-11-09
related: [THM-4985, HYP-9565, HYP-9575, THM-4996, THM-329, THM-4811, THM-4889, HYP-9481]
note: 05-knowledge/results/paleyeven_20261109/README.md
proof: 05-knowledge/results/paleyeven_20261109/PROOF.md
audit: 05-knowledge/results/paleyeven_20261109/AUDIT_P1.md
scripts:
  - 04-computation/experiments/paleyeven_20261109/paley_build.py
  - 04-computation/experiments/paleyeven_20261109/chain_build.py
  - 04-computation/experiments/paleyeven_20261109/mdefect.c
  - 04-computation/experiments/paleyeven_20261109/wit.c (copy of genesis_20261109/wit.c)
  - 04-computation/experiments/genesis_20261109/tprop.c
  - 04-computation/experiments/paleyeven_20261109/bregman_table.py
  - 04-computation/experiments/paleyeven_20261109/d8_explore.py, d87_explore.py, verify_md_py.py
outputs:
  - paley_mdefect.out, qr_mdefect.out, d9_mdefect.out, chain_mdefect.out, circ_mdefect.out
  - p7_tprop.out, p11_tprop.out, d9_tprop.out, chain_tprop.out
  - wit_d9_del_m16.out, wit_d9_H_m17.out, bregman_table.out
  - d8_mdefect.out, d8_md4_tprop.out, d8_144.txt, d8_144_tprop.out, wit_d8_H_m15.out, wit_d8_del_m14_part*.out
  - d87_mdefect.out, d87_md4_tprop.out
witness: "Order 18: 773,858,997 in S(16) and 6,450,103,269 in S(17). Order 16: 144,227,835 in S(15) and the eight deletion values in S(14). All certified."
---

# THM-4995 - Paley doubles: simplicity index (q + 1)/2 and escapers of large simplicity index

## Statement

Let B be a tournament on [q]. Its **double** D(B) has vertices a_x, b_x (x in [q]).
- The sets {a_x} and {b_x} each induce a copy of B.
- Every arc between the copies goes a_x -> b_y (x != y), except the back arcs b_x -> a_x.
- For a Paley tournament QR_p (p prime, p = 3 mod 4), T_p = D(QR_p).

md(T) is the **simplicity index** s(T) of Müller-Pelant (1974): the least number of arc reversals that create a module
M with 2 <= |M| <= n - 1. Equivalently it is the least, over such M, of
sum_{z not in M} min(|N+(z) cap M|, |N-(z) cap M|). THM-4985 called it the decomposability defect.

1. **(PROVED)** If B is doubly regular of order q >= 7, then md(D(B)) = lambda(D(B)) = delta0(D(B)) = (q + 1)/2.
   - The input is md(B) = (q - 1)/2 (Müller-Pelant), together with md(T) <= lambda(T) <= delta0(T) for every
     tournament on >= 3 vertices.
   - D(B) is strong and prime, and Aut(B) acts on it diagonally.
   - D(B) is self-converse when B is; in particular T_p is.
   - For q = 3 the statement fails: md(D(C_3)) = 0.
2. **(FINITE-EXACT)**
   - T_7: md = twin distance = lambda = delta0 = 4. H(T_7) = 4,730,397 and H(T_7 - v) = 686,457 for every v.
   - T_11: md = lambda = delta0 = 6, twin distance 7. H(T_11) = 17,667,825,987,845 and H(T_11 - v) = 1,822,800,617,515.
   - The matched back arcs multiply H(QR_p)^2 by 132.4 and 1953.7.
   - H(T_7) exceeds 3,719,831, the best known order-13 count, and its deletions exceed 531,205, the best known order-12
     count (THM-329). So T_7 is not a certifiable escaper. Whether it is one is undecided.
3. **(PROVED-CONDITIONAL)** Let (S_c) be: for all m >= m_0, every odd h <= c m!/2^(m - 1) other than 7 and 21 is in
   S(m).
   - Under (S_1), D(B) is an escaper for every doubly regular B of order q >= 55 with 2q - 2 >= m_0.
   - Under (S_c) for any fixed c > 0 the same holds for all large q.
   - Hence escapers of arbitrarily large simplicity index exist under (S_c). This settles HYP-9565 (A)
     conditionally.
   - T_11 needs only (S_0.4) at orders 20 and 21: its ratios are 0.3626 and 0.3928.
4. **(FINITE-EXACT, certified escapers)**
   - **(a) Order 18, md 5.** Let C = C_9{1,2,3,5} (x -> x + s for s in {1, 2, 3, 5}).
     - D(C) has md = lambda = delta0 = 5 and is strong and prime.
     - H(D(C)) = 6,450,103,269 is in S(17), and H(D(C) - v) = 773,858,997 is in S(16) for every v.
     - So D(C) is an escaper of order 18 with simplicity index 5. The previous record was 3 (THM-4985 (2)).
   - **(b) Order 16, md 4.** Let B8 be the near-regular 8-tournament with code 1010100100111001101001010101 (upper
     triangle, '1' iff i -> j; scores 3, 4, 3, 4, 4, 3, 3, 4; twin distance 2).
     - D(B8) has md = lambda = delta0 = 4 and is strong and prime.
     - H(D(B8)) = 144,227,835 is in S(15).
     - Its eight distinct deletion values 16,750,389, 17,527,823, 18,351,387, 18,908,079, 19,218,519, 19,720,143,
       20,574,141 and 21,469,761 are in S(14).
     - So D(B8) is an escaper of order 16 with simplicity index 4.
5. **(FINITE-EXACT, negative or undecided)**
   - Denser back-arc patterns lower md: patterns b_x -> {a_x, a_{x+d}} give md 2, and the Fano incidence pattern
     gives md 0.
   - The three-block Paley chain (order 21) has md = lambda = 4 and H = 110,624,292,765.
   - The order-15 chain of three R_5 blocks with back arcs {0, 2} has md 4, H = 29,678,025 = 2.79 E[H(14)] and a
     deletion of 5,081,917.
   - Order-15 unequal doubles (near-regular 8-block plus QR_7) with md 4 (14 of a random 400): each has exactly one
     deletion above 3,719,831, with largest deletions 4.15-4.35·10^6. The number was corrected after the audit
     addendum (Finding A1).
   - In both order-15 cases a deletion exceeds the best known order-13 count. They are not certifiable, but not
     excluded.

## Proof

PROOF.md Part A (Lemma A0, Lemma A1, Theorem A2, Proposition A3, Theorem A4, Theorem A5).
