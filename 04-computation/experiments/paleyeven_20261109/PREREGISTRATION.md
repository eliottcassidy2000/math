# Pre-registration: mac-mini-2026-11-09-paleyeven

These predictions were written and committed before the computations that decide them were run. They are not edited
afterwards; outcomes are appended at the end.

## Objects

- **Paley tournament QR_p** (p prime, p = 3 mod 4): vertex set F_p, x -> y iff y - x is a nonzero square.
- **Paley double T_p**: vertices a_x, b_x (x in F_p). The pairs {a_x} and {b_x} each induce QR_p. Every arc between
  the halves goes a_x -> b_y (x != y), except the matched back arcs b_x -> a_x. So n = 2p, and the back arcs form a
  perfect matching.
- **Module defect md(T)**: the minimum, over vertex sets M with 2 <= |M| <= n - 1, of
  sum_{z not in M} min(|N+(z) cap M|, |N-(z) cap M|). md(T) equals THM-4985's decomposability defect dec(T) for n >= 3.
- **Escaper** (THM-4985): strong, prime, lambda >= 3, H(T) in S(n - 1), and H(T - v) in S(n - 2) for every v.
  Here S(m) is the set of Hamiltonian-path counts of m-tournaments.
- **Dominating-set count** d(X): the number of vertex sets S of a simple graph X with every vertex in S or adjacent to S.
  It is always odd (Brouwer 2009).
- **Royle-even graph**: no automorphism reverses an odd number of edges relative to a fixed reference orientation
  (Royle-Praeger-Glasby-Freedman-Devillers 2023: even graphs and tournaments are equinumerous; THM-4531).

## Predictions

- **Q1** (p = 7, n = 14): md(T_7) = lambda(T_7) = minimum semidegree = 4, and T_7 is an escaper of order 14.
  Predicted YES, probability 0.70. (Previous record: decomposability defect 3 at order 12, THM-4985 (2).)
- **Q2** (p = 11, n = 22): md(T_11) = lambda(T_11) = 6. Predicted YES, probability 0.80. Separately, T_11 is certified
  as an escaper within this session: probability 0.35.
- **Q3** (all p): md(T_p) = (p + 1)/2 for every prime p = 3 mod 4 with p >= 7. Predicted true, probability 0.75; a
  complete proof found this session, probability 0.60.
- **Q4** (dominating sets on even graphs): among the odd numbers <= 255, the values d(X) that are attained by NO
  Royle-even graph X are exactly {3, 7, 21}. Predicted YES, probability 0.30. Weaker parts:
  - 21 is attained by no Royle-even graph: probability 0.60;
  - some odd value <= 255 other than 3, 7, 21 is unattained: probability 0.55.
  (3 and 7 are expected to be unattained, by a partial hand argument made before any computation. Connected graphs on
  at most 4 vertices give d = 3 only for K_2 and d = 7 only for K_3. Ore's bound d >= 2^(m/2) leaves the 5-vertex
  connected graphs for d = 7, which were not checked by hand. K_2 and K_3 have true twins, hence an odd transposition.)

## Outcomes

(Appended after the runs, 2026-11-09. The predictions above are unedited.)

- **Q1 REFUTED** (predicted YES, 0.70). md(T_7) = lambda(T_7) = minimum semidegree = 4 holds (`paley_mdefect.out`).
  T_7 is not a certifiable escaper: H(T_7) = 4,730,397 exceeds 3,719,831, the best known order-13 count (THM-329),
  and its deletions (686,457) exceed 531,205, the best known order-12 count (`p7_tprop.out`).
- **Q2 CONFIRMED** for md(T_11) = lambda(T_11) = 6 (probability 0.80).
  - The separate prediction (certification of T_11 as an escaper within the session, 0.35) was NOT achieved.
  - Its values are 0.363 E[H(21)] and 0.393 E[H(20)], beyond witness search.
  - It is an escaper only conditionally (THM-4995 (3)).
- **Q3 CONFIRMED and PROVED** (0.75 / 0.60): md(D(B)) = (q + 1)/2 for every doubly regular B of order q >= 7
  (PROOF.md, Theorem A2). This covers every Paley prime p >= 7.
- **Q4 REFUTED** (predicted YES, 0.30). The exceptional set is {3, 7, 13}, exact below 512 (`dom_K511_summary.txt`).
  - "21 is unattained" (0.60) is REFUTED: 21 = d(K_1 v (P_3 u K_1)), which is even.
  - "Some other odd value <= 255 is unattained" (0.55) is CONFIRMED: 13, realised only by K_4 - e, which is odd.
- **Not pre-registered, found along the way:**
  - the order-18 double of C_9{1,2,3,5}, a certified escaper of defect 5 (`wit_d9_del_m16.out`, `wit_d9_H_m17.out`);
  - the odd-graph exceptional set {1, 5, 37} below 512.
