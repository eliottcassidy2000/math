---
id: HYP-9575
title: "(A) The Paley double T_11 = D(QR_11) is an escaper of order 22 with simplicity index 6, and every double of a doubly regular tournament of order q >= 11 is an escaper. (B) Even-graph Chan-Pak: the dominating-set count d takes every odd value on even graphs (von Brömssen's notion) except 3, 7 and 13, and is concise there. (C) The least order of an escaper with simplicity index >= 4 is between 13 and 16"
status: "OPEN. Audit P1 PASS WITH CORRECTIONS (applied; MISTAKE-746). (A) PROVED conditional on (S_c) for large q (THM-4995 (3)); T_11 needs only (S_0.4) at orders 20 and 21. (B) FINITE-EXACT below 512 (THM-4996 (3)); constructive evidence to 2^22, with 12,326 unreached values clustered near 2^m - 2^k. (C) The upper bound 16 is a certified escaper of simplicity index 4 (THM-4995 (4b)); for index >= 5 the bound is 18 (THM-4995 (4a)). No lower bound is proved; 13 would follow from k*(12) = 3 (HYP-9565 (C))."
session: mac-mini-2026-11-09-paleyeven
date: 2026-11-09
related: [THM-4995, THM-4996, HYP-9565, THM-4985, HYP-9481, THM-4531, THM-1370]
script: 04-computation/experiments/paleyeven_20261109/
---

# HYP-9575 - Paley doubles as escapers, and the even-graph Chan-Pak

Packet: [note](../results/paleyeven_20261109/README.md), [PROOF](../results/paleyeven_20261109/PROOF.md),
[audit](../results/paleyeven_20261109/AUDIT_P1.md).

md is the simplicity index of Müller-Pelant (1974): the fewest arc reversals that create a nontrivial module.

## (A) Doubles are escapers

THM-4995 (1) proves md = lambda = (q + 1)/2 for doubles D(B) of doubly regular B. Escaper status needs H(D(B)) in
S(2q - 1) and every H(D(B) - v) in S(2q - 2).
- For T_11 these values are 17,667,825,987,845 = 0.363 E[H(21)] and 1,822,800,617,515 = 0.393 E[H(20)].
- Both lie far inside the range the exact censuses show to be solid at orders 10 and 11 (up to about 2 E[H(m)]).
- Certifying them needs witnesses of order 21 and 20, beyond wit.c's reach. The hit rate is about 10^-7 per base
  tournament.

Conjecture: T_11 is an escaper, and so is D(B) for every doubly regular B of order q >= 11.

T_7 is not covered and not certifiable: its count 4,730,397 exceeds the best known order-13 count 3,719,831, and its
deletions 686,457 exceed the best known order-12 count 531,205. It is an escaper iff 4,730,397 is in S(13) and
686,457 is in S(12), and neither is known.

## (B) Even-graph Chan-Pak

d(G) is odd for every graph (Brouwer-Csorba-Schrijver). Even graphs in von Brömssen's sense are equinumerous with
tournaments (RPGFD).

Conjecture: {d(G) : G even} = odd numbers minus {3, 7, 13}, and every value k is attained on O(log k) vertices. This
is the analogue of Chan-Pak Conjecture 6.2 (exceptions {7, 21}) for an intrinsic odd count on the
tournament-equinumerous side.
- **Exact below 512.** The exceptions are exactly {3, 7, 13}; their only realisers K_2, K_3 and K_4 - e have adjacent
  twins. The complementary class of odd graphs misses exactly {1, 5, 37} below 512.
- **Constructions.** Unions, cones K_1 v (G u E_j), E_r-joins and general joins of padded even graphs reach every
  odd value <= 2^22 except about 0.6%. The unreached values cluster near 2^m - 2^k (including 2^m). A proof needs
  families reaching those, for instance:
  - dense even graphs (complements of perfect matchings: d = 2^m - 1 - m);
  - partial hubs. For a hub z with pendant leaves and gadgets (H_i, W_i) adjacent to z on W_i,
    d = prod A_i + prod d(H_i) - prod d_0(H_i, W_i).
- **Conciseness.** By Chan-Pak Proposition 5.14 it would follow from almost completeness for connected even graphs,
  since d >= 2^(n/2) on connected graphs (Ore).

## (C) Least order of an escaper with simplicity index >= 4

- **Upper bound 16 (certified).** D(B8) for a near-regular 8-tournament B8 (THM-4995 (4b)). For index >= 5,
  THM-4995 (4a) gives order 18.
- **Order 15 is not covered by the constructions tried.**
  - The R_5-chain and the unequal doubles (8-block plus QR_7) reach md 4.
  - Each has exactly one deletion above the best known order-13 count 3,719,831 (largest 4.15-4.35·10^6), so they
    are not certifiable.
  - Not excluded.
- **Lower bound.** An escaper has md <= lambda. If every 4-arc-strong 12-tournament has H > 95,095 = max H(11) (that
  is, k*(12) = 3, HYP-9565 (C)), then order 13 is a lower bound.

Conjecture: the least order is 15 or 16.
