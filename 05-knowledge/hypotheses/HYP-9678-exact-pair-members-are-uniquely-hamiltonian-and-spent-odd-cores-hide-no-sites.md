---
id: HYP-9678
title: "Exact slot pairs and hidden transfer sites. (A) REFUTED (pre-registered S3): order-10 strong tournaments with five Hamiltonian cycles (H = 333) share the full slot polynomial of C3 => Y7. (A') Every strong member of a DETACHABLE pair (Q_X = Q_(Y -> TT1)) is uniquely Hamiltonian. (B) Detachable exact pairs, strong X of even order with Q_X = t Q_Y for a strong Y, exist at infinitely many orders; known at 6, 8, 10 with H = 15, 95, 317. (C) Capacity is additive: the maximum number of singleton components over all tournaments with slot polynomial Q_T is s(T) + d(T), where d(T) counts the even strong components of T that have a detachable partner. In particular a sum of odd strong cores hides no site; in particular no odd merge Q(Y1 => Y2) = Q(Z -> TT1) exists. (D) Consequence of (C) for returns: in the Pi clock a neutral rebase of the right module adds at most d sites, so an unbounded return must bring in a fresh detachable even component by a non-neutral move before each refill"
status: >
  OPEN.
  (A): REFUTED. Through order 9 all 7 strong members have C = 1 (complete census). At order 10, pre-registered S3
  searched every strong tournament with at most 5 backward arcs. It found members with C = 5 and H = 333 (partner
  C3 => Y7, Y7 = classgen code 43104), at least 2 isomorphism classes; THM-5151 PROOF Proposition 5(d).
  (A'): OPEN. Every detachable member known has C = 1: orders <= 9 complete, and in S3's range the only detachable
  profile is H = 317, with C = 1.
  (B): FINITE-EXACT pairs at orders 6, 8, 10 (THM-5151 (4), (5)). Order 12 is undecided:
  - pre-registered S1 found no partner among strong order-11 tournaments with at most 5 backward arcs (prediction
    FAILED);
  - a biclique follow-up was also empty;
  - no (3,9) or (5,7) block partner exists within C = 1.
  (C): FINITE-EXACT for every tournament of order <= 9 (198,946 classes, 52,432 profile groups, no violation;
  THM-5151 PROOF Proposition 3'). The odd-merge case is checked completely for |Y1| + |Y2| <= 10. At order 11,
  pre-registered S2 CONFIRMED no match among the 3,847,591 tournaments with at most 5 backward arcs, for all
  (5,7) block pairs.
  (D) is a logical consequence of (C) together with THM-5151 (2) and (9). The unconditional part, that the
  transfer count is at most cap(Q_(N_0)), is THM-5151 (9). It is not a statement about other frames or about
  non-neutral moves.
session: kind-pasteur-2026-11-13-S26
source: 05-knowledge/results/kps_s26_exact_pair_returns_20261113/README.md
related:
  - THM-5151
  - 05-knowledge/results/session310_positive_clocks_20261113/NEUTRAL_TAILS.md
  - 05-knowledge/results/session309_changing_module_returns_20261113/FIXED_BANKS.md
  - 05-knowledge/results/session261_creative_chanpak_20261111/core/PACKING_GRAMMARS.md
---

# HYP-9678 - exact pair members are uniquely Hamiltonian; spent odd cores hide no sites

Q_T(t) = sum_j F_j(T) binom(t, j). Let s(T) count singleton strong components. Call an even strong X
*detachable* if Q_X = Q_(Y -> TT1) for some odd strong Y.

**(A) C = 1 law: REFUTED.** The claim was that if X is strong, Z is not, and Q_X = Q_Z, then X has exactly one
directed Hamiltonian cycle. It holds through order 9, where the census is complete (7 classes). At order 10 it
fails:

    X = (1002,1020,1017,1008,481,704,896,768,544,16)

has C(X) = 5, H = 333 and Q_X = Q(C3 => Y7), where Y7 is classgen code 43104 (H = 111, C = 11). This came from
pre-registered search S3.

**(A') Detachable C = 1 law (OPEN).** If X is strong and Q_X = Q_(Y -> TT1), then C(X) = 1. Every known detachable
member (H = 15, 95, 317) satisfies it, and S3 found no detachable counterexample with at most 5 backward arcs at
order 10.

**(B) Infinitely many detachable pairs.** Known detachable members:

- order 6, H = 15, the one-apex O(01001);
- order 8, H = 95, Session 310's B = O(0011111);
- order 10, H = 317 (prime), not one-apex.

Block-type members, Q_X = Q(Y1 => Y2) with both Y_i odd strong, are known at H = 45 (order 8) and H = 153, 219
(order 10).

**(C) Additive capacity.**

    max { s(T') : Q_T' = Q_T } = s(T) + d(T),

where d(T) counts the detachable even strong components of T. Equivalently, in THM-5151 (2), mu is additive over
strong components. The special case of two odd cores is the no-odd-merge statement:
Q(Y1 => Y2) ≠ Q(Z -> TT1).

**(D) Consequence for returns.** Assume (C). In the Π clock, a neutral rebase of the right module N can expose at
most s(N) + d(N) trailing sites. After a run the right module is a sum of odd cores, so it has d = 0. A return
rule must therefore re-insert a detachable even component by a non-neutral move, such as Session 309's migration
or recharge, or by a change of frame, before each neutral refill. This aligns the LRC "collision-backed handoff"
(refill installed before the last site is spent; THM-4335) with Session 309's fixed-bank obstruction.
