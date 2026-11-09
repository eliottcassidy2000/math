---
id: HYP-9565
title: "Escapers of generation grammars and the density of first creators: (A) the decomposability defect of escapers (non-creators with no first-creator deletion, prime and 3-arc-strong) is unbounded; (B) the share of first creators among the n-tournaments tends to 1; (C) k*(12) = 3, so no 4-arc-strong 12-tournament is a non-creator; (D) every tournament is an induced subtournament of a first creator"
status: "OPEN. Audit G1 PASS WITH CORRECTIONS (applied; MISTAKE-743). (A): dec = simplicity index (Müller-Pelant 1974; credited 2026-11-09, MISTAKE-746); certified escapers of index 4 (order 16) and 5 (order 18), and PROVED conditional on (S_c) for unbounded index (THM-4995); at least 3 at order 12 (THM-4985 (2)); the bands have defect 1, so they say nothing here. (B): 0.719, 0.793, 0.8545, 0.8884, 0.9200 at orders 7-11 (THM-4985 (6)); PROVED conditional on (M_eps). (C): evidence in section (C): f_4(12) <= 229,725 from the explicit C3[R5, R5, TT_2] (lambda = 4), more than twice max H(11) = 95,095; annealing found only 232,001. (D): PROVED conditional on (M_eps) (THM-4985 (6)). Pre-registered Q1-Q3 (commit fb31d2c38a): all three REFUTED as predictions (see the last section)."
session: mac-mini-2026-11-09-genesis
date: 2026-11-09
related: [THM-4985, THM-4986, HYP-9532, THM-4942, HYP-9481, THM-4922]
script: 04-computation/experiments/genesis_20261109/sa12.c; 04-computation/experiments/genesis_20261109/sa13.c; 04-computation/experiments/genesis_20261109/ext13.c; 04-computation/experiments/genesis_20261109/tau_share.c
---

# HYP-9565 - escapers, defects and the density of first creators

Packet: [note](../results/genesis_20261109/README.md), [PROOF](../results/genesis_20261109/PROOF.md). The proved and
finite-exact facts are THM-4985 and THM-4986.

## (A) Unbounded decomposability defect

An **escaper** is a tournament that is not a first creator, has no first-creator deletion, and is strong, prime and
3-arc-strong. These are the tournaments THM-4942's five operations do not reach. Their **decomposability defect** dec
is the fewest arc reversals to a non-strong or non-prime tournament.

- At order 12 the escapers found have dec in {1, 2, 3} (THM-4985 (2)).
- The band escapers have dec = 1, with lambda unbounded under (S-inf). They are THM-4880's ribbon tournaments with the
  constant letter {1..s} (credit after audit G1, M2).
- For n >= 3, dec equals the module defect, since a non-strong tournament is never prime (audit G1, m1).

Conjecture: for every d there are escapers with dec >= d. Then no grammar built from first creators by a bounded number
of reversals of decomposable tournaments works at all orders.

**Update 2026-11-09 (mac-mini-2026-11-09-paleyeven, THM-4995; audit P1, MISTAKE-746).**
- **Priority.** dec is the classical simplicity index s(T) (Müller-Pelant 1974; Boussaïri-Lakhlifi-Talbaoui 2022),
  which was not credited here before.
- **Unconditional records.** Certified escapers of simplicity index 4 (order 16, a double of a near-regular
  8-tournament) and 5 (order 18, the double of C_9{1,2,3,5}).
- **Conditional.** Doubles D(B) of doubly regular B have dec = lambda = (q + 1)/2 (PROVED). Under (S_c) they are
  escapers for all large q (q >= 55 when c = 1), so (A) holds conditionally.
- **The heuristic below was too optimistic at small orders.** The matched back arcs multiply H(A)H(B) by 132 (q = 7)
  and 1954 (q = 11).

The bounded-defect question is Session230's (primes/REPORT.md §5: "whether all prime noncreators admit a bounded-cut or
bounded-defect derivation from creator parents at arbitrary order"; SYNTHESIS §6). (A) makes it precise: escapers are
the prime non-creators that THM-4942's grammar misses, and dec measures the defect. Credit added after audit G1 (M3).

Heuristic: re-connected ordinal sums of quasi-random halves.
- Take A ⊕ B with A and B quasi-random m-tournaments (for instance Paley), re-connected by r back arcs on distinct
  vertices. Then lambda = r, and every vertex pair has about m/2 distinguishers.
- H is at most H(A)·H(B)·m^(O(r)), which is far below max H(2m - k - 2). So under (S-inf) with g(m) = max H(m - 1),
  e.g. HYP-9481 (a), these are non-creators of creator depth > k. Bare (S-inf) is not enough, because H(A)·H(B) is
  superexponential (audit G1, m6).
- The two defect-3 escapers of order 12 have exactly this shape (m = 6, r = 3): six vertices of score 8 send only 3
  arcs back to six vertices of score 3.
- What is missing is a proof that no module larger than a pair is cheap.

NEGATIVE, small cases: random-pattern ribbons of span s <= 6 at n = 14 have lambda <= 2 and dec <= 2, because their
end vertices have few neighbours (`ribbon_defects.out`). Ribbons are no evidence for (A) at small sizes.

## (B) First creators have density 1

The share of first creators among isomorphism classes is 0.719, 0.793, 0.8545, 0.8884 and 0.9200 at orders 7-11.

Conjecture: the share tends to 1, both for classes and for labelled tournaments.
- The labelled version follows from (M_eps) (THM-4985 (6)), which in turn follows from max H(m) = O(m!/2^(m - 1)).
- The class version also needs the classical fact that almost all tournaments are asymmetric, so that class and
  labelled shares have the same limit (audit G1, m7).

## (C) k*(12) = 3

k*(n) is the largest lambda of a non-creator of order n (HYP-9532 (D)); k*(12) >= 3 by THM-4985 (1).

Conjecture: k*(12) = 3. It is equivalent to: every 4-arc-strong 12-tournament has H not in S(11).

Evidence:
- The substitution C3[R5, R5, TT_2] has lambda = 4 and H = 229,725, more than twice max H(11) = 95,095. It is the
  lambda = 4 counterpart of S17's lambda = 3 floor C3[R5, R5, 1] = 11,475, but nothing here shows that f_4(12) equals
  229,725 (audit G1, m5).
- Annealing for the least H with minimum semidegree 4 finds 232,001 at order 12 (`minH_f4_12.out`). The same code reproduces S17's exact f_3 at orders 8-10 but returns 15,213 > 11,475 at order 11, so it is an upper bound only.
- At order 13 the same annealing gives f_4(13) <= 926,401, above 531,205, the largest known count of order 12 (THM-329;
  max H(12) is not known exactly and is at most 1,141,140). This is order-13 data, not evidence for k*(12).
- Certified lower bounds: k*(15) >= 4 and k*(16) >= 4 (band_4(15) and band_4(16) are non-creators with lambda = 4; THM-4985 (3)).

## (D) Down-universality

Conjecture: every tournament is an induced subtournament of some first creator. This is PROVED under (M_eps)
(THM-4985 (6)).

## Pre-registered outcomes (Q1-Q3)

- **Q1 REFUTED** (predicted YES, 60%). Two runs of `sa13.c` (seeds 101 and 202, 400,000 steps each) found no order-13 escaper with lambda = 4 (`sa13_seed101.out`, `sa13_seed202.out` are empty). This is a failed search, not a non-existence proof. Not pre-registered but consistent:
  - among one-vertex extensions of all 193 order-12 escapers with minimum semidegree 4, every one has a deletion with
    H >= 171,039 > max H(11) = 95,095, so none is an escaper (`ext13.c`, `ext13_diag.c`). Those minimising the largest
    deletion have H about 7.8·10^5, above 531,205, the largest known count of order 12 (THM-329). This item was moved
    here from (C) after audit G1 (m5): order-13 data say nothing about k*(12);
  - `minH` gives f_4(13) <= 926,401.
- **Q2 REFUTED** (predicted YES, 45%). No order-13 escaper was found at all.
- **Q3 REFUTED** (predicted a decrease, 65%). The labelled share of new spanning-tree counts is 0.4190, 0.5150, 0.5719 at m = 6, 7, 8 (`tau_share_le8.out`).
