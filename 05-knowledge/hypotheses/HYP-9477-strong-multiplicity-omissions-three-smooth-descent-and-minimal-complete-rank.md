---
id: HYP-9477
title: "Definability lenses, conjectural part. (A) The multiplicity function m(h) of strong counts never takes the sixteen values 9, 13, 21, 22, 34, 37, 41, 44, 50, 51, 55, 57, 60, 62, 64, 67 omitted below the observed floor 73, and {3, 5, 11, 13} is the complete set of counts with a unique strong class (Carmichael exceptional set). (B) Every count h > 63 is attained by at least two non-isomorphic strong tournaments that are not converses of each other. (C) The counts of crossing chains (the first-in-first-out sub-language of G(0,2), rank 14) are exactly the values of an explicit sum of 3-smooth terms indexed by the gaps, and their value set is described by a greedy 3-smooth descent in the style of Blecksmith-McCallum-Selfridge and THM-4719's ladder windows. (D) The minimal Hankel rank of a tournament grammar whose range is cofinite in the odd numbers is at least 10 and at most 17."
status: (A) FINITE-EXACT below 625 (THM-4884 (2)), OPEN beyond; (B) FINITE-EXACT below 625 (m(h) >= 6 for 63 < h < 625, and converse pairing), OPEN beyond; (C) OPEN, the closed form not yet derived; (D) OPEN (upper bound 17 is HYP-9470 (C); lower bound 10 is the failure of the rank-10 nested chains, not a theorem)
session: opus-2026-11-05-S81
date: 2026-11-05
related: [THM-4884, THM-4878, THM-4882, THM-4744, THM-4719, THM-4748, HYP-9470, HYP-9475]
script: 04-computation/experiments/spine_grammar_sublanguages_20261105.cpp; 04-computation/experiments/spectrum_polynomial_congruence_20261105.cpp; 04-computation/experiments/definability_lenses_prereg_20261105.md
---

# HYP-9477 - multiplicity omissions, three-smooth descent, minimal complete rank

## (A) Omitted multiplicities and the Carmichael exceptional set

THM-4884 (2): below 625 the multiplicity function omits exactly the values 9, 13, 21, 22, 34, 37, 41, 44, 50, 51, 55,
57, 60, 62, 64, 67, 79 among 0..80 (every other value up to 80 occurs; the maximum is 327 at h = 615), and
m(h) >= 73 on [375, 623]. Conjectures (separate assertions): m(h) grows without bound along the odd numbers above 623; the
sixteen omitted values below the observed floor 73 are never multiplicities of a strong count; and
the unique-class counts for h >= 3 are exactly 3, 5, 11, 13. The strong inventories at orders 10-12
show m(h) of order 10^2 at h near 600. Growth without bound alone would not prove the omission or
singleton claims: finitely many larger target counts could still realize any omitted small value.
This implication is explicitly withdrawn by the Session 218 scope correction (2026-11-07). No claim is
made about 79 or any omitted value >= 73, which larger h may still attain. (The first draft singled out the six values
up to 38; the audit observed that the same argument covers every omitted value below the floor.) Test: extend THM-4744's ear inventories to H < 1,000 (orders <= 14) and read m(h).
Session 200's open source-fiber conjecture (every n in F(S) other than 3 and 11 has at least two physical profiles)
is the profile-level analogue; THM-4884 (2b) shows the exceptional sets differ (5 collides on one class, 13 is not a
C1 readout at all), so neither conjecture implies the other.

## (B) Two realisations above 63

Conjecture: every odd h > 63 other than those excluded by Chan-Pak is the count of at least two non-isomorphic strong
tournaments which are not converses of each other (a strengthening of Chan-Pak with multiplicity). Evidence: m(h) >= 6
for every attained h in (63, 625). A proof mechanism: any complete family with an involution-free pairing of words
(e.g. interval chains with a non-palindromic composition and its reversal give converse pairs; two different
compositions with the same multiset of parts give the same product and non-isomorphic, non-converse chains when the
chains are distinguished by their ear structure, THM-4745).

## (C) Three-smooth descent for crossing chains

Mode 2 of THM-4884 (3) (crossing chains: each vertex opens at most one reversed arc, terminates at most one, and the
older pending arc always terminates first) has rank 14 and attains 998 of the 1024 odd numbers below 2048. Its count
is a product of 14-dimensional transfer matrices indexed by the gaps between events; conjecture: it has a closed form
as 1 + a sum of 3-smooth terms 2^a 3^b over a structured index set (THM-4748's crossing formula is the two-arc case),
and its value set is the set of odd numbers admitting a greedy 3-smooth representation of that shape. The test is to
derive the closed form from the automaton (the eigenvalues of the gap matrices are 2^g and 3^g) and to compare the
value set with the greedy rule on [1, 2047]; a match would reduce the completeness of crossing chains, and with the
shared-endpoint corrections of THM-4884 (3) the completeness of cutwidth 2, to a representation theorem.

## (D) Minimal complete rank

The one-apex series (rank 4) attains half of the odd numbers, the nested chains (rank 10) attain 440 of 1024 below
2048, the two rank-14 sub-languages miss 10 and 26 values, G(0,2) (rank 17) misses only 7 and 21 to 8191. Conjecture:
the minimal Hankel rank of a tournament grammar (a grade G(k, w) or a sub-language of one) whose range contains every
odd number above a bound is between 10 and 17. (The note's "between 5 and 17" means only that no complete candidate
of rank <= 4 is known, the rank-4 one-apex series attaining half of the odd numbers; the floor 10 is this conjecture,
from the failure of the rank-10 nested chains, not a theorem.) Lower bounds would need a structural argument (a rank-r series has
growth rates among r eigenvalues; a complete range needs two radices and the shared-endpoint pieces).
