---
id: THM-5153
title: "Exact slot pairings of merged tournament classes, with self-converse structure, through order 9 (FINITE-EXACT). Call two merged classes (tournaments up to converse) exactly paired when their full slot polynomials Q agree; such pairs are invisible to every uniform exterior. For n = 3..9 the merged classes number 2, 3, 10, 34, 272, 3528, 97144 = (A000568 + A002785)/2 and fall into 2, 3, 8, 24, 156, 1728, 50511 Q-classes. The Q-classes holding at least two merged classes number 0, 0, 2, 9, 69, 827, 21217. Of these, 0, 0, 1, 3, 15, 34, 158 hold two or more self-converse classes, and 0, 0, 1, 2, 28, 69, 1175 mix a self-converse class with a converse pair, so Q does not determine self-converseness from order 5. Self-converse classes alone in their Q-class number 2, 2, 5, 5, 37, 53, 1348. The odd values below 75 attained by no self-converse class of order <= 9 are exactly 7, 19, 21, 23, 35, 39, 43, 47, 53, 55, 59, 63, which is Session 315's global list. Among the seven strong members of strong/decomposable exact pairs through order 9, the only self-converse one is the rigid order-8 class with H = 45"
status: >
  FINITE-EXACT: every class of orders 3-9 (classgen of first_creators_genealogy_20261107, with self-converse flag and
  converse canonical code; fprof profiles of kps_s26).
  KNOWN, credited:
  - the Q-class counts were recorded by the kps-S26 audit (audit/AUDIT.md), and 24 at order 6 by Session 261
    MODULE_REBASE;
  - A000568 and A002785 are classical;
  - the self-converse gap list below 75 and its global proof are Session 315's; this is an independent finite check
    of it.
  NEW: the merged / self-converse refinement of the exact-pairing census, and the self-converse status of the
  THM-5151 (4) pairs.
  NOT CLAIMED: any order above 9; a recurrence or generating function for Q-classes; Chan–Pak (OPEN).
session: kind-pasteur-2026-11-13-S28
source: 05-knowledge/results/kps_s28_exact_pairing_sc_census_20261113/README.md
proof: 05-knowledge/results/kps_s28_exact_pairing_sc_census_20261113/README.md
related:
  - THM-5151
  - THM-5152
  - 05-knowledge/results/session315_dyadic_pairing_20261113/README.md
  - 05-knowledge/results/session314_metagraph_constraint_algebra_20261113/README.md
  - 05-knowledge/results/session310_positive_clocks_20261113/NEUTRAL_TAILS.md
---

# THM-5153 - exact slot pairings of merged classes, with self-converse structure, to order 9

The census table, its definitions and the two cross-checks are in the packet README. The computation joins two
inputs:
- classgen's canonical class lists, which carry for each class H, strongness, the self-converse flag, |Aut| and the
  canonical code of the converse;
- the fprof full slot profiles.

Each merged class is checked to have equal profiles for T and T^op, as Q(T^op) = Q(T) requires. Replay:
`scpairs.py` in the packet.
