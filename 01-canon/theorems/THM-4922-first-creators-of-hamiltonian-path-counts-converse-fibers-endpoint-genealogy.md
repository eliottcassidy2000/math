---
id: THM-4922
title: "First creators of odd Hamiltonian-path counts: the classes of least order n(h) with H = h are strong and closed under converse (s SC classes + p converse pairs); H(T) <= sum over v of H(T - v) (endpoint injection), so if h > n(h) max H(n(h) - 2) every first creator of h has a vertex deletion that is itself a first creator; FINITE-EXACT through order 10: the 'one class or one converse pair' dichotomy holds for only 69 of 1,520 values (orders <= 9) and 55 of 6,143 (order 10), the only orphan first creators (no first-creator deletion) are the converse pair of first creators of 63 at order 8, and no non-first-creator has every deletion a first creator (the tournament analogue of 'no hypohamiltonian tournament')"
status: >
  PROVED (elementary): strong and converse-closed fibers; the endpoint injection; the genealogy above
  n(h) max H(n(h) - 2). FINITE-EXACT, from an independent class generator (classgen.c and order10.c, which reproduce
  A000568 and klein's level-resolved census exactly at orders 3-9):
  - the dichotomy counts;
  - the orphan set {the two first creators of 63};
  - the absence of the hypo analogue;
  - the k-distributions;
  - the split by parity of n(h) and h mod 4.
  Order 10 was pre-registered (Q1-Q5, 9e07033ca4): Q1-Q4 CONFIRMED, Q5 FAILED, explained by THM-791.
  KNOWN and credited:
  - the least-order census, n(h) and N(n) in Ss(n) (klein-2026-11-06, path_spectrum_two_angles_metagraph_20261106);
  - THM-791 (SC tournaments of even order have H = 1 mod 4, forcing the even-order cell of the split);
  - THM-344 (H = 63 at order 8 in exactly two classes);
  - THM-4870 (Thomassen merge: no hypohamiltonian or hypotraceable tournament; non-critical vertices);
  - THM-1370 (the insertion lemma and 7, 21).
  NOT CLAIMED: the genealogy, orphan or no-hypo statements beyond order 10 (HYP-9512); Chan-Pak OPEN. AUDIT: pending.
session: kind-pasteur-2026-11-07-S16
date: 2026-11-07
source: 05-knowledge/results/first_creators_genealogy_20261107/PROOF.md
related:
  - HYP-9512
  - THM-791
  - THM-344
  - THM-4870
  - THM-1370
  - HYP-9481
  - 05-knowledge/results/path_spectrum_two_angles_metagraph_20261106/README.md
script: 05-knowledge/results/first_creators_genealogy_20261107/classgen.c; 05-knowledge/results/first_creators_genealogy_20261107/order10.c; 05-knowledge/results/first_creators_genealogy_20261107/analyze.py; 05-knowledge/results/first_creators_genealogy_20261107/hypo.py
---

# THM-4922 - first creators of odd Hamiltonian-path counts

Proofs and data: [PROOF.md](../../05-knowledge/results/first_creators_genealogy_20261107/PROOF.md).

A **first creator** of an odd h is a class of order n(h) (the least order attaining h) with H = h.

1. **(PROVED)** Every first creator is strong, and the first creators of h are closed under converse.
2. **(PROVED)** H(T) <= sum over v of H(T - v). Hence, if h > n(h) max H(n(h) - 2), every first creator of h has a vertex
   v such that T - v is a first creator at order n(h) - 1.
3. **(FINITE-EXACT, orders <= 10)**
   - **The dichotomy is false.** The owner's "one class or one converse pair" holds for 69 of 1,520 values (orders
     <= 9) and for 55 of 6,143 values (order 10). It holds only at small values and near each order's maximum. The
     smallest failure is h = 23: two converse pairs at order 6.
   - **Genealogy.** Every first creator except the converse pair of first creators of 63 (order 8) has a vertex whose
     deletion is a first creator. So 63 is the only value whose first creators are not grown from first creators.
   - **No hypo analogue.** No class that is not a first creator has every vertex deletion a first creator.
   - **The two kinds** (blue: some first creator is SC; black: none is). Values 3 (mod 4) with even least order are
     always black (THM-791). Otherwise about 82-84% of values are blue.
