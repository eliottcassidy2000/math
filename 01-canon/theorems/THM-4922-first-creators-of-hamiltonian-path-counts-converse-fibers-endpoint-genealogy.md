---
id: THM-4922
title: "First creators of odd Hamiltonian-path counts (the classes of least order n(h) with H = h): they are closed under converse and have no one-vertex strong component (strong through order 11); by THM-4889's H(T) <= sum H(T - v), above n(h) max H(n(h) - 2) every first creator has a first-creator vertex deletion; at every odd order n >= 3 the self-converse classes split exactly evenly between H = 1 and H = 3 (mod 4); FINITE-EXACT through order 10: the 'one class or one converse pair' dichotomy fails from h = 9 on (69 of 1,520 values hold it to order 9, 55 of 6,143 at order 10), the only orphan first creators at orders 4-10 are the converse pair creating 63, and no non-first-creator of order 3-10 has every deletion a first creator"
status: >
  PROVED (elementary):
  - converse-closed fibers with no one-vertex strong component (Session224 section 3);
  - the genealogy above n(h) max H(n(h) - 2), from THM-4889's counting identity and padding;
  - the odd-order self-converse balance (4). Its proof was supplied by the independent audit (audit/AUDIT.md section 4;
    anti-involution orbit model and an equivariant internal-arc involution), building on Session29's converse
    congruences and OCF (THM-466).
  FINITE-EXACT, from two independent class generators (classgen.c / order10.c, and the audit's own), both reproducing
  A000568 and klein's census:
  - every first creator is strong through order 11 (in general this is HYP-9481 (b), OPEN);
  - the dichotomy counts;
  - the orphan pair of 63 at orders 4-10;
  - no hypo analogue at orders 3-10;
  - the k-distributions;
  - the value-colour table.
  Order 10 was pre-registered (9e07033ca4): Q1-Q3 CONFIRMED (blind); Q4 CONFIRMED but pre-determined by klein's
  order-10 fiber sizes, which the registration quoted; Q5 FAILED (forced by THM-791).
  KNOWN and credited:
  - Session224 (session224_first_occurrence_20261107, 2054513960, same owner prompt, pushed first): the refutation of
    the dichotomy at h = 9, the (s, p) birth ledger through order 8, and the strong-component structure of minimal
    realizations;
  - THM-4889 and klein-2026-11-06 section 4: the counting identity H(T) <= sum H(T - v);
  - klein's least_order_10.txt: the first-creator fiber sizes for n(h) <= 10;
  - THM-791 and Session29 converse congruences: the even-order black cell;
  - THM-344, THM-4748 and the Session24 multiplicities: the two first creators of 63 are the Fibonacci fan F_6 and its
    converse;
  - Rédei and Camion: no hypohamiltonian or hypotraceable tournament (as merged in THM-4870's source packet).
  New:
  - the (s, p) fibers at orders 9-10;
  - the genealogy data with the orphan pair of 63;
  - the no-hypo analogue;
  - the threshold corollary;
  - the value-colour table;
  - the odd-order balance theorem.
  NOT CLAIMED: genealogy, orphan or no-hypo statements beyond order 10 (HYP-9512). Chan-Pak OPEN.
  AUDIT: independent, PASS WITH CORRECTIONS (5 majors, 15 minors, one upgrade; MISTAKE-735).
session: kind-pasteur-2026-11-07-S16
date: 2026-11-07
source: 05-knowledge/results/first_creators_genealogy_20261107/PROOF.md
related:
  - HYP-9512
  - THM-791
  - THM-344
  - THM-4748
  - THM-466
  - THM-4870
  - THM-4889
  - THM-1370
  - HYP-9481
  - 05-knowledge/results/session224_first_occurrence_20261107/README.md
  - 05-knowledge/results/path_spectrum_two_angles_metagraph_20261106/README.md
script: 05-knowledge/results/first_creators_genealogy_20261107/classgen.c; 05-knowledge/results/first_creators_genealogy_20261107/order10.c; 05-knowledge/results/first_creators_genealogy_20261107/audit/aud_enum.c; 05-knowledge/results/first_creators_genealogy_20261107/audit/aud_sc_mod4.py
---

# THM-4922 - first creators of odd Hamiltonian-path counts

Proofs and data: [PROOF.md](../../05-knowledge/results/first_creators_genealogy_20261107/PROOF.md).
Audit: [AUDIT.md](../../05-knowledge/results/first_creators_genealogy_20261107/audit/AUDIT.md).

A **first creator** of an odd h is a class of order n(h) (the least order attaining h) with H = h.

1. **First-creator fibers.**
   - **(PROVED)** The first creators of h are closed under converse.
   - **(PROVED)** For h > 1, none has a one-vertex strong component, and each strong component is a first creator of
     its own count (Session224).
   - **(FINITE-EXACT, orders <= 11; OPEN in general, HYP-9481 (b))** Every first creator is strong.
2. **(PROVED)** H(T) <= sum over v of H(T - v) (THM-4889). Hence, if h > n(h) max H(n(h) - 2), every first creator of h
   has a vertex v such that T - v is a first creator at order n(h) - 1.
3. **(FINITE-EXACT, orders <= 10)**
   - **The dichotomy is false.** The smallest failure is h = 9: two self-converse classes at order 5, neither the
     converse of the other (Session224, PROVED). h = 15 (order 5) and h = 17 (order 6) fail the same way, and h = 23 is
     the smallest value with two converse pairs. "One class or one converse pair" holds for 69 of 1,520 values (orders
     <= 9) and 55 of 6,143 (order 10), mostly at small values and near each order's maximum.
   - **Genealogy.** At orders 4-10, every first creator except the converse pair of first creators of 63 (order 8) has a
     vertex whose deletion is a first creator. Among values h >= 5 with n(h) <= 10, 63 is the only one whose first
     creators are not grown from first creators. The 3-cycle, at order 3, is also an orphan: its deletions have H = 1.
   - **No hypo analogue.** No class of order 3-10 that is not a first creator has every vertex deletion a first creator.
     At order 2 the unique class trivially is one.
   - **Value colours.** A value is blue if some first creator is SC, and black otherwise. These are value colours,
     distinct from the metagraph's line colours. Values 3 (mod 4) with even least order are always black (THM-791).
     Otherwise 82% are blue in aggregate, 0.66-0.86 by order (orders 6-10).
4. **(PROVED; proof supplied by the audit)** At every odd order n >= 3, the self-converse classes split exactly evenly
   between H = 1 and H = 3 (mod 4). An equivariant involution reverses the internal arcs of the pairs of an
   anti-involution, selected by a path-count parity. The counts are 1/1, 4/4, 44/44 and 1376/1376 at n = 3, 5, 7, 9.
