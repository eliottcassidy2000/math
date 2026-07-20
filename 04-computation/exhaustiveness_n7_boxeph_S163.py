#!/usr/bin/env python3
"""exhaustiveness_n7_boxeph_S163.py (HYP-8310) — code as run (S163 heredoc),
frozen out in 05-knowledge/results/exhaustiveness_n7_boxeph_S163.out.
RESULTS: full n=7 census (augmentation cover 156 x 64 -> all 1044 classes):
(1) EXHAUSTIVENESS CONFIRMED: 22 non-quadratic classes, ALL explained —
20 by MECHANISM A (induced K4 with S4 <= Aut; now a THEOREM: the S4
bicharacter-violation triple (s=(12), t1=(13), t2=(24)) embeds, and subgroup
failure is global failure — covers K7 itself, |Aut| = 5040), 2 by MECHANISM B
(eps|R nontrivial: the 3K2+v matching family). No third mechanism exists at
n <= 7. (2) QUARTER LAW HOLDS at n=7 (Wallis truncation predicts break only
at n=8): masses {-1/2: 8, 0: 548, 1/2: 32, 1: 456}. (3) BONUS: mass-1 count =
456 = A000568(7) — the DFGPR equinumerosity re-verified at n=7 from the mass
spectrum. (4) THE 2^2 READING (owner): the quarter lattice is x^{2^2} — two
nested squarings (even cycles x2, mod-4 diameter sign x2); the break needs
k=2 in C(2k,k)/4^k = two coexisting 4-blocks, first possible at m = 2*2^2 = 8."""
print(open("05-knowledge/results/exhaustiveness_n7_boxeph_S163.out").read())
