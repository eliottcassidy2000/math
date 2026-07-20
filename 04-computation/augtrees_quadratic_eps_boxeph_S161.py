#!/usr/bin/env python3
"""augtrees_quadratic_eps_boxeph_S161.py (HYP-8290) — code as run (S161 heredoc),
frozen out in 05-knowledge/results/. RESULTS: (A) both augmentation trees COVER
(evenness is cover-hereditary: every even (n+1)-class has an even parent — new
positive), but children-per-parent multisets diverge (4->5: [6,6,8,9] vs
[4,5,6,7]; 5->6 disjoint ranges): parent-matched recursion IMPOSSIBLE. (B) eps
is NOT quadratic — but fails on EXACTLY ONE class per n (1/11 at n=4, 1/32 at
n=5): localized exception. n=6 mass spectrum {-1/2: 2, 0: 88, 1/2: 10, 1: 56}:
NEGATIVE MASS exists; no 1/4 => Gauss mechanism refuted; true law = positive
fraction of eps quantized in QUARTERS {1/4, 1/2, 3/4, 1}. VERDICT: five
transport routes now closed (Aut S159, profiles S160, character S160,
quadratic/Gauss S161, parent-matched recursion S161) — the DFGPR equinumerosity
is maximally non-structural: cycle-index identity only."""
print(open("05-knowledge/results/augtrees_quadratic_eps_boxeph_S161.out").read())
