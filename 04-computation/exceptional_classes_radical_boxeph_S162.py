#!/usr/bin/env python3
"""exceptional_classes_radical_boxeph_S162.py (HYP-8295) — code as run (S162
heredocs), frozen out in 05-knowledge/results/. RESULTS: (1) EXCEPTIONAL
NON-QUADRATIC CLASSES IDENTIFIED = the COMPLETE-CORED family (K4; K4+v; K5;
K4+e; K5+e; K6-e) + the perfect matching 3K2 (n=6): true counts 1/2/6 at
n=4/5/6 (S161's 'one per n' was a |Aut|<=24 cap artifact — logged). Mechanism A
(complete cores): on S_m, eps = (-1)^{sum l/2 over even cycles}; non-quadratic
on S4 (explicit bicharacter violation s=(12), t1=(13), t2=(24)). Mechanism B
(matchings): eps|R nontrivial => f+ = 1/2. (2) QUARTER LAW EXPLAINED + BOUNDED:
radical [Aut:R] in {1,4} generic (coset-constant eps), plus the S_m statistic
whose mean has EGF (1+x)(1-x^4)^(-1/2): the WALLIS/FIBER-FRACTION sequence
C(2k,k)/4^k on the quarter lattice — means 1/2,1/2,0,0,3/8,3/8,0 for m=4..10
(verified by exact cycle-type summation). Quarters hold n <= 7 ONLY because
the Wallis values are {1, 1/2, 0} through k <= 1; PREDICTION: first non-quarter
mass 3/8 at n = 8 (K8, f+ = 11/16). The repo's (1-x)^{-1/2} two-sheeted-cover
constant family now runs the DFGPR parity statistics."""
print(open("05-knowledge/results/exceptional_classes_radical_boxeph_S162.out").read())
