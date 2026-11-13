---
id: THM-5156
title: "Selmer-refined spectrum below 125. Let nu(T) = dim_F2 ker(M^2 + M + I), the Monsky nullity, where M is the mod-2 adjacency matrix with diagonal chosen so every row sums to 0. (1) nu is additive over ordinal joins, so the nu = 0 tournaments are exactly the ordinal sums of strong nu = 0 tournaments, and their Hamiltonian-path counts are exactly the multiplicative closure of the strong nu = 0 counts. (2) With Busch's bound, the nu = 0 spectrum below f(10) = 125 is decided by strong classes of order <= 9. (3) FINITE-EXACT: the odd numbers below 125 that are the Hamiltonian-path count of no nu = 0 tournament, at any order, are exactly 3, 7, 13, 19, 21. In Monsky's dictionary, where a tournament on distinct primes = 3 (mod 8) has nu equal to the pure 2-Selmer dimension of y^2 = x^3 - N^2 x, the counts 3, 13, 19 force nontrivial pure 2-Selmer, just as 7 and 21 are impossible outright"
status: >
  PROVED: (1) and (2).
  - Additivity of nu over ordinal joins is KNOWN (05-knowledge/results/bsd_component_recursion_20261013c.md) and was
    re-verified on all 1,007 non-strong classes of orders 2-8.
  - Busch's bound is THM-4744.
  - Unlike the chart class of THM-5152, the nu = 0 class is closed under ordinal sums, so the closure is exact in
    both directions.
  FINITE-EXACT: (3), over every class of orders 1-9 (classgen lists). The nu-distributions by order are in the
  packet README.
  CREDITED:
  - the nu-to-Selmer dictionary is Monsky's formula as recorded in 05-knowledge/results/integer32_twoedge_20261113
    (bridges/RESULT.md B2), whose census also defines nu;
  - this session's sweep pre-checked 3, 13, 19 through order 6 with sampling at order 7;
  - the decision to 125 is new here.
  NOT CLAIMED: values from 125 up (127 is undecided); Chan–Pak (OPEN); BSD (OPEN).
session: kind-pasteur-2026-11-13-S30
source: 05-knowledge/results/kps_s30_observer_lattice_20261113/README.md
proof: 05-knowledge/results/kps_s30_observer_lattice_20261113/README.md
related:
  - THM-5152
  - THM-4744
  - THM-1370
  - 05-knowledge/results/integer32_twoedge_20261113/README.md
  - 05-knowledge/results/bsd_component_recursion_20261013c.md
---

# THM-5156 - the Monsky nullity zero sector attains every odd count below 125 except 3, 7, 13, 19, 21

The proof and data are in the packet README, section "Selmer-refined spectrum". Replay: `scripts/monsky.py` in a
directory with the classgen lists cls_n.txt, n = 1..9.
