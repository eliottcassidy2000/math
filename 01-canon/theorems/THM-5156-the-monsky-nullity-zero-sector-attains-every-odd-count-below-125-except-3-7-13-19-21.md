---
id: THM-5156
title: "Selmer-refined spectrum below 125. Let nu(T) = dim_F2 ker(M^2 + M + I), the Monsky nullity, where M is the mod-2 adjacency matrix with diagonal chosen so every row sums to 0. (1) nu is additive over ordinal joins, so the nu = 0 tournaments are exactly the ordinal sums of strong nu = 0 tournaments, and their Hamiltonian-path counts are exactly the multiplicative closure of the strong nu = 0 counts. (2) With Busch's bound, the nu = 0 spectrum below f(10) = 125 is decided by strong classes of order <= 9. (3) FINITE-EXACT: the odd numbers below 125 that are the Hamiltonian-path count of no nu = 0 tournament, at any order, are exactly 3, 7, 13, 19, 21. The audit strengthens this: every odd value from 125 to 2655 is attained by an ordinal sum of strong nu = 0 classes of order <= 9, so the same five are the only gaps below 2657, the first undecided value (a prime whose 14 strong order-9 classes all have nu = 2). Read through Monsky's formula on the Legendre tournament of distinct primes = 3 (mod 8), where nu is the pure 2-Selmer dimension of y^2 = x^3 - N^2 x: H in {3, 13, 19} forces nu = 2 exactly (Selmer, not rank), and H in {1, 5, 11, 23, 35, 43, 47, 63} forces nu = 0"
status: >
  PROVED: (1) and (2).
  - Additivity of nu over ordinal joins is KNOWN (05-knowledge/results/bsd_component_recursion_20261013c.md). It
    was re-verified on all 1,007 non-strong classes of orders 2-8. The audit proves it as an exact block identity
    and checks all 14,410 non-strong classes of orders 2-9.
  - H is multiplicative over strong components (THM-3166).
  - Busch's bound is THM-4744.
  - Unlike the chart class of THM-5152, the nu = 0 class is closed under ordinal sums, so the closure is exact in
    both directions.
  FINITE-EXACT: (3), over every class of orders 1-9 (classgen lists). The nu-distributions by order are in the
  packet README.
  CREDITED:
  - the nu-to-Selmer reading is Monsky's formula on the Legendre tournament, first recorded in-repo in
    05-knowledge/results/bsd_selmer_growth_20261013a.md (sections 2-3) and restated in integer32_twoedge_20261113
    (bridges/RESULT.md B1);
  - opus-S90 (HYP-9657 (C)2) found that the H = 13 order-5 class has positive nullity;
  - this session's census sweep (scratch, uncommitted) pre-checked 3, 13, 19 through order 6;
  - the decision to 125 and, by the audit, to 2657 is new here.
  AUDIT: PASS WITH CORRECTIONS (audit/AUDIT_THM5156.md; MISTAKE-773). Claims (1)-(5) reproduced. The false side
  remark "127 is undecided" is corrected: 127 is attained by six strong order-8 classes.
  NOT CLAIMED: values from 2657 up; any rank or Sha statement beyond what the nu-to-Selmer reading gives; Chan–Pak
  (OPEN); BSD (OPEN).
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
