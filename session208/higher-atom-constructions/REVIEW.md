# Independent higher-atom audit

PASS, 2026-11-06. Independently reconstructed all nine decimal-coded tournament witnesses. All-source reachability confirms strong connectivity. For all three order-nine witnesses, exhaustively enumerated all9! permutations and also used a separate recursive enumeration of directed paths. For orders10–12, used recursive path enumeration without memoization, independent of the forward and reverse subset-DP discovery/verification methods. Every claimed path count and prime-exponent pair matches, including151875 and253125 at order12.

The exact dominance-join product identities and additive vertex counts are correct:1875²*2025³=15^14 on45 vertices;10125^4*9375=15^17 on50;54675*46875=15^8 on22;151875*253125=15^9 on24. The last proves asymptotic rate at most8/3, not optimality. The order-nine1875 atom violates the proposed global extension of prices(16/11,24/11), while leaving the prior finite-dictionary theorem intact.

The disjoint adjacent-position swap action is free; each orbit contains at most one directed Hamilton path because each pair's internal orientation is forced. Thus H<=n!/2^floor(n/2) is valid. The resulting exact inverse threshold and weaker asymptotic lower bound are correct. The lower bound is sublinear and does not imply positive limiting linear efficiency. Search budgets are discovery only, explicitly nonexhaustive, and no global minimum or complete higher-order dictionary is claimed.

Independent checker/output: independent_check.py and independent.out, with sibling-relative input for standalone replay. No direct enumeration of24-vertex joined paths is necessary: the exact SCC product rule with independently verified factors certifies the count.
