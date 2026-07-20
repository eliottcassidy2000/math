#!/usr/bin/env python3
"""gmc2_nogo_elimination_boxeph_S166.py (HYP-8350) — code as run (S166 heredocs),
frozen out in 05-knowledge/results/gmc2_nogo_elimination_boxeph_S166.out.
THEOREM (computer-algebra proof): the radial family
  P = (1+Z)(W - (b+cZ)(mu + lam*ZW)),  Z = X+iY, N = 2 Gaussians,
contains NO GMC(2)-nullcone member: the moment equations E1..E5 are already
inconsistent over C^4. Proof: E1 forces b(mu+2lam) = 2; substituting and
eliminating c via Sylvester resultants gives H1, H2, H3; stripping the spurious
factors (powers of lam, mu, and D = mu+2lam) leaves homogeneous forms with NO
common ray (exact gcd checks); stripped loci die separately: lam = 0 via the
diagonal identity sum_m t^m [s^m](1+s)^m F(s) = (1-t)/(1-2t) F(t/(1-t)), which
forces F(u) = c(1-u) — impossible for the zero-free exponential; mu+2lam = 0
via E1 = 2; mu = 0 via the slice. (Leading-coefficient caveat of the resultant
lands inside the stripped loci — flag for the Lean-grade writeup.)
COROLLARY: the third Gaussian in the GMC(3) counterexample is ESSENTIAL for
this mechanism — the telescope (1/2)_k collapse cannot be simulated radially."""
print(open("05-knowledge/results/gmc2_nogo_elimination_boxeph_S166.out").read())
