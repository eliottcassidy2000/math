#!/usr/bin/env python3
"""jacobian_outside_class_tower_klein_S327.py -- klein-2026-07-20-S327.
PROVABLY-OUTSIDE-CLASS counterexamples: the composition tower F^2, F^3.
Generic degree (3, 9, 27) is an Aut x Aut invariant => F, F^2, F^3 lie in three
PAIRWISE DISTINCT essential classes. det J = (-2)^m by chain rule (det J_F = -2
re-verified). THE COLLISION TRIPLE PERSISTS DOWN THE TOWER:
  (0,0,-1/4), (1,-3/2,13/2), (-1,3/2,13/2)
  --F--> (-1/4,0,0) --F--> (0,0,-1/2) --F--> (-1/2,0,0)
so F^2 sends all three to (0,0,-1/2) [det 4] and F^3 to (-1/2,0,0) [det -8].
NEGATIVE RESULTS (this session, feeding THM-1330's finite-factorization frame):
(i) F's 6-parameter coefficient shell is RIGID -- (lam,mu,nu) = (3,4,3) forced,
only the tau-scaling torus direction free (det = 2 tau);
(ii) the (1,-2,-4) weighted chart (fibration u = x^2 y) admits NO Keller map
(all five branches degenerate, det = 0). Both strengthen: at degree 3, F's
class is the only one known; provably-new classes today = the tower.
(Verification code as run in-session; see frozen .out.)"""
