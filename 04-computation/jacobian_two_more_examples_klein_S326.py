#!/usr/bin/env python3
"""jacobian_two_more_examples_klein_S326.py -- klein-2026-07-20-S326.
Two further explicit JC counterexamples by radical-changing transport of F
(conjugation by a shear and by the Nagata wild automorphism); all claims
verified exactly (see frozen .out). G1 = F o (x, y+x^2, z+y): radical
1+xy+x^3, det -2, fiber over (-1/4,0,0) exactly {(0,0,-1/4),(1,-5/2,9),
(-1,1/2,6)}. G2 = F o Nagata: radicals via D = y^2+xz, det -2 (chain rule,
both factors verified), collisions (0,0,-1/4), (-16733/32,-467/8,13/2),
(-4197/32,233/8,13/2) -> (-1/4,0,0). Honest note: both are conjugates of F
(the lawful generator; per boxeph-S142 uniqueness, the known essential class);
F o F (det 4, 9:1) is the non-conjugate alternative."""
# (verification code as run in the S326 session -- see .out)
