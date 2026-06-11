# THM-469: the sum-free grading dichotomy — why the seam is 2-adic (parity closure, leading digits, and the Schur arity)

**Status:** STUB — CLAIMED by kind-pasteur-2026-06-11-S1. Answers THM-464 D's sharp
open note / THM-465 C's caveat. Proof + computations this session; do not build on
this file until status changes.

## Claimed content

For a prime p, the level sets L_v = {g ≥ 1 : v_p(g) = v} of the p-adic valuation
are sum-free (and difference-stable-free) for ALL v ⟺ the unit group (Z/p)^× is
trivial ⟺ p = 2. ("Parity closure": odd+odd=even is forced; for p ≥ 3 the unit
escape 1+1=2 keeps g₁+g₂ in the same class.) Consequence for the feature games:
the (sign, v₂) algebra's classes are canonical Schur-free building blocks (the
dyadic rungs of THM-446), the (sign, v_p) classes for odd p each contain their own
Schur triples at every scale, strictly enlarging the realizable mono-triangle set
of the quotient — the algebraic property THM-464 D asked for.

Two falsifiable companion predictions (this session's experiments):
1. The asymmetry follows the SCHUR ARITY (triangles = 2-term gap sums), NOT the
   subgrid branching: in the b=3 (ternary-subgrid) game at n=3, (sign,v₂) should
   still beat (sign,v₃). [HYP-2390]
2. The 3-adic algebra is rescued by the leading digit: classes p^v·(u + pZ_p) are
   sum-free for every p, so (sign, v₃, leading-digit) should carry per-size
   witnesses where bare (sign, v₃) is feature-UNSAT. [HYP-2391]

Script: `04-computation/erdos592_parity_closure_kp0611.py` (+ .out).
