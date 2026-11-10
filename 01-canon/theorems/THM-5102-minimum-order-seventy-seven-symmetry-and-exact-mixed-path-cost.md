---
id: THM-5102
title: "Minimum-order strong order77 symmetry and its exact mixed path cost"
status: "PROVED minimum-order classification and independently audited FINITE-EXACT sharp floor on that stratum. All-order77 path minimum OPEN."
session: Session240, 2026-11-10
source: 05-knowledge/results/session240_seven_eleven_symmetry_20261110/mixed/PROOF.md
audit: 05-knowledge/results/session240_seven_eleven_symmetry_20261110/audit_gap/AUDIT.md
---

# THM-5102: the first strong joint carrier needs19 vertices

[Proof](../../05-knowledge/results/session240_seven_eleven_symmetry_20261110/mixed/PROOF.md)
· [Independent audit and direct19-vertex DP](../../05-knowledge/results/session240_seven_eleven_symmetry_20261110/audit_gap/AUDIT.md).

A strong tournament admitting an actual automorphism of order77 has at least
19 vertices. At order19 its sharp minimum Hamiltonian-path count is

    3049061631 = 77 * 39598203.

At that order the action has profile(7)(11)(1). Coprime orbit lengths force
uniform cross arcs; strongness forces the quotient to be C3. Thus the complete
universe is C3[A,B,1], with seven- and eleven-circulants A,B. Below19 only the
profile(7)(11) can carry order77, and its quotient is a non-strong pair.

For a=H(A), x=F2(A), b=H(B), y=F2(B), the inherited THM-3121
(`path-cover-walk-content-substitution-kernel`) gives

    H(C3[A,B,1]) = 3ab+xb+ay+xy.

The exact core census has simultaneous minima a175,x1638,b92411,y1571526.
There are6*10 factor choices and two quotient orientations, hence exactly120
minimizing labelled tournaments for a specified(7)(11)(1) permutation.
A concrete witness uses connection sets{1,2,3} mod7 and{2,3,4,5,10} mod11.
The complete512 strong invariant orientations are counted by the kernel;
a separate vertex-level C++ DP checks the witness without that kernel.

This is not an all-order77 path floor. It requires an element of order77,
not merely77 dividing the automorphism-group order or H. No isomorphism-class
count, full automorphism classification, or Chan-Pak completeness is asserted.
