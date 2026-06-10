# THM-464: the b-ary subgrid games — classical Ramsey numbers as the height-1 row, and where the 2-adic seam actually lives

**Status:** parts A–B PROVED (one-liners, machine-verified); part C computed
(complete verifier, cross-validated); part D = the corrected seam statement
(pointing at THM-465 C). Scripts: `erdos592_ternary_seam_macmini_s3.py`.
**Source:** mac-mini-2026-06-10-S1 (HYP-2373)

## A. The games

Q_b(n,t): does a triangle-free graph on the leaves [t]^n of the t-ary height-n
tree exist with NO independent b-ary subgrid (choose b children at every node;
b^n leaves)? R_b(n) := least t where the answer is no. The THM-453 D1 König
bridge holds verbatim for every fixed b (full grids contain b-ary subgrids):
Q_b(n,t) SAT for all t ⟹ ω^n ↛ (ω^n,3)²; positive relations force R_b(n) < ∞.

## B. The height-1 calibration row (PROVED)

At n = 1 the b-ary subgrids of [t] are ALL b-subsets, so a witness is exactly a
triangle-free graph with independence number ≤ b−1, and the largest t carrying
one is R(3,b) − 1 by definition of the Ramsey number. Hence

> **R_b(1) = R(3, b)** — the classical Ramsey numbers R(3,2)=3, R(3,3)=6,
> R(3,4)=9, R(3,5)=14, … are the height-1 row of the (branching × height) table
> whose height-2 binary entry is R_2(2) = 5 (THM-459).

Machine-verified at b=3: SAT at t=5 (the 5-edge witness is C₅, the Ramsey-
critical graph, found by the solver), UNSAT at t=6.

## C. The grading matrix at n=2 (computed)

Binary game (b=2): free = invariant = v₂-graded = **v₃-graded** = 5. The v₃
control (this session) joins S2's results: at n=2 EVERY tested gap-grading has
the same cutoff — coarse gradings are free of charge, and nothing is special
about 2 at this level (the column space is flat; only signs matter there).
Ternary game (b=3): free SAT through t = 6 at minimum (sweep beyond t=6 running
at write-up; witnesses 3, 12, 31, 70 edges at t=3..6). R_3(2) > 6 > R_2(2):
killing rarer, larger subgrids is a weaker demand, so the cutoff climbs with b.

## D. Where the seam actually is (corrects HYP-2373's first guess)

The seam-follows-branching prediction FAILS at n=2 in the grading form (all
gradings equal) but HOLDS at n=3 in the algebra form (THM-465 C): the bi-dyadic
feature algebra (sign + v₂ per coordinate) carries witnesses through t=6 while
the equal-sized triadic algebra (sign + v₃) is feature-UNSAT already at t=4 with
zero CEGAR rounds. The 2-adic structure is special exactly where the BINARY
subgrid demands meet a column space with depth — i.e., the seam follows the
branching through the algebra needed to win, not through which quotient of the
row-gaps you grade by. Sharp open note: identify the algebraic property of the
v₂ class structure (the forced parity closure odd+odd→even = the ℤ→ℤ/2 quotient)
that the v₃ classes lack with respect to the Schur/composition constraints.
