# Independent Session 57 geometry and literal-variant audit

Verdict: the substantive claims in the two supplied notes pass independent review. No bare K4,4 rational point, nondegenerate magic square, prize solution, novelty claim, or classification of a resolution of X44 follows from these calculations.

## Proof review

- The m=(r-1)(s-1) monic quadrics form a regular sequence: after each previous relation, the next interior variable is still an indeterminate over the existing quotient. The monic-polynomial argument works even if that quotient has zero divisors. The free boundary-coordinate presentation proves dimension and the complete-intersection degree. Adjunction yields omega=O(2m-rs), without requiring smoothness. Thus X44 is six-dimensional, degree 512, with dualizing sheaf O(2), not O(0).
- The second additive-base presentation is a different finite free presentation, of rank 2^(rs). Distinct nonassociate linear prime factors have independent square classes by valuations. The multiquadratic generic algebra is a field, and freeness injects the whole algebra into it. This proves geometric integrality in characteristic different from 2.
- The projective squaring morphism is also finite flat, of degree 2^(rs-1). On a target chart z_e != 0, normalize q_e=1 and adjoin all other rs-1 coordinate ratios by monic equations u_f^2=z_f/z_e. These give a free basis of squarefree monomials. Restricting to the additive linear base preserves this freeness. The projective cover degree 32768 for X44 is not its embedding degree 512.
- In the Jacobian C diag(q), dependence of rows is equivalent to a linear relation among additive-base coordinate forms supported on zero coordinates. The graphic incidence rank theorem identifies dependence with a zero-edge cycle. At an actual additive-base point, every connected zero block is complete bipartite, so a cycle forces a zero rectangle. Four zero entries in a rectangle impose exactly three independent linear conditions. Finite surjectivity then gives singular codimension exactly 3, except X22 where the projective zero locus is empty. The complete intersection is S2, and codimension-three singularities imply R1, hence normality. This argument does not calculate discrepancies.
- The arithmetic gauge argument is valid: its interval has width equal to the smallest squared root entry, strictly positive. Avoiding finitely many cross-part collisions is possible over Q. Multiplication by the square of a common denominator clears both vertices and roots. This proves exactly the stated open-point existence equivalence, not the existence of such points.
- The decorated base has affine dimension 3 and projective dimension 2. Repeated forms contribute seven independent relative root signs, so there are exactly 128 geometrically integral components over Q. Fixing one root for each distinct form identifies each component with the nine-form root surface in P8; choosing three linearly independent forms shows it is a six-quadric complete intersection. Its degree is 64 and dualizing sheaf O(3). The degree total is 128*64=8192=2^13, agreeing with the full thirteen-quadric pullback. Global sign does not halve 128: these are relative signs within repeated-form classes. The generic product algebra is reduced and the free module over Q[C,x,y] is torsion-free, excluding nilpotents globally. Different components meet only where some of those repeated-form roots vanish.

## Independent exact controls

The new audit.py uses explicit RuntimeError checks, not assert, so all 72,097 checks run under python -O. It independently checks incidence/forest ranks for every edge subset of K22, K23, K33 and K44, including all 65,536 K44 edge subsets. It exhausts the nonzero additive bases over F3 for X22 and X44 and F5 for X23 and X33, verifying cycle versus zero-rectangle equivalence. Finite checks support, but do not replace, the preceding proofs.

It also symbolically expands both MM^T=N^2 I and M^T M=N^2 I, verifies the quaternion substitution and column permutation, verifies nine distinct positive squared entries, and recomputes all eight line totals with exact integers. The totals are six copies of 16008001=4001^2, followed by 26532801=5151^2 and 8526400=2920^2. Therefore the witness solves the unequal-square-line-sums interpretation and is expressly not magic.

Compiling and rerunning a local copy of the supplied search.cpp at bound 120 reproduces (48,32,23,12), matching indices 1 and 3, and counts 233046 / 160271 / 1229. The program stops at this witness and is not an exhaustive search through bound 120. Its 64-bit arithmetic is safe for that run; no unbounded-input safety claim is made.

The equal-common-square-sum obstruction is valid. Summing the four lines through the center gives T=3c. If c=e^2 is nonzero, 3e^2 is not a square. If c=0, nonnegative squared entries and T=0 force the entire grid to be zero. Thus no nonzero integer-square magic grid has square common sum.

## Engineering recommendation

The original check_geometry.py and certificate.py use assert, so python -O suppresses their checks. Their ordinary-mode results are valid, but these files should either be hardened to explicit failures or documented as requiring assertions. The new independent audit has no such limitation. No original file or remote was changed.

## Primary-source spot checks, 25 October 2026

- Auel and Singer, arXiv:2609.09351v1, https://arxiv.org/html/2609.09351v1 , Section 1.1 and Lemma 1.1: six quadrics in P8, degree 64, 256 ordinary double points, dualizing O(3), and a general-type resolution. Section 1.3 displays the singular-point orbits. This supports the comparison for the decorated components only, not X44.
- O. M. Cain, arXiv:1908.03236v2, https://arxiv.org/html/1908.03236v2 , abstract and Sections 2–4: quartic/factorization formulation and Gaussian-integer specialization. Theorem 2.1 independently states the center/common-total relation. The source does not identify every bare K4,4 with the decorated surface.

Both versioned pages were actually retrieved. No literature-status or novelty inference beyond these checked claims is needed for the witness or the proofs.
