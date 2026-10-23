# Independent Coulomb-powering audit

Status: PASS, major errors0. Covers COULOMB_POWERING_ADDENDUM.md and powering_checks.py/json, separately from the earlier frozen Coulomb report. No Git mutation.

The descended polynomial P_m(c)=13−12T_m((13−c)/12) satisfies C(t^m)=P_m(C(t)) and composition P_m∘P_n=P_mn. The second-kind factor is exactly U_(m−1), without an omitted factor2. Indeed Δ(C(t))=36(t−t^−1)^2, so the stated discriminant identity follows. As a polynomial identity it extends to t±1.

For c rational and below1, t>1 and U_(m−1)>0 for every positive integer m. Thus the multiplying square is rational and nonzero; both discriminants define the same nonzero square class. Since t^m>1 the positive mirror branch remains available. The c13,m2 counterexample correctly shows why the unrestricted square-class assertion fails when the factor vanishes. m0 is also correctly excluded.

For t=φ², the even-index Fibonacci/Lucas identities produce c_m=13−6L_(2m) and Δ(c_m)=180F_(2m)^2. All constants, initial terms and affine recurrence c_(m+1)=3c_m−c_(m−1)−13 are correct. F_(2m)>0 for m≥1, so square class5 never collapses. This proves quadratic irrational reciprocal parameters and transcendental mirror bases for the whole family, via the already audited Gelfond–Schneider application. Negative integer c_m excludes rationality only of the simple Virasoro algebra, as stated.

Recovered session32 CAT_MAP_PATH_FORESTS.md from the Git object. Its integral conjugacy uses U=[[1,−1],[0,1]], with AU=UB. Rooted-forest determinant det(I+L(Path_m))=F_(2m) applies for m≥1; m0's formal recurrence value is not the empty-graph count. The new note does not extend the graph formula to m0. Trace and forest count remain distinct observables. At m4 they are47 and21, giving c−269 and discriminant180·21², as stated.

Primary source spot-checks: Chebotarev–Shamis https://arxiv.org/abs/math/0602575 explicitly identifies det(I+L) with total spanning rooted forests; Chebotarev https://arxiv.org/abs/math/0612792 gives the path/Fibonacci framework. The inherited graph theorem and the displayed conjugacy supply the bridge; no tournament Hamiltonian or Monster conclusion is inferred.

The four integer c<1 with rational parameter, {0,−2,−7,−24}, follow exhaustively from same-parity factor pairs of144. An integer that is a rational square is an integer square, justifying the factorization argument. Parameters3/2,2,3,6 match. The smaller rational VOA L(c_(2,5),0)⊗L(c_(3,5),0) has charge−5 and is regular by the verified Dong–Li–Mason results cited in the addendum.

Independent executable powering_independent_checks.py uses SymPy exact polynomial expansion, rather than the author's rational-grid tests:12 discriminant polynomial identities,25 composition identities, and20 direct path Laplacian determinants. It independently verifies the integral conjugacy and golden trace/forest equations. These are finite independent controls; the all-m proof was checked algebraically above.
