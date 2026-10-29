# Bounded slopes but an infinite fourth intercept moment

**Status: PROVED.** The mean identities below are inherited; the integer-moment boundary and second-moment formulas refine them. No every-integer orbit claim.

## 1. Scope and inherited results

Use the cheaper (3,3) controller and normalized odd Haar source law from `CARRY_CLOCK_POLES.md`. Write the stopped affine map as n ↦ U n+V, with U=3^J/2^D and V=C/2^D. Different complete subtrees have the product source law after exact original-source pullback. This statement concerns the Haar law, not resampling a fixed integer.

The pinned `affine_tree_recurrences_20261028/MOMENTS.md`, sections 2–3, already proves EU=1 and EV=7/12 for this rule, including monotone limits of its nonnegative generating functions. Its section 5 proves a *different*, critical-parity slope-moment threshold; that result must not be transported to the cheaper rule.

If K=k is the number of branch vertices, J=3k+1 and D>=5k+1. Therefore

    0<U <= (3/2)(27/32)^k <= 3/2.                         (1)

All nonnegative slope moments exist. Set m_r=E U^r, r>=0, and

    B_r=3^r/(2^(3(r+1))-1),
    L_r=3^r/(2^(r+1)-1)-B_r.

The root decomposition gives

    m_r=L_r+B_r m_r^3.                                   (2)

It is the least nonnegative root, by summing complete trees of bounded height and increasing the height. At r=0,1 the least root is one. Boundedness (1) removes any issue of infinite slope moments.

## 2. The theorem: finite through order three, infinite from order four

For integer r>=1, put y_r=E V^r, possibly infinite. Then

    y_1,y_2,y_3 < infinity,
    y_r=infinity for every integer r>=4.                 (3)

The precise recursive criterion, provided all lower integer intercept moments exist, is

    y_r<infinity iff alpha_r=(1+m_r+m_r^2)/7 <1.          (4)

For a branch root let b=2^(-a), and let (U_i,V_i), i=1,2,3, be independent stopped child maps in chronological order. Then

    U=3b U_1 U_2 U_3,
    V=b U_1 U_2 U_3+V_1 U_2 U_3+V_2 U_3+V_3.            (5)

For a leaf, U=3b,V=b. The root branch mass is p=1/7. In the multinomial expansion of V^r, the three pure child-intercept terms contribute exactly alpha_r y_r: their coefficients average to p m_r^2, p m_r and p. Every remaining term either contains no child intercept or has every child-intercept exponent strictly less than r. Since U_i<=3/2 and b<=1/8 on a branch, these terms have finite expectations whenever all lower y_j are finite. Correlations within (U_i,V_i) do not matter for this bound; one does not factor U_i from V_i.

To justify sufficiency without presuming y_r finite, set V^(h)=V on completed trees of height at most h, and zero otherwise; define U^(h) similarly. These increase to V,U almost surely because the controller's Galton–Watson law is subcritical. Apply the same expansion with the indicator that all children complete within height h-1. Pure terms are bounded by alpha_r E[(V^(h-1))^r], because their accompanying truncated mass and slope moments are bounded by the full ones. Every cross term is bounded uniformly in h by the finite lower moments. Hence

    E[(V^(h))^r] <= alpha_r E[(V^(h-1))^r]+constant.

When alpha_r<1 this is uniformly bounded, and monotone convergence proves finiteness. The base r=1 follows by the identical argument with only bounded root forcing (or the inherited mean formula).

For necessity, if y_r were finite then nonnegativity and the pure terms in (5) would give y_r>=alpha_r y_r+c_r, where c_r is the strictly positive leaf contribution. This is impossible for alpha_r>=1. Thus the equality case also diverges.

For r>=1, B_r<=B_1=1/21. Consequently f_r(x)=L_r+B_r x^3-x is strictly decreasing on [0,2], since f_r'(x)<=12/21-1<0. It follows that m_r<2 precisely when f_r(2)<0, and m_r>2 when f_r(2)>0. The test is

    L_r+8B_r = 3^r/(2^(r+1)-1)+7*3^r/(2^(3(r+1))-1).

For r=1 it is below 2. For r=2 it equals 9/7+63/511<2; for r=3 it equals 27/15+189/4095<2. For r=4, even L_4=81/31-81/32767 exceeds 2. This proves the first three finite moments and the divergent fourth. For r>4, divergence follows from V^r>=V^4 on V>=1, since the contribution of V<1 to EV^4 is at most one.

The fixed-K normalized moments in `CARRY_CLOCK_POLES.md` are all finite, but mixing over unbounded tree sizes creates this fourth-moment divergence. This is distinct from raw integer-carry moments, already infinite at each fixed K>=1.

## 3. Exact second and mixed moments

Let m=m_2 be the least nonnegative solution

    m = 648/511 + (9/511)m^3,

and put S=1+m+m^2, y=7/12, h=E(UV), t=E(V^2). Then

    h = (m/3)/(1-S/21),
    t = [m/9+(2/63)h S+(2/7)y h(m+2)]/(1-S/7).           (6)

For h, multiplying the two expressions in (5) gives the pure mixed-child contribution (S/21)h. Root-only and leaf terms sum to m/3 by (2). For t, the pure square terms give (S/7)t, the root-child cross terms give (2/63)hS, and the three child-child cross terms give (2/7)yh(m+2). Root-only and leaf squares sum to m/9. Finiteness proved above authorizes these expectation equalities; bounded U also gives h<infinity.

The symbolic checker independently expands (5) and integrates each monomial, retaining U_i,V_i correlation via h. It also certifies the threshold comparisons with rational arithmetic and gives numerical values as secondary approximations.
