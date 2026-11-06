# A one-arc repair of 775, with its exact failure boundary

Status: VERIFIED finite witness and PROVED exact family identity, independently audited. No eventual threshold or all-u extraction is claimed. No novelty claim is made for the inherited transitive-clone expansion.

In E(129), with names r=0,c=1,p_j=j+2, reverse c->p_0. Add x with outgoing mask11 on the old names. Then add y with outgoing mask83 on the old names plus x. The resulting rows are

    [2,0,4091,3059,2019,4035,1923,3843,3587,3075,2059,83].

This actual twelve-vertex tournament has775=6*129+1 Hamilton paths. It is nonstrong. The old-arc reversal turns E(129) into a transitive tournament; the later births are not a formal subtraction of sectors. This is an upper bound of one on the required old-arc reversals for this particular physical recipe, not a proof that zero reversals are globally impossible for this target.

`verify775.py` constructs all edges and counts all spanning directed paths by an independent forward subset DP. C++ discovery instead uses ending/starting path profiles and exposed two-path cuts. Normal and optimized Python runs agree.

## Do not extrapolate the masks

Apply the same named reversal and birth masks to E(2^k+1), k>=4. The exact count is

    H(k)=26*3^(k-4)+10*2^(k-4)-7.

In particular H(6)=267 rather than391, H(7)=775, and H(8)=2259 rather than1543. It hits the target6*(2^k+1)+1 exactly when k=7.

To prove the identity, freeze the nine-vertex graph at k=4. Every additional p_j is a transitive clone with fixed outgoing mask387 on that core. For an ordering of the nine core vertices, let a be the number of gaps that can contain clones. A gap between two core vertices x,y is allowed precisely when x->clone->y; the first and last gaps have the corresponding one-sided tests. Let t be the number of core steps with the wrong directed orientation. Every such step must be replaced by a nonempty allowed gap. If any is not allowed, that ordering contributes zero. Otherwise assigning the m=k-4 named clone vertices to allowed gaps, each internally in increasing order, gives

    sum_(j=0)^t (-1)^j binom(t,j) (a-j)^m.

This is inclusion-exclusion for the mandatory gaps. Each resulting Hamilton word has a unique core order and clone assignment, so this counts every path once. At m=0 use0^0=1; the retained zero-base coefficient is0 anyway. Exact enumeration of all9! core orders yields coefficient vector

    [0,-7,10,26,0,0,0,0,0,0,0].

The finite core calculation proves the formula for every m, rather than interpolating it from sample clone lengths. `check_repair_polynomial.py` independently computes the gap coefficients and checks the literal DP controls.

For the unique-hit assertion put m=k-4. Equality with the target is equivalent to

    13*3^m-43*2^m=7.

The left side at m=0,1,2 is -30,-47,-55, at m=3 it is7, and it is strictly increasing for m>=2 because its next difference is26*3^m-43*2^m>0. Thus only m=3 works.
