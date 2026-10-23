# When finite angular guards really do certify planar rigidity

This is a bounded-degree positive repair, not a proof of the planar Jacobian conjecture.

Let F=(P,Q) be a polynomial map over C, equivariant for the cyclic action
(x,y)->(zeta_L x,zeta_L^(-1)y), with the same weights on its output. Suppose its total degree is at most D and L>D+1. Then F is equivariant for the full action (x,y)->(t x,t^(-1)y), t in C*.

Proof. A monomial x^a y^b has integer charge a-b in [-D,D]. Every monomial of P has charge congruent to1 modulo L, and every monomial of Q has charge congruent to-1. Its difference from the required charge has absolute value at most D+1<L. Thus the congruence is an equality. This proves continuous equivariance monomial by monomial, with no hidden cancellation between distinct monomials.

If additionally det JF is a nonzero constant, then F=(a x,b y) for nonzero constants a,b. Indeed continuous equivariance gives P=x A(xy), Q=y B(xy). Writing u=xy, the determinant is (u A(u)B(u))'. A nonzero constant derivative forces uAB=c u (its constant term is zero), so AB is constant and both A,B are nonzero constants.

The charge-lifting threshold is sharp. At degree D=L-1, the map

(x+y^(L-1),y)

is C_L-equivariant and has determinant1, but is not continuously equivariant. For L>=3 it is nonlinear. At L=2 this particular boundary example is linear but still not continuously equivariant, so it demonstrates failure of the charge lift, not a nonlinear degree-one map. These maps are polynomial automorphisms and are not Jacobian counterexamples.

In particular, simultaneous C8 and C9 equivariance is C72 equivariance: the charge congruences combine modulo72. A Keller map of degree at most70 with these EXACT symmetries is therefore diagonal linear. Degree71 already admits the nonlinear shear (x+y^71,y). The numbers8 and9 here specify actual group actions on polynomial variables, not an arbitrary residue guard on a tournament statistic.

This illustrates a legitimate compiler contract: a finite guard is exact when its modulus exceeds a proved integer-charge bandwidth. Without the degree bound or the equivariance hypothesis, no rigidity follows. Gaussian power tests raise the bandwidth with the power, so a single finite modulus does not cover all powers.
