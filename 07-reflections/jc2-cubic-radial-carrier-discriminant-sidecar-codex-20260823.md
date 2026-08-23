# Cubic radial carriers: one polynomial controls both gates

**Status:** independent proof rederivation plus **FINITE-EXACT** controls for
proved THM-3771.  This sidecar agrees with the canonical theorem and adds one
compact representation.  It does not prove JC(2), which remains **OPEN**.

For

    z=XT,       U=Xu(z),       W=U+3z+r,       Q=U phi(W),

in the genuine quadratic-carrier sector deg(phi)=2, define

    V(S)=(S-r)u((S-r)/3)phi(S).                       (1)

The roots of V are exactly the values of W on the irreducible components of
Q=0: the axis value r, the profile values r+3 beta, and the carrier roots.
Thus V restores precisely the component label lost by the target coordinate
Q.

The independent derivation gives three equivalent smoothness tests:

    Q is smooth;
    V is squarefree;
    u(0)phi(r)!=0, u and phi are squarefree,
    and gcd(u(z),phi(3z+r))=1.                        (2)

For d=deg(u)>=1, this has the coefficient-level certificate

    Disc(V)=3^(-d(d+3)) Disc(u) Disc(phi) u(0)^2 phi(r)^2
            Res_z(u(z),phi(3z+r))^2.                 (3)

For constant u=u0, the boundary formula is
Disc(V)=u0^4 Disc(phi)phi(r)^2.

The log-canonical chart has

    J(U,W)=3U,       k(X,T)=k(U,W)=k(Q,W),
    P0=-cW/(3Q),     J(P0,Q)=c,

and every rational mate is P0+H(Q).  On the component with W-address s, its
possible simple-pole coefficient is h_(-1)-cs/3.  Simultaneous regularity is
therefore exactly

    V(S) divides h_(-1)-cS/3.                         (4)

In the quadratic sector deg(V)=d+3>=3, whereas the right side is nonzero
linear.  This independently recovers the vertical no-entry mechanism and
also shows that a smooth zero fibre has exactly d+3 irreducible components.

The exact companion checks the universal chart identities and rational mate;
the discriminant formula on profile degrees zero, one, and two; unit gradient
ideals for those three positive controls; and explicit critical points for
u(0)=0, phi(r)=0, repeated roots of either profile, and a profile/carrier
collision.  Normal, optimized, and frozen executions must byte-match.

This packaging connects THM-3771 to THM-3770 without enlarging either scope:
the address set in THM-3771 is Roots(V), while (4) is the simple-pole
specialization of THM-3770's vertical equalizer.  The viable Jacobian frontier
therefore lies beyond radial product dressing: a nonbirational log-canonical
map with a synchronized component spectrum remains the exact open target.
