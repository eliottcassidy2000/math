# The all-finite and finite-seven placements of the (7,1) boundary class

**Status: PROVED + FINITE-EXACT + INDEPENDENTLY AUDITED.**
This concerns the complete global square-prefix class on the fixed DG
surface. The finite-simple/infinity-seven placement is not included.
JC(2) remains open.

## 1. Statement, inherited mechanisms and complete global entry

Keep

    W=(P1_x x P1_z) minus {z=x^2}, t=1/(z-x^2),
    x=1/r, t=-r^2-r^4b_D, omega=dx wedge dt=r^2 dr wedge db_D.

Let `H in L2`, `L in L1`, `F=H^2+L`, with `deg_t H=2`.
Suppose the leading binary octic has exactly two distinct roots,
of multiplicities seven and one.

**Theorem.** If the sevenfold point is finite, `F`
has no polynomial mate of any degree. This covers both placements
of the simple point and every original finite position. Actual
constant-`L` rational mates in each placement are retained below.
No assertion that rational mates always force constant `L` is made.

The closest suppliers are [the complete global filtration](planar_jc48_sep08_dg_quadratic.md),
[M-unit local forms](planar_jc48_sep08_quartic_common_root.md),
[the finite first-jet obstruction](planar_jc48_sep08_shared_roots.md),
and [leading inverse exactness](planar_jc48_sep08_leading_exactness.md).
The polynomial weight argument is reproduced in Section 2, so the
proof is self-contained at that step; the broader
[all-m weight gate](planar_jc48_sep08_polynomial_weight_gate.md)
records its separate generalization. No reserved assertion is a
necessary dependency of the argument reproduced here.

Live concepts are original polynomial support, actual rational field
maps, the retained volume, generic moving residues and degree at
infinity of the residue identity. The source-to-rational-curve map
preserves the field and relative differential but loses polynomial
regularity. All polynomial restrictions are therefore established
before that map; subsequent obstructions are genuinely rational.
The canonical hostiles are the two actual rational families in
Section 7. The corrected expansion in Section 4 retains a lower
term of `B0` which does not alter its leading-degree obstruction.

At a finite sevenfold point an M-unit has regular relative forms:
seven is never the balanced multiple `3j` of a positive integral
normal order. A shared normal unit is regular too. At the other
simple point the complete M-unit/shared analysis is regular,
including if that point is the original infinity, where the
canonical factor `r^2` only improves regularity. If the sevenfold
point were inactive, a rational primitive would be holomorphic
and constant on every compact normalized generic component,
contrary to the nonzero actual relative form on its source part.
If the shared values vanish but first jets do not, the finite
first-jet lemma gives nonzero logarithmic residues. Consequently
any polynomial mate necessarily enters

    P(0)=P'(0)=M(0)=M'(0)=0                         (1)

at the finite sevenfold point `u=x-p=0`. This componentwise
entry does not assume geometric generic integrality.

In the all-finite case, actual source/surface scalings normalize
the nonzero separation and the leading scalar to
`N=u^7(u-1)`. The complete original global rows are

    P=2u^6-(4p+2)u^5+a u^4+b u^3+c2 u^2+c1 u+c0,
    Q=u^4-(4p+1)u^3+(a+4p^2)u^2
                         +(b-2pa+4p^2)u+d.         (2)

If the simple point is infinity, a leading scalar normalization
gives `N=u^7`, and the complete rows are

    P=2u^5+a u^4+b u^3+c2 u^2+c1 u+c0,
    Q=u^3+a u^2+(b-2pa-4p^2)u+d.                   (3)

In both cases `H=N t^2+P t+Q`. For every lower section write

    M=sum_(i=0)^4 m_i u^i,
    R=m4 u^2+(m3-2p m4)u+e, L=M t+R.               (4)

These are full coefficient spaces: the three original quadratic
numerator rows `Q,P-2Qx^2,N-Px^2+Qx^4`, with `x=u+p`, all have
degree at most four. The complete fixed-leading section space
has dimension six; the six independent coordinates
`a,b,c2,c1,c0,d` in (2) or (3) give exactly that space. The
linear rows `R,M-Rx^2` similarly give (4). The source verifies
the full symbolic all-p rows used below. Writing `u=x-p` is
an ordinary polynomial coordinate change, not a claimed surface
translation. The actual second chart and its volume are retained.

## 2. A polynomial first-jet gate before rational coordinates

Under (1), and since `ord_u N>=6`, put `w=u^2t`. Then
`F` lies in the polynomial ring `C[u,w]` and has the expansion

    F=f0(w)+u f1(w)+O(u^2),
    f0=(P2 w+Q0)^2+M2 w+R0,
    f1=2(P2 w+Q0)(P3 w+Q1)+M3 w+R1.                (5)

Here the subscripts mean original coefficients in `u`. For an
original polynomial `G(u,t)`, its image is a finite Laurent sum
`sum u^j g_j(w)`. If its least exponent `l` is negative, the
nonzero polynomial `g_l` is divisible by `w`: each original
monomial has exponent `j=i-2k` and carries `w^k` with `i,k>=0`.
The actual bracket is `u^2 J_(u,w)`.

If `f0` is nonconstant, the first bracket term for `l<0` is
`-l f0' g_l u^(l+1)`, nonzero. Only `l=-1` could give a
constant bracket, but its coefficient cannot be one because
`g_l` is divisible by `w`. For `l>=0` every nonzero bracket
term has positive `u` order, including the vanishing leading
constant case `l=0`. Thus `f0` must be constant, forcing

    P2=M2=0.

Now suppose `f1` is nonconstant. For `l<0` the first bracket
coefficient is `f1 g_l'-l f1' g_l`, at order `u^(l+2)`.
If the respective positive degrees are `s,j`, its leading
coefficient has the nonzero factor `j-ls`, and its degree is
`s+j-1>=1`. It cannot vanish or equal a nonzero constant at
the only possible order-zero case `l=-2`. For `l>=0` all
terms again have positive order. Hence `f1` must be a constant
`beta`. If `beta=0`, both original partial derivatives of `F`
vanish on `u=0`, so no constant-Jacobian mate exists. Therefore

    P2=M2=0,
    M3=-2P3 Q0,
    beta=2Q0 Q1+R1!=0.                              (6)

This is an original polynomial obstruction. Rational functions
do not satisfy its divisibility constraint and are not excluded
by it. Later rational charts do not transport polynomiality.

For the all-finite case, T2 exactness at the finite simple root
also requires `M(1)=0`. Combining (1),(4),(6) gives

    P=2u^6-(4p+2)u^5+a u^4+b u^3,
    M=lambda u^3(u-1),
    R=lambda[u^2-(1+2p)u]+e,
    lambda=2bd.                                    (7)

If `lambda=0`, `L` is constant and the final polynomial factor
obstruction applies. Otherwise `b,d,lambda` are all nonzero.

For finite-seven/infinity-one, the entire survivor is

    P=2u^5+a u^4+b u^3,
    M=u^3(lambda u+mu),
    R=lambda u^2+(mu-2p lambda)u+e,
    mu=-2bd.                                       (8)

The parameters in (7),(8) have not been normalized by a source
translation or by assuming that `p` vanishes.

## 3. All-finite case: the actual rational generic fibre

Put

    v=u+u^3t-2p, z=uv, A=a-4p, k=d+2bp.

Direct expansion of the entire original rows gives

    H=z^2+A z+k+v(b-z),
    L=lambda(z-v)+e-2p lambda.                      (9)

Translate the target `F` by the last constant. The source field
is exactly `C(z,v)`: `u=z/v` and
`t=(v+2p-u)/u^3` are its rational inverse. Its actual volume is

    omega=(v^2/z^3) dz wedge dv.                    (10)

Set `h=H,f=F,D=h^2-f,S=A+b`. Then the complete fibre relation
can be solved without an algebraic extension:

    Q=bD+lambda(k-h), E=D-lambda S,
    z=Q/E,
    v=z+D/lambda=W/(lambda E),
    W=D^2-lambda A D+lambda^2(k-h).                 (11)

Conversely substitution into (9) gives precisely `H=h,F=f`.
The map is dominant: at fixed `h` its original derivative is
`lambda[k+b(A+b)-h]/(b-z)^2`, nonzero. Thus these are actual
field coordinates, not an unproved quotient of a generic
component. The denominators are nonzero rational functions,
and `lambda!=0` was paid in (7).

Substituting the actual volume (10), or differentiating both
maps in (11), gives the relative differential

    eta=-W^2 dh/(lambda^2 Q^3).                     (12)

A rational mate would make (12) rational-exact over `C(f)`.
No polynomial regularity is inferred from the inverse maps.

## 4. All-finite case: a nonzero generic residue at every parameter

Since `b!=0`, define the fixed degree-two polynomial

    P(h)=h^2-(lambda/b)h+lambda k/b,
    w=P(h)-f,
    B1=2(h-k)/b-A,
    B0=(h-k)^2/b^2-A(h-k)/b-(h-k).                  (13)

Then `Q=bw` and the entire numerator is

    W=w^2+lambda B1 w+lambda^2 B0.                  (14)

The term `-A(h-k)/b` in `B0` is essential for the exact identity;
an initial scout omitted it, and the independent auditor restored
it before the primary/source freeze. It is lower degree, so the
leading obstruction below is unchanged. The source tests (14)
without discarding that term.

The polynomial part of `W^2/w^3` is rational-exact. Up to one
common nonzero scalar, the remaining differential has terms

    (B1^2+2B0) dh/w
       +2lambda B1 B0 dh/w^2
       +lambda^2 B0^2 dh/w^3.                       (15)

At a generic root of `P(h)=f`, the local variable `w` is a
parameter. Put `D_P=(1/P')d/dh`. Its residue is

    R(h)=(B1^2+2B0)/P'
       +D_P[2lambda B1 B0/P']
       +(1/2)D_P^2[lambda^2 B0^2/P'].               (16)

Exactness on the generic fibre requires this rational function
to vanish identically. Indeed substituting `f=P(h)` at the
generic root makes `h` transcendental over the original constant
field; a nonzero rational function of `h` cannot vanish there.
Critical values of the fixed quadratic `P` affect only finitely
many special levels and do not affect this argument.

But at infinity in this rational `h` line,

    R(h)=3h/b^2+O(1).

The first term of (16) has that leading term, the second is
bounded and the third tends to zero. Since `b!=0`, the residue
identity is impossible. Thus every nonconstant-`L` polynomial
candidate surviving (7) has no rational mate at all.

## 5. Finite-seven/infinity-one: full quadratic trace

Using the original rows (3),(8), put

    v=u+u^3t-2p, A=a+4p, k=d+2bp,
    h0=bv+k.

After an actual constant target translation, the complete
family becomes

    H=u v(v+A)+h0,
    L=lambda u v+mu v+e.                            (17)

Here `e` is the translated arbitrary constant. The inverse
`u=(h-h0)/[v(v+A)]` gives the actual field `C(v,h)`; the
denominator is not identically zero even when `A=0`. Its volume
and generic fibre are

    omega=v^2(v+A)^2 dh wedge dv/(h-h0)^3,
    f=h^2+delta(h-h0)+mu v+e,
    delta=lambda/(v+A).                            (18)

The extension over `C(f,v)` is genuinely quadratic. Its
discriminant is `delta^2+4(f+delta h0-mu v-e)`, linear in the
transcendental `f` with nonzero coefficient four, hence nonsquare
in `C(v)(f)`. No geometric connectedness assertion is needed.

Assume first `b!=0`. Put

    P(v)=h0^2+mu v+e, w=P(v)-f,
    B=2h0+delta, C(v)=v^2(v+A)^2.

In `y=h-h0` the fibre equation is `y^2+B y+w=0`.
The full quadratic trace, without normalization, satisfies

    Tr[1/(y^3(2y+B))]=1/w^2-B^2/w^3.

Consequently a rational mate forces the full trace differential

    C(v)[1/w^2-B^2/w^3]dv                           (19)

to be exact. A normalized trace is one half of (19), with the
same conclusion. This identity is checked independently by the
literal two-by-two companion matrix; neither a branch nor a
factor of two is suppressed.

Since `b!=0`, `P` has degree two. Let

    L1=C B^2/P', L2=C/P', D_P=(1/P')d/dv.

Generic residues in (19) give
`(1/2)D_P^2 L1-D_P L2=0` as a rational identity in `v`, by
the same generic-root argument as Section 4. Integrating once
in `C(v)` and multiplying by `P'` gives

    L1'-2C=c_* P'                                  (20)

for a constant `c_*`. But `L1~2v^5`, so the left side is
`8v^4+O(v^3)`. The right side has degree one or is zero.
This is impossible. The argument includes `lambda=0`, every
`mu`, every `A`, and both possible roots of the actual quadratic.

## 6. The remaining b=0 case and the polynomial conclusion

If `b=0`, the original polynomial gate (8) makes `mu=0`.
If `lambda=0` too, `L` is constant. Otherwise `lambda!=0`,
and (17) has the actual rational field parametrization

    z=uv=(f-h^2)/lambda,
    v=(h-k)/z-A.

With the original volume, its relative differential is

    eta=(h-k-Az)^2 dh/(lambda z^6).                 (21)

At a generic level `f=rho^2`, `rho` is transcendental. After
removing the common nonzero factor `lambda^3`, the residue of
(21) at `h=rho` is

    -[80A^2 rho^4+140Ak lambda rho^2
                  +63k^2 lambda^2-7lambda^2 rho^2]
                         /(512rho^11).             (22)

The top coefficient forces `A=0`, after which the coefficient
of `rho^2` forces `lambda=0`, a contradiction. This is a generic
residue identity, not a test at selected fibre values.

All nonconstant-`L` polynomial candidates in both placements
are therefore excluded. For constant `L`, the original polynomial
identity `J(F,G)=2H J(H,G)` has a nonconstant factor and cannot
equal one. This finishes the polynomial theorem without a degree
bound on the mate and without invoking a rational-fibre
Jacobian-conjecture theorem.

## 7. Actual rational sharpness and remaining scope

For every original finite `p`, put `v=u+u^3t-2p`. The two
families

    H=u(u-1)v^2+d,
    G_H=-(8u^2+4u+3)/(15u^3v),                      (23)

and

    H=u v^2+d,
    G_H=1/(5u^3v)                                  (24)

have respectively the exact leading sections `u^7(u-1)` and
`u^7`. Each is an actual global `L2` section for every `p`,
and each satisfies `J(H,G_H)=1`. Hence `H^2+e` has rational
mate `G_H/(2H)`. Its poles are essential; these examples do
not contradict the polynomial theorem.

The finite-simple/infinity-seven placement has leading `N=u`
after its scalar normalization and remains a separate OPEN
problem here. The all-finite formal inverse hierarchy also has
a separate scope: exact coefficient rows alone do not prove
rational algebraization. The present proof uses original
polynomial support followed by actual generic-fibre residues.

## 8. Exact controls and audit boundary

The [source](../../04-computation/planar_jc48_sep08_seven_one_two_locations.py)
checks both entire original-position carriers, actual field inverses,
all volume factors, the complete numerator expansion, both generic
degree obstructions, the full quadratic trace, the exact exceptional
residue and both sharp families. It does not search over a bounded
degree of mates. Reproduce with

```sh
python3 -B 04-computation/planar_jc48_sep08_seven_one_two_locations.py
python3 -B -O 04-computation/planar_jc48_sep08_seven_one_two_locations.py
```

Both modes pass **46 gates** and reproduce the identical
[322-byte frozen output](planar_jc48_sep08_seven_one_two_locations.out).
The source has5,925 bytes, SHA256
`62aac3b8ce69b6781b0434c32bf40d92a21fd36d30249866285db25d746c8dd1`.
Output SHA256:
`ee0d52487c625eb31ea960a20b5fed7f515f9696da162e72464b65c9a68b0431`.
Semantic gate SHA256:
`d6888a5c8c4ad87d99944d0224fe9b9098458fb83c8dbd23af7bc1d54a9b6d66`.

The complete independent analytic/source audit remains pending.
This candidate is outside the proved dependency graph until accepted.


## Accepted independent audit

The [complete independent audit](planar_jc48_sep08_seven_one_two_locations_audit.md)
accepts the analytic proof, full coefficient scope and actual source
coordinates. Both normal and optimized replays reproduce all46
gates and the frozen output. The candidate is now in the proved graph.
