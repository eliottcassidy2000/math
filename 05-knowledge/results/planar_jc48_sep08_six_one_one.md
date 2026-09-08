# The complete quartic boundary type 6+1+1 has no polynomial mate

**Status: PROVED + FINITE-EXACT + INDEPENDENTLY AUDITED.**
For the all-finite class, a rational mate forces L constant. An actual
same-partition constant-L family has rational mates, so the polynomial
boundary is sharp. No bound on the unknown mate's degree is imposed.
JC(2) remains OPEN.

## 1. Statement, inheritance, and actual coordinates

Keep the actual DG source

    W=(P1_x x P1_z) minus S, S={z=x^2}, t=1/(z-x^2),
    omega=dx wedge dt,
    x=1/r, t=-r^2-r^4 b_D, omega=r^2 dr wedge db_D.

Let H in L2 and L in L1 be global, deg_t H=2, and F=H^2+L. Suppose
the complete binary-octic leading section of H has three distinct zeros
of multiplicities 6,1,1.

**Theorem.** If all three zeros are finite and a rational G in C(x,t)
satisfies J(F,G)=1, then L is constant. Consequently no polynomial mate
of any degree exists. The polynomial conclusion also holds when one
of the three roots is infinity.

The closest proved mechanisms are
[leading exactness](planar_jc48_sep08_leading_exactness.md),
[the complete radical classification](planar_jc48_sep08_boundary_exactness.md),
and the trace/value/first-jet argument in
[five_one_one_one](planar_jc48_sep08_five_one_one_one.md).
The sixfold M-unit exception is retained using
[quartic_common_root](planar_jc48_sep08_quartic_common_root.md) and
[the componentwise pole-degree gate](planar_jc48_sep08_pole_degree.md).
The corrected near miss is to copy the fivefold regularity assertion at
a balanced sixfold point, where second-kind poles really occur. The
least-used sidecars are the actual D zeros, the full shifted global
rows, and an ordinary field trace before any geometric normalization.

The five live concepts are global sections, local pole budgets, shifted
carriers, formal exactness, and field trace. The leading-section map
loses lower rows; this proof restores all of them. A rational chart may
omit an original source line, so no source-critical conclusion is
transported through that chart. The nonconstant-L obstruction below is
rational throughout. The final contradiction is an explicit Euclidean
sequence, not a mate-degree census.

## 2. Complete global rows and the first reductions

Write

    H=N(x)t^2+P(x)t+Q(x), L=M(x)t+R(x).

The proved [DG filtration](planar_jc48_sep08_dg_quadratic.md) gives

    deg N<=8,
    P=2N8 x^6+2N7 x^5+sum_(i=0)^4 p_i x^i,
    Q=N8 x^4+N7 x^3+(p4-N6)x^2+(p3-N5)x+q0,
    M=sum_(i=0)^4 m_i x^i, R=m4 x^2+m3 x+r0.          (1)

In particular M=0 forces L constant.

For the all-finite case write u=x-p at the sixfold point and

    N=alpha*u^6*(u^2-s*u+h).

Distinctness gives alpha*h!=0. Leading exactness requires
du/sqrt(N) exact, and the proved degree-eight classification gives

    3s^2=4h.                                         (2)

Thus s!=0. The actual global scaling x=lambda X,z=lambda^2 Z,
t=lambda^(-2)T with lambda=s/2 changes s to 2 and h to 3. Its Jacobian
is a nonzero constant. Dividing H by its new nonzero leading scalar and
F by its square removes alpha. These are valid source/target scalings
for both rational and polynomial mates. The sixfold point remains an
arbitrary finite p. Henceforth

    N=u^6 C, C=u^2-2u+3.                              (3)

Writing u=x-p does not assert that translation extends to W. Section 3
computes the actual shifted rows.

The universal inverse coefficient T2=-M/(4N), proved in Section 5, has
an exact differential in C(x)(sqrt(N)). Normalized trace to C(x) makes
-M dx/(4N) exact there. Its two simple-root residues force C|M.

Use the full local numerators, with s0=z-x^2:

    Ncal=N+s0 P+s0^2 Q, Mcal=M+s0 R,
    E=Ncal^2+s0^3 Mcal-f s0^4,
    eta=s0^2 du/E_(s0).                              (4)

The two simple boundary points are regular for eta on every normalized
branch, including shared points: the tangential derivative of Ncal is
a unit. Infinity is not a zero of N|S in this all-finite case.

Suppose M(p)!=0. The complete local sixfold M-unit analysis leaves only
normal order j=2 as a possible pole source. At a simple leading cubic
root eta is regular. At its possible double root the generic critical
value order is at most four, because its f-derivative has exact order
four. The total possible pole degree of a rational primitive there is
at most max(lambda-2,0)<=2. The remaining branch is regular. The
j=0,1 and j>=3 cases are regular. This includes the odd normalized
branch, the nonzero logarithmic case, and residue-free order-four
case; none is discarded by assuming generic smoothness.

But H|D is quadratic with leading coefficient N8=1, so F|D is quartic
with leading coefficient one. At a generic D point of a compact fibre,
the actual chart gives ord eta=2. A primitive would have local mapping
degree three there. Its poles can only be on S: eta is regular on
the smooth generic fibre in W, and any pole of a rational function
gives a pole of its derivative. The total boundary budget is at most
two. The component meeting this D point therefore cannot have degree
three. This componentwise contradiction uses no geometric
irreducibility hypothesis. Hence M(p)=0.

If M(p)=0 but P(p)!=0, the normal-unit shared-root calculation is
regular on all branches. Its centre s0=psi(u) has order six, and the
generic correction order is at most six. All other boundary points are
regular too, so eta is holomorphic on every relevant compact generic
component. It is nonzero on each component meeting W away from D,
and cannot have a rational primitive. A generic component cannot equal
the fixed D or the boundary S; possible exceptional denominator
components of G occur at only finitely many fibre values. Thus P(p)=0.

The proved [finite first-jet logarithmic obstruction](planar_jc48_sep08_shared_roots.md)
at this singular shared point now gives P'(p)=M'(p)=0. Its local
tangent calculation uses s0=uZ. Since [u^2]N=0, the tangent quartic is
Z^2[(P'(p)+Q(p)Z)^2+M'(p)Z+(R(p)-f)Z^2]. If P'(p)!=0 it has two
simple nonzero roots for generic f; if P'(p)=0 but M'(p)!=0 it has one.
Each actual normalized branch has nonzero residue Z0^2/P_f'(Z0).
Thus this necessity is rational and makes no assumption on the other
branches. Together with C|M and deg M<=4,

    M=kappa*u^2*C.                                   (5)

If kappa=0, the theorem's rational conclusion already follows. Assume
kappa!=0 for the remaining argument.

## 3. All shifted coefficients and the next local log

The complete remaining family has the compact rational expression

    q=1+u^2 t-2p/u,
    H=u^2 C q^2+(b u^2+c u+a)q+k+2pa/u,
    L=kappa C q+ell+6kappa p/u.                      (6)

The apparent poles cancel in the original source. Expanding gives

    P=2u^6-(4+4p)u^5+(6+8p+b)u^4+(c-12p)u^3+a u^2,
    Q=u^4-(2+4p)u^3+(3+8p+4p^2+b)u^2
        +(c-12p-2pb-8p^2)u+a+k-2pc+12p^2,
    R=kappa u^2-kappa(2+2p)u+kappa(3+4p)+ell.         (7)

After x=u+p these satisfy every global row in (1). Conversely (1)
and P(p)=P'(p)=0 give precisely (7). The coefficient map from
(a,c,b,k) to (P2,P3,P4,Q0) in u has determinant one, and ell is free.
Thus this is the entire space, not a translated subset.

If a=P2!=0, set s0=u^4 Z in (4). Write

    n(u)=N/u^6, p2(u)=P/u^2,
    B(u,Z)=n(u)+p2(u)Z+u^2 Q(u)Z^2.

Here B(0,Z)=n6+aZ with n6=3. Its derivative in Z is a unit near
Z0=-n6/a. Thus it has a unique moving analytic zero Zbar(u), with
Zbar(0)=Z0. This retains the first and all higher numerator units;
in particular the zero need not stay equal to Z0. Equation (4) gives

    E/u^12=B(u,Z)^2+u^2 Z^3[M/u^2+u^2 Z R]-f u^4 Z^4.

Put Z=Zbar(u)+uY. Since M/u^2=m2+O(u), m2=3kappa, the constant
equation after division by u^2 is

    a^2 Y^2+m2 Z0^3=0.

Both roots Y0 are nonzero and simple. They give two actual branches
with u as normalized parameter. On each,

    eta=[Z0^2/(2a^2Y0)] du/u + regular terms.         (8)

The residue is nonzero; higher coefficients and the generic fibre
value do not cancel it. Thus a=0 is necessary. This local exclusion
is rational and works at every finite p.

## 4. The actual quadratic field and its trace

Now a=0. Set w=H and Z=f-ell-w^2. Eliminating q from (6) gives

    E2=A(w)u^2+B(w)u+C0(w)=0,

    A=Z^2+kappa b Z-kappa^2(w-k),
    B=kappa(c-12p)Z+kappa^2[2(w-k)-6pb],
    C0=kappa^2[-3(w-k)+36p^2-6pc].                   (9)

Equivalently E2=(uZ-6kappa p)^2+
kappa(bu+c)(uZ-6kappa p)-kappa^2 C(w-k).
The inverse recovers q=(Z-6kappa p/u)/(kappa C).

The map to (F,H) is dominant: the q^2 coefficient of
J_(u,q)(L,H) is kappa*u*C(2u-6), not zero. Moreover the leading f^2
coefficient of the discriminant in (9) is

    kappa^2[(c-12p)^2+12(w-k)-144p^2+24pc],

a nonconstant linear polynomial in w. It is not square in C(w).
A square polynomial in f over C(w) must have square leading
coefficient, so this discriminant is not square in C(w)(f).
Thus C(x,t) is genuinely quadratic over C(F,H). This field argument
does not assume geometric generic irreducibility.

Since omega=du wedge dq/u^2, the actual derivative identity gives, on
E2=0,

    eta=-kappa dw/(u^2 E2_u).                        (10)

Off the fibre, the exact identity has the additional term
C'(u)E2/(kappa C(u)); the source retains and checks it. Put
Y=2Au+B, so Y^2=B^2-4AC0. Then

    1/u^2=(B^2-2AC0+B Y)/(2C0^2).

The normalized trace of (10) is therefore

    -kappa B dw/(2C0^2).                            (11)

A rational mate makes this the derivative of its normalized trace
in C(f)(w). Trace commutes with the uniquely extended derivation
holding f fixed. The linear polynomial C0 has root
w0=k+12p^2-2pc. Vanishing of the rational residue in (11) requires
B'(w0)=0, giving

    (c-12p)(k+12p^2-2pc)=kappa.                     (12)

In particular c-12p!=0. The conclusion includes B(w0)=0: cancellation
of a double pole does not remove its remaining residue condition.
No primitive on a chosen geometric branch is assumed to descend
without this actual field trace.

## 5. Universal inverse coefficients before specialization

Let A0=sqrt(N) and solve F(x,T)=v^(-4), with leading
T=A0^(-1)v^(-1), in K((v)), K=C(x)(A0). Every rational source
function embeds: a nonzero polynomial denominator has a unique lowest
Laurent term after substitution. For a rational mate,

    partial_x G(x,T)=(1/4)v^5 partial_v T.            (13)

Thus every positive-index T_j=[v^j]T used below has an exact
differential in K: the derivative of [v^(j+4)]G is (j/4)T_j.
There is no mate-degree bound.

Retain the unspecialized definitions D=P^2-4NQ and Z0=R-MP/(2N).
Center y=t+P/(2N). The whole polynomial becomes

    F=(Ny^2-D/(4N))^2+My+Z0.

Formal integration by parts and the change T'(v)dv=dt give

    T_j=-(1/j)[y^(-1)]F^(j/4), j>0.                 (14)

The minus sign converts the residue at t=infinity to its Laurent
coefficient. Centering leaves that residue unchanged. Expanding
all terms of orders y^(-2),y^(-3),y^(-4) yields

    T1=D/(8N^(3/2)),
    T2=-M/(4N),
    T3=-D^2/(128N^(5/2))-Z0/(4sqrt(N)),
    T5=D^3/(1024N^(7/2))+(DZ0-M^2)/(32N^(3/2)),
    T7=-5D^4/(32768N^(9/2))-3D^2Z0/(512N^(5/2))
            -3Z0^2/(32sqrt(N))-3DM^2/(256N^(5/2)).   (15)

These are universal identities, before any branch condition. The
source verifies them by the complete multinomial expansion in (14)
and, independently, by recursively solving every inverse coefficient
through T7 with independent A0,D,M,Z0. It then reconstructs each actual
family residue directly from these universal formulas, independently
of the collected expressions used for the elimination.

For T2, normalized trace K/C(x) gives the rational primitive of
-M dx/(4N) used in Section 2, including the degree-one square case.
At a simple root its residue is -M(q)/(4N'(q)), so M(q)=0.
For the odd rows, the two points above u=0 are unramified, with
sqrt(N)=+/-u^3 sqrt(C). Their ordinary Laurent residues must vanish.
The sign changes the odd residues by one common nonzero scalar.

In the a=0 family,

    D2=(bu+c)^2-4kC, D=u^6 D2,
    Z0=ell-kappa b/2-kappa(c-12p)/(2u).              (16)

Translate F by a constant to take ell=kappa b/2. This is an actual
target translation, preserves every mate, and changes neither H nor
(12). The constant part of Z0 is then zero.

The T1 residue is a nonzero scalar times

    b^2+2bc+c^2/3=0.                                (17)

Reducing the T3 residue modulo (17) gives

    c^3(6b+c)=12kappa(c-12p).                        (18)

Since kappa and c-12p are nonzero, c!=0 is now proved. The apparent
b=c=0 branch is excluded, not divided away.

## 6. Homogeneous ratios and an explicit Euclidean certificate

There is no geometric normalization c=1. Instead the necessary
equations are homogeneous for

    wt(b)=wt(c)=wt(p)=1, wt(k)=2, wt(kappa)=3.

The T1, trace, T3,T5,T7 equations have weights 2,3,4,6,8. The source
checks each identity before dividing. Because c!=0, take the
coefficient ratios b/c,p/c,k/c^2,kappa/c^3, divide each equation by
its power of c, and rename the ratios. Put d=1-12p. Equations
(12),(17),(18) become

    3b^2+6b+1=0,
    kappa=(6b+1)/(12d),
    k=kappa/d+(1-d^2)/12.                            (19)

The only parameter denominator is d, already nonzero by (12).
The factor 6b+1 is also nonzero, but is never divided out.

Put V=d^2. Substituting (19) in the exact T5,T7 residues and reducing
modulo 3b^2+6b+1 leaves, after nonzero numerical factors and powers
of d only,

    E5=(96b+16)V^2+(1296b+240)V+900b+165=0,
    E7=(96b+16)V^4+(1296b+240)V^3+(180b+33)V^2
                       +(14952b+2744)V+10098b+1853=0. (20)

The source clears the full rational denominators and verifies the
actual reduced numerator factors -d^2 and d^2. No additional parameter
factor is removed. The contradiction even allows V=0.

Literal resultants, with their orientations, are

    Res_b(E5,3b^2+6b+1)=3 R5(V),
    Res_b(E5,E7)=-6 R57(V),

    R5=256V^4+3072V^3-2208V^2-2880V+225,
    R57=192V^4-4320V^3+2656V^2+3252V-255.             (21)

Any common b makes both resultants zero, even if a specialized
leading coefficient vanishes. Only this necessary direction is used;
no converse root lifting or leading-coefficient division occurs.

Starting from R5,R57, their primitive Euclidean remainders over Q are

    26496V^3-17248V^2-21648V+1695,
    1875232V^2-4290180V+305385,
    379612V-27815,
    1.                                               (22)

Each is a nonzero rational multiple of the remainder of the preceding
two. Only fixed rational numbers are divided out in primitive
normalization. A separate standard-library Fraction implementation
replays this entire sequence. Hence the two quartics have no common
complex root, contradicting (21).

Every denominator branch was paid before this reduction: kappa=0
gives L constant; a!=0 gives (8); c-12p=0 contradicts (12); c=0
contradicts (18). Coefficient ratios are not claimed to be surface
automorphisms. This proves rational exclusion for every kappa!=0
at every finite p.

## 7. Polynomial conclusion and both infinity placements

In the all-finite case, a rational mate forces L constant. A polynomial
mate is then impossible, since

    J(H^2+L,G)=2H J(H,G)

cannot be a nonzero constant in the original polynomial ring. This
does not bound deg G or transport polynomiality through a rational map.

If the sixfold point is infinity, finite N has degree two and two
simple roots; dx/sqrt(N) has nonzero residues at infinity. If a simple
point is infinity, finite N has degree seven and partition 6+1:
its finite primitive pole budget is four, but the unique infinity
point forces local degree five. Thus neither is exact. These are
already proved cases of leading exactness and the radical
classification. They exclude rational mates without moving infinity
or changing the source volume. The full polynomial theorem follows.

## 8. Sharp rational control and exact reproduction

For any finite p let u=x-p, C=u^2-2u+3, and put

    v=u+u^3t-2p,
    H=Cv^2+k, F=H^2+ell,
    G_H=(u+1)/(12u^2v), G_F=G_H/(2H).                (23)

These are the a=b=c=0,kappa=0 members of the complete global family.
The whole second-chart substitution is polynomial. Literal source
differentiation gives J(H,G_H)=J(F,G_F)=1. The leading section u^6C
has exactly the required distinct finite roots. Thus the rational
exception is realized in the exact class, for every finite p.
It does not contradict the polynomial exclusion.

The standalone source imports no inherited mathematical implementation.
There is no coefficient census. It checks universal inverse formulas
by two independent paths, full shifted section completeness, the
local logarithmic branch, exact quadratic elimination and on-curve
volume transport, quadratic-field nonsquareness, residue and trace
identities, homogeneous ratio substitutions, full reduced factors,
both resultants, an independent Fraction Euclidean sequence, and
the same-class global rational control.

Run from the repository root:

    python3 04-computation/planar_jc48_sep08_six_one_one.py
    python3 -O 04-computation/planar_jc48_sep08_six_one_one.py

Every gate raises an explicit exception; none uses Python assert.
The normal and optimized producer replays are byte-identical to the
frozen output: **63 exact gates, 772 output bytes**.

* Source: `04-computation/planar_jc48_sep08_six_one_one.py`, 12,727 bytes,
  SHA256 `8448a49b7a52a2b3a35c13c8666425144781c76ed2343eccd21add012b96c196`.
* Output: `05-knowledge/results/planar_jc48_sep08_six_one_one.out`,
  SHA256 `ecd81171fe3f0dbcb353eecb240b877a0e00e6488d279971dd63efaafe1a907d`.
* Semantic trace:
  `8faee9ffd19ef9d77353e8d6f3156d7ebd2d1831ea948eea08bef2c7e397b7a9`.

The [independent full audit](planar_jc48_sep08_six_one_one_audit.md)
accepts every denominator branch, actual quadratic field and trace,
all universal residues and exact coprimality. Independent direct-formal
extraction and original-source checks agree; both63-gate replays
match the frozen output. The theorem is accepted into the proved graph.
