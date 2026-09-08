# Constant boundary values from two moving inverse indices on every W_m

**Status: PROVED + FINITE-EXACT + INDEPENDENTLY AUDITED.**
This is a necessary condition on the specified global square-prefix
family. It does not exclude all rational mates, does not claim a
polynomial classification on W_m, and makes no Jacobian-conjecture claim.

## 1. Statement and complete source carrier

For every integer `m>=2`, use the proved surface and filtration in
[the all-m DG theorem](planar_jc48_sep08_dg_genus.md):

    W_m=(P1_x x P1_z)\{z=x^m},  t=1/(z-x^m),
    x=1/r, t=-r^m-r^(2m)b, D={r=0},
    omega=dx wedge dt=r^(2m-2) dr wedge db.           (1)

Let `H in L2`, `L in L1` be arbitrary global functions, with

    H=N t^2+P t+Q, L=M t+R, F=H^2+L,
    deg N=4m-2.                                    (2)

No assumption is made on the positions, multiplicities or number of
finite roots of N, or on its being square. The complete numerator boxes
are precisely

    A=N-x^m P+x^(2m)Q, B=P-2x^mQ,
    deg A,deg B,deg Q<=2m,
    S=M-x^m R, deg R,deg S<=m.                      (3)

Equivalently `N=Qx^(2m)+Bx^m+A`, `P=2Qx^m+B`, and
`M=Rx^m+S`. All coefficients in these boxes are retained subject to
the exact degree in (2).

**Theorem.** If any rational `G in C(x,t)` satisfies
`J(F,G)=1`, then `H|D` and `L|D` are both constant. Equivalently,
`F|D` is constant. The same statement holds for any nonzero constant
Jacobian after scaling the mate. There is no bound on its poles or degree.

The proof uses the two inverse indices

    j_1=2m-3,  j_2=6m-7.                             (4)

For m=2 these are exactly T1 and T5 from the
[degree-six predecessor](planar_jc48_sep08_degree_six_constant_d.md).
The proof here derives the all-m statement directly from the proved
filtration and formal chain rule; it does not depend on the status of
that stronger all-root degree-six sidecar. For m>=3 the same two fixed
indices need not see the boundary slopes. The moving indices are paid
by an exact weight argument, not inferred from finite experiments.

The inherited mechanism is the
[formal leading-field exactness map](planar_jc48_sep08_leading_exactness.md).
The live concepts are full numerator boxes, infinity pole degree, the
index of an inverse coefficient, pure monomial weight, and rational
exceptions. The source is the complete rational Jacobian equation; the
target is one labelled residue in its coefficient field. The map retains
that equation coefficientwise but discards algebraization and all finite
residues. The latter loss is exhibited by both controls in Section 5.

One proposed transfer formula had the wrong volume exponent `m-2`.
Direct differentiation of (1) gives `2m-2`, recovering order two when
m=2. This was corrected before any statement was used. The actual volume
is retained in the map contract, although the proof below computes in
the original x-coordinate and needs no compact-fibre volume estimate.

## 2. Uniform degree bounds and the actual boundary restrictions

Normalize N to be monic. If its leading coefficient is n!=0, replace
`(H,L,F,G)` by `(H/n,L/n^2,F/n^2,n^2G)`. This preserves globality,
the mate equation, the exact degree, and constancy on D.

The degree condition forces

    deg Q<=2m-2.                                    (5)

Indeed a nonzero coefficient of Q above that degree would contribute a
term of degree at least `4m-1` in `Qx^(2m)`. Neither Bx^m, of degree
at most 3m, nor A, of degree at most 2m, can cancel it for m>=2.
For m=2 there is a possible tie at degree `3m=4m-2`; it only changes
the relation between the leading Q and B coefficients and has no effect
on (5). No value of that Q coefficient is assumed below.

Put

    b0=[x^(2m)]B, r0=[x^m]R,
    D0=Q-P^2/(4N), E0=R-MP/(2N), kappa^2=N,
    C0=M/kappa.                                     (6)

The full numerator identities give

    D0=(4AQ-B^2)/(4N),
    E0=(BRx^m+2AR-2QSx^m-BS)/(2N).                  (7)

Consequently, at original x=infinity,

    D0=-b0^2 x^2/4+O(x), E0=O(x^2),
    C0=r0 x+O(1), kappa=x^(2m-1)(1+O(x^-1)),        (8)

where the chosen branch has leading sign +1. Once b0=0, the sharper
bounds are

    D0=O(1), E0=O(x), C0=r0 x+O(1).                 (9)

For D0, both numerator terms in (7) then have degree at most 4m-2.
For E0, its highest possible term BRx^m has degree at most 4m-1;
all remaining terms have degree at most 4m-2 or less. These statements
include every lower coefficient and every finite root of N.

The full numerator charts give the actual values

    H|D=b0 b+A_(2m),
    L|D=-r0 b-S_m.                                  (10)

There is no quadratic term in H|D because of (5). Thus F|D is constant
if and only if b0=r0=0: its quadratic coefficient is b0^2, and after
that vanishes its linear coefficient is -r0. This proves the stated
equivalence of boundary restrictions independently of the residue proof.

## 3. Every nonzero-index coefficient and its exact Lagrange formula

Choose the field `K=C(x)(kappa)`. If N is square this is C(x), with one
chosen polynomial square root; the reducible equation is not treated as
a connected double cover. Otherwise it is a quadratic field. In either
case the derivation extends uniquely and the formal inverse exists:

    F(x,T)=v^-4,
    T=kappa^-1 v^-1+sum_(j>=0) T_j v^j in K((v)).

Rational substitution is injective, since the highest nonzero t-degree
of a polynomial has uniquely smallest v-valuation. The chain rule gives

    partial_x G(x,T)|v=(1/4)v^5 T_v,
    ( [v^(j+4)]G(x,T) )'=(j/4)T_j.                  (11)

Thus `T_j dx` is exact in K for every j!=0. Only positive j from
(4) are needed here.

Center the quadratic by `t=y-P/(2N)`. Define

    C(z)=1+2D0 z^2+C0 z^3+(D0^2+E0)z^4.             (12)

Writing `q=1/(kappa y)` gives `v=q C(q)^(-1/4)`.
For every j>=1, formal coefficient extraction yields

    T_j= -1/(j kappa) [z^(j+1)] C(z)^(j/4).          (13)

For completeness, change variable in the formal residue for
`[v^j]q^-1`. The Jacobian factor is
`C(q)^(-1/4)(1-q C'(q)/(4C(q)))`.
The resulting coefficient is

    [z^(j+1)]C(z)^(j/4)
      -(1/4)[z^j]C'(z)C(z)^(j/4-1)
    =-(1/j)[z^(j+1)]C(z)^(j/4),

where the second equality differentiates the formal power. This proves
(13) without an analytic convergence or finite-index premise. Centering
changes only T0.

## 4. Unique top weights produce the two residues

First take `j=2m-3` and `k=m-1`. Under (8), assign bounds
`deg_infinity D0<=2`, `deg_infinity E0<=2`, and
`deg_infinity C0<=1`. In a monomial of (12), D0 uses two z-degrees
for at most two x-degrees, C0 uses three z-degrees for at most one,
and E0 uses four for at most two. The D0^2 term is still pure D0.
Hence the only contributions to x-degree `j+1=2m-2` in (13) come
from the pure-D0 part. Every other contribution has strictly smaller
x-degree.

Setting C0=E0=0 solves the centered inverse exactly as
`kappa v y=sqrt(1-D0 v^2)`. Therefore the highest numerator term of
T_j is

    binom(1/2,m-1)(-D0)^(m-1)/kappa.

By (8), its x^-1 coefficient is

    binom(1/2,m-1) (b0^2/4)^(m-1).                  (14)

This generalized binomial is nonzero for every m>=2. No higher numerator
power can mix with lower coefficients of kappa to change (14), because
`2m-2` is the maximum numerator degree. Exactness in (11) forces its
residue at infinity to vanish, hence b0=0.

Now (9) holds. Take `j=6m-7`, so `j+1=6m-6`.
A term containing q copies of C0 and ell copies of E0 has x-degree
at most q+ell and uses at least `3q+4ell` z-degrees. Any D0 adds
positive z-degree without increasing its x-degree. Attaining the
maximum x-degree `2m-2` therefore forces exactly

    q=2m-2, ell=0, and no D0 factor.

Indeed if ell>=1 then `3(q+ell)<=6m-6-ell`, so q+ell<=2m-3;
if ell=0, adding a D0 factor also prevents q=2m-2.
Thus the unique possible highest contribution in (13) is

    -binom((6m-7)/4,2m-2) C0^(2m-2)/[(6m-7)kappa].

Its x^-1 coefficient is

    -binom((6m-7)/4,2m-2) r0^(2m-2)/(6m-7).         (15)

The binomial is nonzero: its upper argument is an odd integer divided
by four, and none of its finitely many falling factors is zero.
Exactness forces r0=0. Equations (10) now prove the theorem.

The residue used in both steps is a genuine point of the field K.
For nonsquare N of even degree 4m-2, its two infinity points are
unramified, with parameter r=1/x and leading kappa signs +/-1. Choose
the + branch. For square N choose its polynomial square root with that
sign; K has just one infinity point. In both cases an x^-1 coefficient
is minus the residue. Repeated finite roots and a square N do not alter
the argument.

There is an exact stopping boundary for this particular infinity test.
After b0=r0=0, (7) gives E0=O(1), and M has degree at most 2m-1,
so C0 is bounded too. Formula (13) then shows that every nonzero-index
coefficient (T_{-1} and T_j for j>=1) is O(x^(-(2m-1))). Thus all
of their infinity residues vanish automatically. Further inverse indices
at this same place cannot strengthen the constant-D conclusion; another
place or invariant is needed. This is not an exactness or algebraization
claim for the complete differentials.

## 5. Sharp positive and negative constant-D controls for every m

Set

    y=x^(m-1)(1+x^m t), H=y^2, L=lambda y.

The proved basis shows y in L1 and H in L2. Their restrictions to D
are both zero, and `deg_t H=2` with `N=x^(4m-2)`.
For every m>=2 and every lambda,

    J(y,1/[(2m-2)x^(2m-2)])=1,
    J(y^4+lambda y,
      1/[(2m-2)x^(2m-2)(4y^3+lambda)])=1.           (16)

Thus constant D permits actual rational mates, including nonconstant L.
The mate is not claimed polynomial or globally regular.

Conversely, put

    H=y^2-t^2, L=t.

These are again global with constant restrictions to D and leading
`N=x^(4m-2)-1`. At its simple root x=1 the universal coefficient
`T2=-M/(4N)=-1/(4N)` has residue

    -1/[4(4m-2)]!=0.                                (17)

No rational mate exists. Hence constant D is not sufficient. These
controls are actual original-source families, not abstract local models.

## 6. Exact universe, replay and scope

The source derives the general numerator identities and formal coefficient
formula, checks the two top-weight integer universes for m=2..9, and
reconstructs full global boxes and boundary values for m=2..5 with all
coefficient parameters retained. It verifies (16)--(17) in original source
coordinates for the named finite range. The proof of every m is the
uniform degree and weight argument above; finite controls are independent
hostile checks, not an all-m enumeration claim.

Both normal and optimized runs pass **439 exact gates**, with
byte-identical output (**329 bytes**):

```bash
python3 04-computation/planar_jc48_sep08_all_m_constant_d.py
python3 -O 04-computation/planar_jc48_sep08_all_m_constant_d.py
```

Frozen source: **7,034 bytes**, SHA256
`7b2956263b1c03298c348e180b782cb3761f9e93e03cfcd86bbd4b818fa3b00d`.
Frozen output: SHA256
`28ddd67d3d7aac3ddb7abd508aea1be1bb857b3f7f5e0953d1d2e406570fb8f9`.
Semantic record SHA256:
`1552692247545b143e22e9391d5ad8d0ec21f67938462a72c36542bba652bd0f`.

The [full independent audit](planar_jc48_sep08_all_m_constant_d_audit.md)
accepts all439 gates,1189 alternative controls and the explicit W3
fixed-index hostile. The proof removes a mistaken
fixed-index extrapolation and replaces it by two proved moving indices;
it does not infer that their vanishing supplies a rational mate.
