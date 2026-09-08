# The sevenfold infinity point: an exact rational-mate classification

**Status: PROVED + FINITE-EXACT + INDEPENDENTLY AUDITED.**
This is a complete fixed-DG square-prefix coefficient theorem. It is not
an exclusion for arbitrary quartic source polynomials or JC(2).

## 1. Statement, complete coefficients and connection contract

Use the fixed surface and the original source volume

    W=(P1_x x P1_z)\{z=x^2},  t=1/(z-x^2),
    x=1/r,  t=-r^2-r^4b_D,  omega=dx wedge dt=r^2 dr wedge db_D.

Let `H in L2`, `L in L1`, `deg_t H=2`, and `F=H^2+L`. Suppose
its leading binary octic has a sevenfold zero at original infinity and
one simple finite zero `x=p`. Scale the nonzero leading scalar and
write `u=x-p`; this is a polynomial source coordinate, not an asserted
translation automorphism of W. The full original global coefficient
rules from [the proved filtration](planar_jc48_sep08_dg_quadratic.md) give

    H=u t^2+P t+Q,
    P=a u^4+b u^3+c u^2+d u+e,
    Q=a u^2+(b-2pa)u+k,
    L=M t+R,
    M=alpha u^4+beta u^3+gamma u^2+delta u+m0,
    R=alpha u^2+(beta-2p alpha)u+ell.                 (1)

Every parameter in (1) is arbitrary over C. In particular, all induced
constant and linear rows are retained. Completeness is also visible in
the original x equations `P5=P6=0`, `Q2=P4`, `Q1=P3`,
`Q_j=0` for `j>2`, and `R2=M4`, `R1=M3`, `R_j=0` for `j>2`;
expanding `x=u+p` gives exactly (1).

**Theorem.** A rational `G in C(x,t)` with `J(F,G)=1`
exists if and only if

    a=b=c=d=alpha=beta=gamma=delta=m0=0.             (2)

In that case

    H=u t^2+e t+k,  L=ell,
    G_H=-1/t,  J(H,G_H)=1,
    G=-1/(2Ht),  J(H^2+ell,G)=1.                    (3)

Consequently no member of the complete class (1) has a polynomial
mate of any degree: under (2), its polynomial Jacobian has the
nonconstant factor `2H`. The statement applies equally to any nonzero
constant Jacobian after scaling the proposed mate.

The inherited mechanisms are the
[formal leading-field chain rule](planar_jc48_sep08_leading_exactness.md),
[finite simple-root and normal-unit regularity](planar_jc48_sep08_shared_roots.md),
and [the complete M-unit Newton classification](planar_jc48_sep08_quartic_common_root.md).
The closest weighted-infinity predecessor is
[the two (5,3) infinity placements](planar_jc48_sep08_five_three_infinity.md).
The corrected near miss is to transport unweighted finite valuations
through inversion, or to assume a rational mate remains polynomial.
The sharp hostile (3) retains both the full global carrier and every p.

The live concepts are the full section space, the actual canonical
factor `r^2`, branch exhaustion, logarithmic residues, and the
rational/polynomial boundary. The map is restriction of the original
relative differential to the normalization of each compact generic
fibre component. It preserves exactness and labelled local orders;
forgetting the volume or an entire component destroys the needed
predicate. Section 5 supplies a separate birational-field check of the
last nonconstant-L residual, without replacing the preceding local
exhaustion. Only the needed T2 inverse coefficient is used; no claim
that more formal rows suffice for algebraization is made.

## 2. Generic-fibre consumer and the exact infinity chart

Under a rational mate, `M du/u` is rational-exact. To recall the exact
coefficient argument, solve `F(u,T)=v^-4` in `C(u)(sqrt(u))((v))`.
Completing the quadratic square, or substituting the first four inverse
terms, gives

    T2=-M/(4u),
    ([v^6]G(u,T))'=-M/(8u).                         (4)

This is the same formal chain rule as in the proved leading-field note,
not an assumption about convergence. Normalized field trace to C(u)
commutes with differentiation, so `M du/u` has zero residue at u=0.
Thus `m0=0`. This also follows by its doubled residue on the ramified
quadratic cover.

The finite simple boundary point has regular relative differential on
every normalized generic branch, including shared-root degenerations.
The surface W is smooth; its generic fibre is smooth in W, and omega is
regular there. Hence the relative form is regular everywhere in W,
including its possible zeros on D. An ambient rational mate restricts
to a meromorphic function on each compact normalized generic component,
whose differential is the relative form up to one fixed sign. A pole
of that function creates a pole of its differential. If all boundary
forms are regular, the function is constant on each compact component,
which is impossible: no generic component is D or the deleted section,
so every component meets the original affine source where omega is
nonzero. No connectedness or geometric irreducibility assumption is
used in this consumer. A single nonzero logarithmic residue is already
a contradiction before applying it.

For infinity use the actual chart `b_D=1/s`. The local equation and
relative form are exactly, up to the harmless convention sign,

    calN=s^2 H(1/r-p,-r^2-r^4/s),
    calM=s L(1/r-p,-r^2-r^4/s),
    Phi=calN^2+s^3 calM-zeta s^4,
    eta=(unit) r^2 s^2 dr/Phi_s.                     (5)

Here zeta is the original transcendental fibre value. Write

    calN=r^7(1-pr)+s B(r)+s^2 C(r),
    calM=D(r)+s E(r).                                (6)

The exact coefficient identities, not just their orders, are

    B=-r^4 P(1/r-p)+2r^5-2pr^6,
    C=-a p^4r^2+4a p^3r-3a p^2
      +b p^3r^2-3b p^2r+2bp
      -c p^2r^2+2cp r-c+dpr^2-dr-er^2+k+r^3-pr^4,
    D=-r^4 M(1/r-p),
    E=-alpha p^4r^2+4alpha p^3r-3alpha p^2
      +beta p^3r^2-3beta p^2r+2beta p
      -gamma p^2r^2+2gamma pr-gamma+delta pr^2-delta r+ell.
                                                               (7)

All are polynomials. In particular, successive vanishing of a,b,c
makes the normal order `j=ord B` successively 0,1,2 and then at least 3;
successive vanishing of alpha,beta,gamma does the same for
`n=ord D`. Infinite order is allowed. The factor `r^2` in (5) is
essential. It comes directly from `omega=-r^2 s^-2 dr wedge ds`,
not from a guessed normalization at a finite boundary point.

## 3. Exhausting the first normal and lower-section orders

We prove that a rational mate forces

    a=b=c=alpha=beta=gamma=0.                        (8)

In the local arguments set `A=r^7 a0(r)`, `a0(0)!=0`, and allow all
higher coefficients in B,C,D,E from (6), or arbitrary analytic ones.
All orders below include the actual r^2 factor. If `j=0`, the generic
local Weierstrass degree in s is two. If `j>0,n=0` it is three. If
`j,n>=1`, it is four because its fourth coefficient at r=0 contains
`-zeta`. Thus the stated numbers of determinations pay every branch.

First, `alpha!=0` gives an M-unit. The unweighted classification is
regular in all unbalanced cases, and `7!=3j` for every integral j.
This includes j=0 and infinite j; adding r^2 only improves the orders.
All compact forms would therefore be regular, a contradiction.
Hence alpha=0. If `a!=0`, normal-unit regularity applies for every n,
so the same argument gives a=0.

For reference, when `m=7`, a cancelling pair centred at
`s0(r)~r^(m-j)` has generic contact

    ell_c=ord_r(calM(r,s0)-zeta s0)<=m-j.             (9)

The centre is the analytic solution of calN=0 and is independent of
zeta. Its nonzero leading coefficient pays this upper bound from the
term `-zeta s0`. In each application below the perturbation has higher
order than s0. Since `N_s~r^j`, the two determinations have

    eta=(unit) r^((m-3j-ell_c+4)/2) dr.              (10)

Indeed `N^2=-s^3(calM-zeta s)` gives
`ord N=(3(m-j)+ell_c)/2`, and `2N N_s` strictly dominates
Phi_s. A nonnegative exponent in (10) remains regular after quadratic
normalization, including the extra order of dr.

### 3.1 Normal order j=1

Suppose b!=0. There are two low smooth branches `s~r`, each of
weighted order one. Their nonzero roots are those of

    (B1+C0 Z)^2+D1 Z+(E0-zeta)Z^2,                 (11)

with B1!=0. Its constant coefficient is nonzero and its discriminant
has a nonzero zeta coefficient; both roots are nonzero and simple
for the generic fibre. The remaining two determinations cancel calN
at a centre of order six. They have `1<=ell_c<=6`; (10) gives
exponent `(8-ell_c)/2>=1`. Thus all four determinations are regular,
and the compact consumer forces b=0.

### 3.2 Lower-section order n=1

Now j>=2. Suppose beta!=0. There is one low smooth branch `s~r`,
weighted order one; its nonzero face is
`D1+(C0^2+E0-zeta)Z`, with D1!=0.

If j=2, a second smooth branch has `s~r^3`, with nonzero root of
`B2^2+D1 Z`, and weighted order one. The last two determinations
cancel calN at order five. Their contact is exactly one, so (10)
gives exponent two. This counts all four determinations.

If j>=3, the remaining three determinations have balance
`s~r^(13/3)` with simple face `a0^2+D1 Z^3`. They form one
normalized branch with `(ord r,ord s)=(3,13)`. Its denominator
Phi_s has order 29 and its numerator, including dr, has order
`6+26+2=34`; hence eta has order five. All forms are again regular,
forcing beta=0.

### 3.3 Normal order j=2 with n>=2

Suppose c!=0. Two smooth low branches have `s~r^2`; the same
quadratic (11), now with B2,D2, gives simple nonzero roots.
Their weighted order is zero. The remaining two determinations cancel
calN at order five and have `2<=ell_c<=5`. Formula (10) gives
exponent `(5-ell_c)/2>=0`. These four determinations are regular,
forcing c=0.

### 3.4 Lower-section order n=2 with j>=3

Suppose gamma!=0. One low smooth branch has `s~r^2`, weighted
order zero. The remaining three determinations have `s~r^4`.
For j>=4 their cubic face is `a0^2+D2 Z^3`, simple, and eta
has order zero on every branch. For j=3 the full face is

    C_*(Z)=(a0+B3 Z)^2+D2 Z^3.                       (12)

It has nonzero roots and no triple root. A double root, if present,
has `Z0=-3a0/B3`, `D2=4B3^3/(27a0)` and nonzero second
derivative. The remaining root is simple. After `s=r^4 Z`, divide
Phi by r^14. Near the double root, parameterized analytic Morse
reduction has equation `v_M^2=R(r,zeta)`. Differentiating its
critical value in zeta gives a unit times r^2, because the term
`-zeta s^4` has order sixteen and `Z0!=0`. Therefore the generic
contact is at most two; higher analytic terms cannot erase this
coefficient identically in zeta.

The weighted relative form is a unit times `dr/v_M`. At contact one,
`r=tau^2`, `v_M~tau`, and its normalized order is zero. At contact
two there are two smooth branches with a nonzero `dr/r` coefficient.
Those nonzero logarithmic residues forbid any rational primitive.
Consequently, under a hypothetical mate, the only allowed case would
have all forms regular, which the compact consumer also forbids.
This forces gamma=0 and proves (8).

Each rescaling used a nonzero simple face, the paid analytic centre,
or the explicit nontriple cubic (12). All omitted terms have higher
weight in the corresponding chart; (9) and the critical-value derivative
pay the only possible further contacts. Thus this is an all-coefficient
classification of the specified regimes, rather than a list of selected
Newton monomials.

## 4. The remaining tangent logs and sharp sufficiency

After (8), the actual original polynomials are

    H=u t^2+(du+e)t+k,  L=delta u t+ell.              (13)

Their exact infinity coefficients are

    B=-d r^3+(dp-e)r^4+2r^5-2pr^6,
    C=k-dr+(dp-e)r^2+r^3-pr^4,
    D=-delta r^3+p delta r^4,
    E=ell-delta r+p delta r^2.

Set `s=r^3 Z`. The leading face of `Phi/r^12` is

    P_zeta(Z)=Z^2[(-d+kZ)^2-delta Z+(ell-zeta)Z^2].   (14)

If d!=0, the bracketed quadratic has two simple nonzero roots for
generic zeta: its constant is d^2 and its discriminant has nonzero
zeta coefficient `4d^2`. If d=0 but delta!=0, its nonzero part
has the simple root `Z=delta/(k^2+ell-zeta)`.
For each such root, the analytic implicit-function theorem gives an
actual smooth normalized branch parametrized by r. Formula (5) gives

    eta=(nonzero constant) dr/r+regular terms,
    coefficient = (unit at 0) Z0^2/P_zeta'(Z0).      (15)

This residue cannot vanish. We do not need to classify the other
roots Z=0 once one actual logarithmic branch has been exhibited.
Thus every rational mate forces d=delta=0.

Conversely, the resulting family is exactly (3), whose displayed
Jacobian is the rational identity 1 for every p,e,k,ell. Both H and L
are globally regular in W: the complete source coefficient equations
already prove this, and direct substitution gives a polynomial in
r,b_D. The rational function (3) need not be globally regular, and
cannot be polynomial because `J(H^2+ell,Gpoly)=2H J(H,Gpoly)`.
This proves necessity, sufficiency, and the unrestricted polynomial
exclusion without confusing their scopes.

## 5. Independent field consumer for delta!=0

The last nonconstant-L residual has a second exact obstruction that
retains the original volume. Put

    v=u t,  h=H=t(v+e)+dv+k.

Then

    t=(h-dv-k)/(v+e),  u=v/t,
    C(u,t)=C(v,h),
    omega=dv wedge dh/(h-dv-k).

These are identities of function fields; neither denominator is
identically zero, for any coefficient choice. If delta!=0, putting
`f=h^2+delta v+ell` gives `C(u,t)=C(f,h)`. Under our Jacobian
convention the primitive on the generic fibre would have differential

    dh/[d h^2+delta h-d(f-ell)-delta k].              (16)

For d=0 this has a nonzero simple residue. For d!=0 the generic
quadratic denominator has two distinct simple roots and its numerator
is nonzero; again both residues are nonzero. This separately excludes
all delta!=0 cases. It does not replace the weighted-log proof for
constant L with d!=0, and it is not a surface automorphism claim.

## 6. Exact controls, explicit universe, and remaining scope

The companion source verifies the full original coefficient identities
and the complete infinity numerator; the exact faces for every regime;
all cancelling-contact bounds in their named finite integer ranges;
the balanced cubic's double/no-triple identity and generic-value order;
the remaining logarithmic residues; and independent original-source
Jacobian controls for both (16) and the sharp rational family (3).
It retains nonzero p and every lower coefficient symbolically.

The finite universe is the displayed face regimes, contact ranges
`j=1:1..6`, `j=2:2..5`, the two Morse contacts, and the stated literal
positive and hostile coefficient controls. These computations support
the analytic all-parameter proof; they do not replace it with a finite
coefficient census or a bound on the mate degree. The sharp positive
family proves that replacing the polynomial conclusion by unconditional
rational exclusion would be false. A d!=0 constant-L control supplies
a literal logarithmic failure just outside that family.

Both normal and optimized replays pass **113 exact gates** and produce
byte-identical output (**317 bytes**):

```bash
python3 04-computation/planar_jc48_sep08_seven_one_infinity.py
python3 -O 04-computation/planar_jc48_sep08_seven_one_infinity.py
```

Frozen source: **10,179 bytes**, SHA256
`f1db3c2338a9f49703ec250bd8ee52d9de432cde46a149446ddcb0d56c1593cc`.
Frozen output: SHA256
`faae700d8a80809f2b6155285b15b17a01fda57a580d814f6db29f1c082e6043`.
Semantic record SHA256:
`730eb5fee6ebc4b630fa9546698e0229566d2f6b9c42ae5844e905e1b930b5a6`.

The [full independent audit](planar_jc48_sep08_seven_one_infinity_audit.md)
accepts all113 gates,26 alternative original-coordinate controls,
complete normalized branch exhaustion and the exact rational iff.
Both this theorem and the separate finite-sevenfold theorem are now
audited, completing every (7,1) placement for polynomial exclusion.
