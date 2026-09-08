# Independent audit: the sevenfold point at original infinity

**Status: INDEPENDENT FULL ANALYTIC / SOURCE / REPLAY AUDIT PASS.**
The [primary proof](planar_jc48_sep08_seven_one_infinity.md), standalone
producer and frozen output are accepted without correction. This sidecar
does not modify their mathematical content or promote their status.
Parent owns status promotion and checkpointing. JC(2) remains open.

## 1. Exact accepted scope and inherited mechanisms

The object is the fixed surface

    W=(P1_x x P1_z) minus {z=x^2},
    t=1/(z-x^2), x=1/r, t=-r^2-r^4 b_D,
    omega=dx wedge dt=r^2 dr wedge db_D.

For `H in L2`, `L in L1`, `deg_t H=2`, and `F=H^2+L`, suppose
the leading binary octic has a sevenfold zero at original infinity and
one simple zero at the arbitrary finite point p. After the permitted
constant leading-scalar normalization, let `u=x-p` and use the full
coefficient family in primary (1). The accepted statement is the iff

    rational Jacobian mate exists
      iff a=b=c=d=alpha=beta=gamma=delta=m0=0.

The constants e,k,ell and the original p remain arbitrary. The exact
positive family is `H=u t^2+e t+k`, `L=ell`, with rational mate
`-1/(2Ht)`. Its original polynomial Jacobian always has a factor `2H`,
so the complete class has no polynomial mate, of any degree.

This is a statement about the complete displayed fixed-DG square-prefix
class. It does not say that an arbitrary source quartic has such a
global square prefix. It does not make the rational mate globally
regular. The separate sevenfold-finite theorem is a different placement
and has its own accepted audit; this note does not replace it or assert
an all-partition compiler theorem.

I used the current proved suppliers for the
[full section filtration](planar_jc48_sep08_dg_quadratic.md),
[formal inverse chain rule](planar_jc48_sep08_leading_exactness.md),
[simple and normal-unit shared roots](planar_jc48_sep08_shared_roots.md),
and [M-unit local regimes](planar_jc48_sep08_quartic_common_root.md).
Their relevant conclusions were reread during this audit wave. The
local exhaustion added here is checked directly below. In particular,
the argument retains the original canonical factor at infinity instead
of importing an unweighted finite-point pole count.

## 2. Complete global coefficients, T2 descent, and the actual chart

For fixed leading `N=u`, the full quadratic layer has six free
coefficients. In original x coordinates, the filtration forces
`P5=P6=0`, `Q2=P4`, `Q1=P3`, and no higher Q coefficients. Re-expanding
at `x=u+p` gives exactly the primary's P,Q. The analogous linear-layer
conditions give its M,R. These identities leave six parameters for H
and six for L before the first exactness condition.

I checked the three quadratic numerator rows

    Q, P-2Qx^2, u-Px^2+Qx^4

and the two linear numerator rows directly. Their degree bounds are
four and two, respectively, for arbitrary p and all displayed
coefficients. Completeness comes from solving the whole filtration,
not merely from exhibiting a subspace. The producer independently
starts with the complete coefficient boxes, solves the unique forced
block, and recovers the same six-dimensional spaces. Its solved block
has polynomial coefficients in p and the displayed unit pivots; no
special-p rank loss is hidden in a generic matrix-rank calculation.
Writing u is an ordinary polynomial coordinate substitution. No
translation of the compactification is assumed.

The formal inverse coefficient argument is valid for rational G.
Center `H=u s^2+D` and write `F=(u s^2+D)^2+M s+E`. With
`v=F^(-1/4)` and `alpha0^2=u`, inversion begins

    s=alpha0^-1 v^-1-D/(2alpha0) v-M/(4u) v^2+O(v^3).

The original centering is independent of v and hence does not change
T2. The exact coefficientwise chain rule

    Gtilde_u=(1/4)v^5 T_v

gives `([v^6]Gtilde)'=-M/(8u)`. Formal substitution is injective on
the original rational field, including its denominators, by the proved
inverse-field argument. Trace of a primitive from the separable
quadratic coefficient field commutes with differentiation. Since M/u
already lies in `C(u)`, normalized trace makes `M du/u` rational-exact.
Its residue at zero is m0, so m0 vanishes. The equivalent local
ramified calculation multiplies the residue by the ramification index
and gives the same obstruction; it does not make a simple pole exact.

The compact infinity chart is actual, not a formal substitute. In
`q=1/z`, the boundary graph is `q=r^2`, and `s=q-r^2=1/b_D`.
Thus `(r,s)` are smooth compact coordinates at the intersection of
the deleted section and the original infinity divisor. Direct
differentiation gives

    omega=-r^2 s^-2 dr wedge ds.

I reconstructed `calN=s^2H` and `calM=sL` by literal substitution of
`u=1/r-p`, `t=-r^2-r^4/s`. The constant row of calN is precisely
`r^7(1-pr)`, its normal row is the stated B, and the complete remaining
rows C,D,E agree with the primary and source. In particular, the
successive leading coefficients of B are `-a,-b,-c,-d`, and those
of D are `-alpha,-beta,-gamma,-delta` after the preceding ones
vanish. The extra lower coefficients involving p are all retained.

On the actual generic fibre

    Phi=calN^2+s^3 calM-zeta s^4=0,

the relative form is, up to the fixed orientation sign,

    eta=r^2 s^2 dr/Phi_s.

The formulas below keep both this r squared and the order of dr on
a ramified normalization.

## 3. Compact-component consumer and local branch counts

The finite simple root has regular relative form on every normalized
generic branch, whether M is a unit or shares the point. A generic
fibre is smooth in W; excluding the finitely many critical values
makes the relative form regular there as well. A rational mate
restricts to a meromorphic function on each compact normalized generic
component, and its derivative is eta up to a fixed nonzero constant.
That equality forbids a pole wherever eta is regular.

If every boundary form is regular, the proposed primitive is
holomorphic and hence constant on each compact component. Every
generic component meets the original affine source, where the volume
is nonzero. Neither the deleted section nor the original infinity
divisor can be a component at the transcendental generic level.
Thus eta is nonzero there and exactness is impossible. This reasoning
is componentwise and does not require geometric generic integrality,
connectedness, or an irreducible polynomial presentation. A single
nonzero logarithmic residue gives an even earlier contradiction.

Let m=7, `j=ord B`, and `n=ord D`. Infinite order is allowed.
The local Weierstrass degree in s is two when j=0, three when
j>0,n=0, and four when j,n are positive. In the last case the
fourth coefficient at r=0 contains `-zeta` and is nonzero at the
generic level. These statements also cover vanishing constant C or E.

An M-unit gives the previously proved unbalanced analysis: 7 is never
3j for integral j. All unweighted forms are regular, and the r squared
factor cannot spoil that. Hence alpha vanishes. If a is nonzero,
normal-unit regularity applies to the shared point for every n,
including infinite n. It likewise gives a contradiction, so a vanishes.

For later cancelling branches, solve `calN(r,s0(r))=0` at the high
centre of order `n0=7-j`. The leading derivative B is nonzero in the
rescaled equation, so this analytic centre exists and is independent
of zeta. Put

    ell_c=ord_r(calM(r,s0)-zeta s0).

Generically `ell_c<=n0`, since the coefficient of order n0 receives
the nonzero affine term from `-zeta s0`. At most one special level
can cancel that coefficient. On the pair of determinations,

    ord_r calN=(3n0+ell_c)/2,
    eta=(unit) r^((7-3j-ell_c+4)/2) dr.

I checked the displacement from s0 has strictly higher order than
s0 in each used range. The `2 calN calN_s` term dominates Phi_s:
the remaining `3s^2(calM-zeta s)` and `s^3(calM_s-zeta)` terms
have strictly higher orders there. Perturbing s away from the centre
also cannot alter the leading order of the contact. If half-integral
powers occur, passing to the actual quadratic normalization adds the
order of dr. Nonnegative displayed exponents are therefore regular.

## 4. Independent verification of every nonfinal Newton regime

I reconstructed each leading face from the full local equation, keeping
all terms of the same weight. The following table records the simple
branches; the orders in the last column include the source r squared
and the differential of the normalized parameter.

| Regime | Branch valuation `(ord r,ord s)` | `ord Phi_s` | `ord eta` |
| --- | --- | ---: | ---: |
| j=1, low roots | (1,1) | 3 | 1 |
| n=1, low root | (1,1) | 3 | 1 |
| n=1,j=2, middle root | (1,3) | 7 | 1 |
| n=1,j>=3, high branch | (3,13) | 29 | 5 |
| j=2,n>=2, low roots | (1,2) | 6 | 0 |
| n=2,j>=3, low root | (1,2) | 6 | 0 |
| n=2, high simple roots | (1,4) | 10 | 0 |

For j=1, the two low roots arise from

    (B1+C0 Z)^2+D1 Z+(E0-zeta)Z^2.

Its constant is nonzero and its discriminant has zeta coefficient
`4B1^2`. It therefore gives two distinct nonzero roots generically,
also when D1=0. The other two determinations have centre order six
and contact between one and six. Their exponents are
`(8-ell_c)/2>=1`. This accounts for all four determinations and
forces b=0 by compact exactness.

For n=1,j>=2, the low face is

    Z^3[D1+(C0^2+E0-zeta)Z].

It gives one nonzero simple low root. For j=2, the middle face is
`Z^2(B2^2+D1 Z)`. It adds a simple order-three branch; the final
pair has centre order five and contact exactly one, giving exponent
two. For j>=3, use `r=tau^3`, `s=tau^13 Z`: the high face is
`a0^2+D1 Z^3`, with three nonzero simple Puiseux determinations.
They form one normalized branch because gcd(3,13)=1. The denominator
has order 29; the numerator has order `6+26+2=34`. These are all
remaining determinations. All forms are regular, forcing beta=0.

For j=2,n>=2, the low quadratic is the previous quadratic with
B2,D2. It has nonzero constant and discriminant slope `4B2^2`.
There are two simple order-two branches. The high centre has order
five and contact between two and five; its exponent
`(5-ell_c)/2` is nonnegative. These four determinations are regular,
forcing c=0.

For n=2,j>=3, the low face is

    Z^3[D2+(C0^2+E0-zeta)Z].

It gives one simple order-two branch. At the high scale `s=r^4 Z`,
j>=4 gives the simple cubic `a0^2+D2 Z^3`, accounting for the
other three determinations. For j=3 the whole high cubic is instead

    T(Z)=(a0+B3 Z)^2+D2 Z^3.

It has no zero root. A repeated root must be

    Z0=-3a0/B3, D2=4B3^3/(27a0).

Direct evaluation gives `T(Z0)=T'(Z0)=0` and
`T''(Z0)=-2B3^2/3`, a nonzero value. Thus it cannot be triple;
the other cubic root is simple. No root multiplicity is inferred
solely from a finite discriminant sample.

Let `Q(r,Z,zeta)=Phi(r,r^4 Z)/r^14`. It is analytic, including
every higher-order local coefficient, and exactly

    partial_zeta Q=-r^2 Z^4.

At the nondegenerate critical centre `Q_Z=0`, implicit differentiation
therefore gives `partial_zeta Qcritical=-r^2 Zcentre^4`. Its order-two
coefficient has nonzero slope `-Z0^4`. The generic critical contact
is one or two, never larger. This is the needed all-coefficient
contact bound; it is not a guess from the cubic face alone.

Analytic Morse coordinates give, up to a unit and an orientation sign,
`v_M^2=R(r,zeta)`. Here `Phi_s=r^10 Q_Z`, whereas the numerator
is `r^10 Z^2 dr`. Hence eta is a holomorphic unit times `dr/v_M`.
At contact one the normalization has `r=tau^2`, `v_M` of order
one, and the form is regular of order zero. At contact two there
are two smooth branches and a nonzero logarithmic coefficient on
each. The unit cannot remove their residues. The simple high root
and the low root have regular forms. Thus either all forms are
regular, or a nonzero residue already excludes exactness. This
forces gamma=0.

All higher terms omitted in the face calculations have strictly higher
weight. The analytic high centres and the exact critical-value
derivative handle the only possible further contacts. The source's
representative integer normal orders do not truncate the proof:
higher or infinite orders only increase the omitted weights.

## 5. Final tangent logs, actual rational field, and sharpness

After the preceding forced vanishings and m0=0, the exact original
family is

    H=u t^2+(du+e)t+k, L=delta u t+ell.

Substituting the full p-dependent compact coefficients and setting
`s=r^3 Z` gives the leading polynomial

    P_zeta(Z)=Z^2[(-d+kZ)^2-delta Z+(ell-zeta)Z^2].

For d nonzero the bracketed quadratic has nonzero constant and
discriminant slope `4d^2`; generically it has two simple nonzero
roots. For d=0,delta nonzero, its nonzero root is
`Z0=delta/(k^2+ell-zeta)`, again simple. In either case the
implicit-function theorem produces an actual smooth branch with
parameter r. The numerator order is eight and the denominator
order is nine. Its nonzero logarithmic coefficient is, up to the
fixed orientation sign,

    Z0^2/P_zeta'(Z0).

This gives a genuine residue obstruction. It is legitimate to leave
the remaining zero-root determinations unclassified once one such
branch has been established. Thus d=delta=0 is necessary. This
includes the constant-L but d-nonzero regime, which a nonconstant-L
field argument alone would miss.

I also independently checked the supplemental exact field map for
delta nonzero. With `v=u t`, `h=H`, its rational inverse is

    t=(h-dv-k)/(v+e), u=v/t,
    omega=dv wedge dh/(h-dv-k).

Neither denominator is identically zero for any coefficient choice.
The inverse is literal, so the field is actually `C(v,h)`. Setting
`f=h^2+delta v+ell` then gives the actual field `C(f,h)` and
relative differential

    dh/[d h^2+delta h-d(f-ell)-delta k].

For d=0 it has residue `1/delta` at h=k. For d nonzero its
discriminant is linear in the transcendental f with nonzero
coefficient `4d^2`; the two poles are simple generically and their
residues are nonzero. There is no geometric-integrality or selected
branch assumption. This separately confirms the delta-nonzero
exclusion while retaining the actual source volume.

When all the claimed necessary coefficients vanish,
`H=u t^2+e t+k`, and direct original differentiation gives

    J(H,-1/t)=1, J(H^2+ell,-1/(2Ht))=1.

The full global coefficient rules already show this H belongs to L2;
direct second-chart substitution is polynomial in r,b_D for every p.
H is nonconstant, so no original polynomial mate survives the factor
`2H`. The rational positive family proves the sharp iff rather than
merely marking a residual stratum as undecided.

## 6. Source universe, replays and independent controls

I read the complete
[standalone producer](../../04-computation/planar_jc48_sep08_seven_one_infinity.py).
It imports no inherited mathematical implementation. All gates use
always-active exceptions. Its finite symbolic universe comprises:

* the complete H,L coefficient boxes and forced global sections;
* the short unsaturated T2 inverse calculation and trace-residue factor;
* full arbitrary-p infinity numerators and source volume;
* all leading faces that change with j,n, including the actual (3,13)
  normalization, with higher orders handled analytically;
* contact ranges `j=1:1..6`, `j=2:2..5`, the middle contact one,
  and both possible Morse contacts;
* the nontriple cubic, exact generic-value derivative, final tangent
  logarithms, original-source rational field derivation and sharp iff.

The finite contact checks support, but do not replace, the analytic
centre and generic-value proofs. In particular, a one-line arithmetic
valuation gate is not being treated as an independent branch-existence
theorem. Likewise the displayed full matrix rank supplements the
explicit all-p forced coefficient solution.

Literal hostiles include `H=u t^2+u t`, `L=0`, for which the
generic-level-one tangent polynomial is `Z^2(1-Z^2)` and the
two stated residue coefficients are `-1/2,1/2`; and the
nonconstant-L control `H=u t^2`, `L=u t`, with a nonzero tangent
residue. The positive family retains arbitrary e,k,ell and p.
No finite bound on a mate's degree or pole divisor is assumed.

Fresh independent normal and optimized runs both exited zero and
passed **113 gates**. Both are byte-identical to the **317-byte**
frozen output. Verified pins are:

| Artifact | Bytes | SHA256 |
| --- | ---: | --- |
| Primary before status promotion | 15487 | `766fe916be58eeffb371217ba4f2899a96df1976c24a2c11f0366f6663aedfeb` |
| Source | 10179 | `f1db3c2338a9f49703ec250bd8ee52d9de432cde46a149446ddcb0d56c1593cc` |
| Frozen, independent normal, independent optimized outputs | 317 each | `faae700d8a80809f2b6155285b15b17a01fda57a580d814f6db29f1c082e6043` |

Semantic record SHA256:
`730eb5fee6ebc4b630fa9546698e0229566d2f6b9c42ae5844e905e1b930b5a6`.

Reproduce from the worktree root:

```sh
python3 -B 04-computation/planar_jc48_sep08_seven_one_infinity.py
python3 -B -O 04-computation/planar_jc48_sep08_seven_one_infinity.py
```

An additional independent temporary reconstruction, written before
reading the frozen producer, passed **26 controls**. It starts from
the original H,L, differentiates the original chart volume, checks
the five global numerator rows, reconstructs the Newton faces, checks
the cubic double and second derivative, differentiates the scaled
equation in zeta, reconstructs the exact final tangent polynomial,
and verifies the actual rational inverse, its volume and both sharp
Jacobian identities. The calculations and their mathematical formulas
are described above; the temporary script is not a new repository
dependency. No producer source, output or primary was changed.

**Final acceptance:** complete fixed-location coefficient space,
primitive-field descent, all normalized local determinations,
generic simple/double faces, contact bounds, retained source volume,
compact-component exactness, rational iff, polynomial consequence,
standalone source, and both replay modes PASS. No correction remains.
