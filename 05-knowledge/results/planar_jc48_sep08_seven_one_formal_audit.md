# Independent audit of the entire (7,1) inverse hierarchy

**Status: PROVED / INDEPENDENT AUDIT PASS.**
The complete analytic proof, frozen source, and independent normal and
optimized replays of [the primary](planar_jc48_sep08_seven_one_formal.md)
are accepted. I independently reconstructed the universal opposite-pair
resolvent identity and a rational binomial truncation of the actual
hostile's formal mate. No mathematical or source correction was needed.
Producer artifacts were not edited; root owns promotion and integration.

## 1. Accepted equivalence and its exact scope

The object is the complete actual global square-prefix family
`F=H^2+L`, `H in L2`, `L in L1`, with two distinct finite
leading binary roots of multiplicities seven and one. In the
normalized polynomial coordinate `u=x-p`, the primary retains every
coefficient of its full displayed rows, including arbitrary `P0,P1`.
It imposes no active-point, submersion, mate, or finite-jet hypothesis.

Write `N=u^7(u-1)`, `kappa=sqrt(N)`, and
`K=C(u)(kappa)`. For the unique inverse with

    F(u,T(u,v))=v^-4,
    T=kappa^-1 v^-1+T0+T1v+T2v^2+...,

the following are equivalent:

1. `M(1)=0`.
2. `T_j du` is exact in `K` for `j=-1` and every integer `j>=1`.
3. There is a formal Laurent mate in `C(u)((1/t))`.

A formal mate can then be chosen in `t^-3 C(u)[[t^-1]]`.
There is no requirement on exactness of `T0`; its multiplier in
the formal differential equation is zero. The explicit global
hostile satisfies the entire hierarchy and admits such formal
mates but no rational mate, and consequently no algebraic formal
mate either. None of these statements is promoted to a rational-
mate criterion for the complete stratum or to a Jacobian-conjecture
theorem.

## 2. Complete global rows and both infinity points

I checked the complete original-coordinate coefficient equations
against the displayed shifted rows. In particular,

    P=2u^6-(4p+2)u^5+a u^4+b u^3+c u^2+d u+e,
    Q=u^4-(4p+1)u^3+(a+4p^2)u^2
                    +(b-2ap+4p^2)u+f,

and the complete linear rows for `M,R` retain arbitrary `p` and
all lower parameters. This is an intersection with the actual
global section space in original `x`, not a claim that the shift
acts on the compactification. No hypothesis sets `e` or `d` to
zero.

Center by `t=y-P/(2N)` and define

    D0=Q-P^2/(4N), E0=R-MP/(2N),
    C0=M/kappa, B0=M^2/N.

The full global equations give exactly

    deg(P^2-4NQ)<=8,
    deg(2NR-MP)<=8,
    deg(M^2)<=8.

Thus `D0,E0,C0,B0` are bounded at both points above infinity;
`kappa~+u^4` and `kappa~-u^4` are separate unramified branches.
The first cancellation removes degrees twelve through nine and
the second removes degrees ten and nine. These bounds use the
induced `Q,R` rows and would not hold for a free quartic prefix.

The exact centered equation for `W=kappa v y` is

    W^4+2D0 v^2 W^2+C0 v^3 W+(D0^2+E0)v^4=1,
    W(0)=1.

At every order the new coefficient occurs as `4w_n` plus a
polynomial in previously determined coefficients and `D0,E0,C0`.
Its invertible constant coefficient proves by induction, without
an index cutoff, that every `w_n` lies in `C[D0,E0,C0]`.
For `j=-1` and all `j>=1`, `T_j=w_(j+1)/kappa`, so

    T_j=O(u^-4)

at both infinities. Centering changes only `T0`. All other
possible poles of these coefficients are above zero or one,
since the recurrence introduces no denominator away from those
of `N,kappa`.

## 3. Ramified parity and all odd indices

The actual field is rational and connected, with

    z^2=(u-1)/u, u=1/(1-z^2), kappa=u^4 z.

The two finite roots are ramified: `u=1` is `z=0` and `u=0`
is `z=infinity`. The two points over `u=infinity` are `z=+1,-1`.
The deck involution is `z->-z`; it fixes `u` and changes the
sign of `kappa`.

Uniqueness of the formal inverse gives `tau(T(v))=T(-v)`.
The leading terms agree after this substitution, and both solve
the same equation, which pays the use of uniqueness. Hence
`tau(T_j)=(-1)^j T_j`. For every odd index, including `-1`,
the differential `T_j du` is anti-invariant. At a ramified point
choose a parameter `xi` with `tau(xi)=-xi`. If
`alpha=a(xi)d xi` is anti-invariant, then `a(-xi)=a(xi)`.
Its Laurent coefficient has only even exponents and cannot
contain exponent minus one. The residue at either ramified
point is therefore zero.

At both unramified infinities the differential is regular by
the preceding all-index bound, and there are no other poles.
Ordinary partial fractions on the rational `z` line now give a
primitive in `K`. This proves every odd-index exactness claim
independently of `M(1)`.

I independently checked the literal leading primitive

    d[2z-(4/3)z^3+(2/5)z^5]=du/kappa.

The two infinity points are neither merged nor discarded, and
anti-invariance is used for differentials, including the sign
of `d xi`, rather than only for their scalar coefficients.

## 4. The opposite-pair recurrence and all even indices

The four centered branches are `y(v),y(-v),y(iv),y(-iv)`.
Their leading terms are distinct and enumerate the four roots
of the centered quartic. Their sum is zero, so the centered
coefficient of every index divisible by four vanishes. For
positive such indices this is `T_(4ell)=0`.

I independently derived the full opposite-pair identity, rather
than inferring it from the first few coefficients. The monic
centered quartic is

    y^4+B y^2+C y+D=0,
    B=2D0/N, C=M/N^2,
    D=(D0^2+E0-v^-4)/N^2.

Factor it into opposite-root quadratics. If `S` is one pair sum,
their constant terms satisfy

    U+V-S^2=B, S(U-V)=C, UV=D.

Thus `X=S^2` satisfies the exact resolvent

    X(B+X)^2-C^2-4DX=0.

Set `s=v^4` and `S=-M v^2 V(s)/(2N)`. After substitution,
the resolvent becomes exactly

    (1-E0s)V^2+(D0 B0/4)s^2V^4
                   +(B0^2/64)s^3V^6=1,
    V(0)=1.

This yields

    T_(4ell+2)=-M/(4N) [s^ell]V(s).

The algebraic reduction can first be performed with `M!=0`;
the resulting universal polynomial coefficient identities extend
to `M=0`. In that case the opposite pair sum is zero directly,
and every displayed even coefficient also has the zero prefactor.
No allowed zero lower row is lost through division in an
intermediate derivation.

The recurrence coefficient of each new `V` coefficient is two.
Induction at every order therefore puts it in

    C[E0,D0 B0,B0^2].

This particular ring is essential: a naked power of the possibly
singular `D0` cannot enter. When `M(1)=0`, the simple-root
orders are

    ord_1 D0>=-1, ord_1 E0>=0, ord_1 B0>=1.

All three generators of that ring are regular at one, as is
the prefactor `M/N`. Thus every even coefficient is a rational
function of `u` regular at one. Its differential is regular
at infinity by Section 2 and can have a pole only at zero on
`P1_u`. The residue theorem forces its remaining residue to
vanish, so rational partial fractions give a primitive already
in `C(u)`. This proves sufficiency for every even index.

Conversely the literal coefficient `T2=-M/(4N)` has rational
residue `-M(1)/4` at one. Pullback to the ramified normalized
point doubles it to `-M(1)/2`. Exactness in `K` forces this
to vanish. I independently recomputed this pullback for the
minimal failed-condition control `M=1`. Consequently no weaker
condition replaces `M(1)=0` in the stated equivalence.

## 5. Formal composition, exactness, and descent

For any mate written in inverse coordinates, differentiating
`F(u,T)=v^-4` at fixed `v` and using the Jacobian convention
gives exactly

    partial_u Gtilde |_v=(1/4)v^5 partial_v T.

The sign follows from `F_t T_v=-4v^-5` and
`Gtilde_u=-1/F_t`. Its coefficient at `v^(j+4)` is
`(j/4)T_j`. The term with `j=0` vanishes, so there is no
hidden exactness condition on `T0`.

When all required differentials have primitives in `K`, choose
those primitives coefficientwise and form

    Gtilde in v^3 K[[v]],
    g_(j+4)'=(j/4)T_j, g4=0.

There is no infinite compatibility constraint between independent
coefficients beyond these equations. The formal inverse

    v=F(u,1/q)^(-1/4) in q K[[q]],
    v=kappa^-1 q+O(q^2), q=1/t,

has nonzero linear term. Composition is well-defined: only
finitely many terms contribute to each coefficient. The formal
chain rule gives a mate `Ghat` in `q^3 K[[q]]`.

The two derivations are `partial_u` at fixed `q` and
`partial_t=-q^2 partial_q`. The extension of the first to `K`
is unique in characteristic zero and commutes with its deck
involution; the second commutes with it coefficientwise. Hence
the normalized coefficientwise trace

    (Ghat+tau(Ghat))/2 in q^3 C(u)[[q]]

still has Jacobian one. This proves descent to the original
coefficient field. It uses neither analytic convergence nor
an unjustified interchange of infinite coefficient sums.

Conversely an arbitrary formal Laurent mate has a finite lower
bound in `q`. Substitution of `q=1/T(u,v)` is therefore valid
in `K((v))`, including negative exponents. The same chain
identity and coefficient extraction supply primitives of every
required `T_j`. This proves the full equivalence, not only a
necessary direction.

If such a formal mate were algebraic over `C(u,t)`, it would
generate a finite separable field extension. Both derivations
extend uniquely and agree with the formal derivations on that
element. Normalized field trace commutes with them and would
produce a rational mate. Thus the later rational obstruction
also excludes algebraic formal mates; transcendence is a
proved consequence, not a statement about the appearance of
a computed series.

## 6. The actual global hostile

For

    H=u^3(u-1)(1+u^2t)^2, L=(u-1)t,

the displayed original global rows hold exactly at `p=0`.
The normal coefficient at zero has order five, not an
assumed active order, and `M(0)=-1` is a unit. The condition
`M(1)=0` holds. The exact second-chart formulas are

    H=(1-r)b_D^2,
    L=-r(1-r)(1+r^2 b_D), F|D=b_D^4.

The finite sevenfold point is the unbalanced M-unit case. Its
three small determinations `s~u^(14/3)` form one normalized
branch with `(ord u,ord s)=(3,14)`. Since `Phi_s` has order
28, `s^2 du/Phi_s` has order two. The simple shared root is
regular by the proved local lemma. There is no root of the
leading octic at infinity, while all generic intersections
with retained `D` have zeros from the actual canonical factor
`r^2`. Ordinary points on the smooth generic fibre have unit
relative form.

The complete relative form is therefore holomorphic and
nonzero on every compact normalized generic component. A
rational primitive could have no pole, since differentiating
a pole creates a pole. It would be constant on that compact
component, a contradiction. This argument is componentwise;
it imposes no unproved generic-connectedness hypothesis.
The example consequently has all the formal properties above
but no rational or algebraic formal mate.

## 7. Independent finite controls and frozen replays

The full producer source was read. Its explicit universe is
symbolic: the complete shifted global rows, all cancelled top
coefficients, the entire `M(1)=0` hyperplane, rational
parametrization, universal recursion controls, exact chart
identities for the hostile, and a descended truncation through
`q^7`. Every gate raises an exception and remains active under
optimization. The finite recursions are controls; the
all-index conclusions are proved by the unrestricted ring
inductions and parity/residue arguments above.

Separate referee calculations, without importing the producer,
reconstructed the complete resolvent substitution, the leading
primitive, and the ramified simple-root residue. I also
constructed the hostile's rational truncation directly by
binomial expansion. Put

    rho=(u-1)/u, s=u^-2, ell=(u-1)/N^2,
    R3=-(8u^2+4u+3)/[30u^13(u-1)],
    R6=1/(48u^6N^3),
    R7=u^-28[1/(8rho^2)-9/(40rho)+9/56-rho/24].

The coefficients of `q^3,...,q^7` are respectively

    R3,
    -3sR3,
    6s^2R3,
    (-10s^3-3ell/4)R3+R6,
    (15s^4+21s ell/4)R3-6sR6+R7.

Direct original-coordinate differentiation of this rational
truncation gives constant coefficient one and coefficients
zero at `q^1,q^2,q^3,q^4`, independently verifying
`J(F,Gtrunc)=1+O(q^5)`. This calculation is a finite control
of the declared descent, not the proof of its all-order
existence.

Independent full replay commands were

```sh
python3 -B 04-computation/planar_jc48_sep08_seven_one_formal.py > /tmp/seven-one-formal-audit-normal.out
python3 -B -O 04-computation/planar_jc48_sep08_seven_one_formal.py > /tmp/seven-one-formal-audit-optimized.out
```

Both completed successfully with **93 always-active gates** and
are byte-for-byte identical to the frozen output. The full
primary, source, output, and final footer pins were reread
and checked before acceptance.

| Artifact | Bytes | SHA256 |
| --- | ---: | --- |
| Primary before promotion | 14,236 | `af64f60d83fa910e6a346ca206faf13b630bc4e6b52d9e44625d459fef100da4` |
| Frozen source | 7,342 | `723dcb0a947295ceec0483755fcc8487476139ac69513063dd81319700fa8453` |
| Frozen output and each independent replay | 407 | `dec45c98a41a9b4ab1443a16137ec231a9f44722149517e2c22fb76e98d6e522` |

The semantic record digest is
`af2ddce9f9de8d962728c9bdeb5ab5d4a0a4cad6e5ba21c06f6ea5c23f92ff2e`.
No correction remains. This audit accepts the all-index formal
equivalence, the coefficient-field descent, and the actual
global failure of rational or algebraic algebraization, with
precisely the primary's stated scope.
