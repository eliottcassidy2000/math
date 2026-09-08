# Independent audit: the fourfold DG unit response has exact primary order three

**Verdict: ACCEPTED, with the stated affine response-module scope.**
The new corollary in `incoming_synthesis.md`, section 4, is correct for every
`a,h in C*` and arbitrary `B0 in C`. Its polynomial witness gives the upper
order; its two unequal component principal parts give the all-degree lower
order. The canonical connection raises the order by exactly one at every
step. These conclusions do not assert anything about the response quotient
formed in `O(W)`, a moving-source response, or JC(2).

Audited input: `C:/w/continuing10_20260908_incoming/incoming_synthesis.md`,
section 4, and `fourfold_unit_bridge.py/.out`. The referee used no producer
imports and did not infer the all-degree conclusion from its finite checks.
The independent executable is `04-computation/continuing10_20260908_fourfold_unit_bridge_audit.py`.

## 1. Supplier definitions and paid hypotheses

I read the actual complete torsion and connection proof in
`05-knowledge/results/planar_jc48_sep06_torsion.md`, sections 2--4, including
the pole lemma, scalar-coefficient recursion, diagonal kernel, CRT
surjectivity, and connection formula. Its antecedents are **THM-3770**,
`01-canon/theorems/THM-3770-vertical-principal-part-equalizer-and-log-canonical-dressing-gate.md`,
and **THM-3412**,
`01-canon/theorems/THM-3412-hamiltonian-principal-part-differential-and-prufer-torsion-arms.md`.
The supplier itself credits classical relative-exactness antecedents; this
audit makes no literature-priority claim.

I also checked the actual source family, two charts and fourfold special
fibre in `continuing10_20260907_dg_linear_carrier.md`, sections 2, 4 and 6,
and compared `planar_jc48_sep08_unit_torsion.md`, sections 1--4. The latter
is an order-two comparison family, not a substitute proof for this order.

Let `z=x-h`, `c=B0-3ah^2`, `K=z^3t+z-2h`. Then

    F=c+g,  g=azK,
    D=F_z partial_t-F_t partial_z,
    C_F=C[z,t]/D(C[z,t]),  theta=[1].

Here `F_t=az^4`, and `F_z|z=0=-2ah` is nonzero. Thus the gradient has
no affine zero, equivalently its ideal is the unit ideal. More explicitly,
the following polynomial Bezout field (with scalar denominators) pays
that hypothesis and the connection input directly:

    V=A partial_z+B partial_t,
    A=-(1+z/h+z^2/h^2+z^3/h^3+2z^3t/h)/(2ah),
    B=(4h^2t^2z^2+4h^2t+2htz+2tz^2+1)/(ah^4),
    V(F)=1.

The field equality and derivation are literal:

    t=((F-c)/a-z(z-2h))/z^4,
    C(z,t)=C(F)(z),
    D=-az^4 partial_z |_F.

In characteristic zero the last field has constant field exactly `C(F)`.
Consequently any two rational primitives of the same response differ by
one rational function `H(F)`, not by independent functions on the special
components. Both hypotheses of the complete torsion theorem are now paid.

## 2. All affine fibre components and all affine poles

At the special value c there are exactly two components, `E1={z=0}` and
`E2={K=0}`. The polynomial K is irreducible: as a linear polynomial in t,
its coefficient `z^3` and constant term `z-2h` are coprime. Its source
coordinate ring is `C[z,z^-1]`. The explicit identity

    1=(z^3t+z-K)/(2h)

shows that the two component ideals are comaximal. Both factors occur once,
so the fibre is reduced and its two components are disjoint.

Every other fibre `F=d`, `d!=c`, is irreducible: in its t-linear
presentation, `az^4` is coprime to `az(z-2h)+c-d`. Therefore c is the only
reducible affine fibre. This checks all components, not only two selected
branches in a potentially larger factorization.

The rational primitive

    q=1/(3az^3),  Dq=1

has precisely one affine polar divisor, E1. It is regular on all other
affine divisors, including E2. No unmentioned horizontal pole or pole on
another special fibre needs repair. The module theorem therefore identifies
the entire torsion submodule with a single principal-part arm
`(C^2/C diagonal) tensor g^-1 C[g^-1]` supported at c. This is a statement
about torsion, not the entire possibly nontorsion quotient C_F.

## 3. Independent upper and lower order proofs

The full scalar principal parts in the actual uniformizer g are

    pp_E1(q)=alpha/g^3+beta/g^2,
    alpha=-8a^2h^3/3 != 0,  beta=-2ah,
    pp_E2(q)=0.

I independently recovered the first two coefficients by the limits
`lim_z->0 g^3q` and `lim_z->0 g^2(q-alpha/g^3)`. The next scalar limit
`lim_z->0 g(q-alpha/g^3-beta/g^2)` is zero. More strongly, simplifying
the whole remainder gives a denominator specializing to `-24ah^3` at
z=0. It is therefore regular along E1 for every admitted parameter,
not just for generic a,h. Thus there is no omitted simple scalar term.

For the upper order, the literal polynomial

    P3=g^3 q=(a^2/3)K^3

satisfies `DP3=g^3`, so `g^3 theta=0`.

For the lower order, suppose a polynomial P satisfied `DP=g^2`. Then
the paid constant field forces `P=g^2q+H(F)`. The rational primitive
`g^2q` has a simple pole at E1, whose scalar coefficient in g is the
nonzero alpha; equivalently its coefficient of z^-1 is `4ah^2/3`.
Cancelling it forces H to have a simple pole at c with coefficient
`-alpha`. A higher-order pole cannot work, since g^2q has no such pole
there. At E2, g has order one while g^2q is regular, so precisely that
same H introduces an uncancelled simple pole. This contradicts P being
polynomial. Poles of H at other target values cannot alter this local
contradiction. Hence `g^2 theta!=0`, proving exact primary order three.

Equivalently, theta maps to

    [e_E1] tensor (alpha g^-3+beta g^-2),

with a nonzero leading coefficient modulo the common diagonal. The
lower bound is an analytic argument with no degree bound on P.

## 4. Canonical derivatives and the global-source distinction

For the paid polynomial Bezout field above, the inherited canonical
connection is `nabla[f]=[Vf+(div V)f]`. It is independent of the chosen
Bezout field and differentiates component scalar parts with respect to g.
Consequently, for every integer j>=0,

    nabla^j theta maps to (-1)^j [e_E1] tensor
      (alpha (j+2)!/(2g^(j+3)) + beta (j+1)!/g^(j+2)).

The leading coefficient is nonzero for every j. Thus its exact primary
order is j+3, and the derivative sequence is linearly independent over C.
The executable additionally checks `D(Vq)=div V` and the whole regular
remainder after differentiating the stated principal part; the analytic
supplier proves the unbounded j statement.

The stronger assertion that F has no critical point on all W is also
correct. On the second chart its expression is

    F=c-3ah^2+4ah^3r-ah^4r^2-ab(1-hr)^4,

so `F_b|r=0=-a`. E2 acquires the boundary point `(r,b)=(0,-3h^2)`;
with w=1/z its global parametrization has

    t=2hw^3-w^2,
    r=w/(1+hw),
    b=-2h^5w^3-7h^4w^2-8h^3w-3h^2,
    q=w^3/(3a).

Hence q is regular at this added point as well. These checks pay the
claimed global-submersion comparison, but the module C_F and its polynomial
Bezout field remain defined on the specified affine plane. In particular,
generic fibres on that affine plane are punctured lines; the generic A1
description in the source note refers to W with the boundary point added.
No change of response ring is silently made.

## 5. Hostile parameter boundary and reproduction

The condition h!=0 is essential. At h=0 one has
`g=az^2(1+z^2t)`, a repeated critical component, and already

    g^2 q=(a/3)z(1+z^2t)^2

is polynomial with derivative g^2. Thus exact order three cannot extend
to that boundary. A source pole of order three is not a primary order
of three when the target parameter vanishes twice. The referee never
applies the smooth-fibre torsion theorem to this singular control.

Run `04-computation/continuing10_20260908_fourfold_unit_bridge_audit.py` normally and with `-O`. The two runs
pass **47 always-active exact gates** and have identical raw LF output.
The checks are independent symbolic reconstruction; the all-degree lower
bound, full component accounting and unbounded derivative order are the
proofs above. No repository files or Git state were changed by this audit.
