# Independent review of the fixed-time and trace bridges

2026-10-19. Reviewed independently against the packet 127 draft
`planar_monomial_carrier_flow_trichotomy_20261019.md`, particularly sections 1–5.

Reviewed final draft SHA-256:
`98ba600bf22b72f62f605952adfe8a2521be6b2a0c9970ff959dbf8288cad26a`.

**Verdict: accepted.** The argument from algebraicity of a nonzero scalar-time
source flow to exactness of the invariant-curve time differential is valid.
The rational-function exclusion and the separate resonant case `k=3` are
also valid. Together they prove individual transcendence for every integer
`k≥0` outside `{4,5,6}`. No computational test or external classification of
binomial integrals is needed for this conclusion.

This review checks the universal exclusion proof. It does not independently
re-run the companion checker or certify the exceptional coordinate degrees
in sections 6–7. Those are separate constructive claims.

## 1. The scalar source exponential is well-defined

Use the bracket `{f,g}=f_x g_t-f_t g_x`. For

    H=x^k t²(1+x²t),   δ={−,H},

one has `δx∈tK[x,t]` and `δt∈t²K[x,t]`. Therefore δ raises the t-adic order
of every polynomial or t-adic series by at least one, apart from a zero
derivative. In particular,

    ord_t(δ^n x)≥n,   ord_t(δ^n t)≥n+1.

The scalar exponential exists coefficientwise in `K[x][[t]]`. The formal-time
exponential is compatible with the same coefficient array and has its
τ-coefficients in `K[x,t]`. The case `k=0` causes no negative powers of x in
the actual derivation: the apparent term `k x^(k−1)` is zero there.

Invariance gives `H(X,T)=H(x,t)`. If one coordinate is algebraic, the other
satisfies a nonzero positive-degree polynomial over the corresponding
algebraic extension of `K(x,t)`. This applies to both coordinates for every
`k≥0`; when solving for X at `k=0`, use

    T³X²+T²−E=0.

The coordinate series are nonzero. Thus exclusion of simultaneous
algebraicity really does establish individual transcendence.

## 2. Fixed-time algebraicity implies formal-time algebraicity

Assign dilation weights `wt(x)=1` and `wt(t)=−2`. Both monomials of H have
weight `a=k−4`, and δ increases the weight by `b=a+1=k−3`. Consequently

    ρ_z(δ^n x)=z^(1+nb) δ^n x,
    ρ_z(δ^n t)=z^(−2+nb) δ^n t,

where `ρ_z(x,t)=(zx,z^(−2)t)`. The two displayed dilation identities in the
draft therefore hold coefficientwise, not merely heuristically as an
analytic conjugation statement.

Here is an explicit ring version of the completion argument. Inside
`K[x][[t,τ]]`, let

    A=K[x,τ][[t]] ∩ K[x,t][[τ]].

Both formal-time coordinate series lie in A. For `λ≠0` and `b≠0`, the map

    K[x,τ][[t]] → K[x,z,z^(−1)][[t]],
    τ ↦ λz^b,

is injective: a nonzero polynomial in τ remains nonzero after substituting
`λz^b`, for either sign of b, and one applies this separately to every
t-coefficient. This explicitly handles the negative dilation weights
`k=0,1,2`.

Start from a polynomial relation over `K(x,t)` for a fixed scalar-time
coordinate. Dilation gives algebraicity of the corresponding coordinate at
time `λz^b` over `K(x,t,z)`. Its subfield

    K(x,t,τ),   τ=λz^b,

has finite index in `K(x,t,z)`; the degree is `|b|`. Transitivity of
algebraicity supplies a nonzero relation over `K(x,t,τ)`. Clearing
denominators gives `Q∈K[x,t,τ,Y]`, and evaluating Q on the formal-time
coordinate belongs to A. Its image under the injective substitution above
is zero, so it is zero in A. It is therefore also zero in
`K[x,t][[τ]]`.

No specialization at a value where a denominator vanishes is used. No
unjustified continuous map from a t-adic completion to a τ-adic completion
is used. The argument compares the same polynomial coefficient arrays.

## 3. Constant extension and the actual finite cover

Let `F=K(x₀,t₀)` be the field of initial values, and regard all of F as
constants for formal-time differentiation. Take an algebraic closure C of
F. Formal algebraicity over `F(τ)` remains algebraicity over `C(τ)`, so

    L=C(τ,X(τ),T(τ)) ⊂ C((τ))

is a finite extension of `C(τ)`.

Assume `a≠0`, put `E=H(x₀,t₀)∈C*`, choose `η∈C*` satisfying `η^a=E`, and
define

    v=X²T,   w=η/X.

Then `w^a=v²(1+v)` and

    dv/dτ=−aEX=−aEη/w.

This derivative is nonzero. Hence v is a nonconstant formal series and is
transcendental over the algebraically closed field C. In particular the
curve map is not constant, even though the initial point belongs to C.

For positive a, `W^a−v²(1+v)` is Eisenstein at `v+1` over `C[v]`, or over
the corresponding discrete valuation ring. For negative a, setting
`s=w^(−1)` gives `s^(−a)=v²(1+v)` and the same Eisenstein argument. Thus

    M=C(v,w)

is precisely the indicated irreducible curve function field, of degree
`|a|` over `C(v)`. The geometric irreducibility survives the constant
extension because the simple zero at `v=−1` survives it.

Both M and L have transcendence degree one over C. Since L is finitely
generated over C and contains M, the extension `L/M` is algebraic and
finitely generated, hence finite. Characteristic zero makes it separable.
These points establish the hypotheses needed for trace; finiteness over
`C(τ)` alone would not have been enough without a nonconstant curve map.

## 4. Trace gives an exact differential, not merely a formal primitive

As an equality in the space of differentials of L over C,

    dτ=−w dv/(aηE).

For a finite separable extension `L/M`, the differential trace satisfies

    Tr(dα)=d Tr(α),
    Tr(ω)=[L:M]ω   for ω∈Ω_(M/C).

For an elementary verification here, choose the separating variable v and
extend the derivation `d/dv` uniquely to the finite separable fields. It
commutes with all embeddings into an algebraic closure, and hence with
their sum, the trace. Applying this to τ gives

    d Tr_(L/M)(τ)=−[L:M] w dv/(aηE).

Thus

    B=−aηE Tr_(L/M)(τ)/[L:M] ∈ M

satisfies `dB=w dv`. Division by the degree is legitimate in characteristic
zero. No descent of B to F or K is necessary: the subsequent impossibility
proof is valid over the larger C itself.

## 5. Character projection and the rational-function obstruction

Let `n=|a|`. The cyclic Galois group of `M/C(v)` sends `w↦ζw` for
`ζ^n=1`. Its eigenspace of this character is `wC(v)`. This includes `n=1`,
when `w∈C(v)*` and `wC(v)=C(v)`. Projecting a primitive B into this
eigenspace is valid because differentiation commutes with the group
action. Therefore exactness implies a primitive `wR(v)` with `R∈C(v)`.
Differentiation gives exactly

    a v(1+v)R′+(2+3v)R=a v(1+v).

For `a≥3` or `a<0`, a rational R has no finite poles:

- At a point other than 0 and −1, differentiation produces an unmatched
  pole one order higher.
- A pole of order r at 0 requires `ar=2`.
- A pole of order r at −1 requires `ar=1`.

The last two equations are impossible in the stated ranges, so R is a
polynomial. A polynomial of degree `d≥2` requires cancellation of its
highest term, namely `ad+3=0`. For integral a in these ranges the only
possibility is `a=−1,d=3`. Otherwise, substituting `R=Av+B` gives

    2B=0,
    (a+2)A+3B=a,
    (a+3)A=a,

which forces `A=B=0` and then `a=0`, a contradiction. For `a=−1`, one can
use the draft's explicit cubic comparison, or simply observe that

    w dv=dv/[v²(1+v)]

has residue −1 at 0 and therefore cannot be a derivative of a rational
function. This verifies the exceptional resonance as well.

For reference, rational solutions exist at the two omitted values:

    a=1: R=v(4+3v)/(12(1+v)),
    a=2: R=2(1+v)(3v−2)/(15v).

These confirm the exactness test's exceptional set, although sufficiency
for the claimed fixed-time coordinate degrees is established separately.

## 6. The resonant scalar-time case k=3

The dilation has zero time weight here, so the previous scaling lemma
cannot be applied. The separate argument in the draft correctly avoids
that gap.

The source substitution `t=u/x²` is injective into `K(x)[[u]]`, and the
restricted flow has

    U=exp(λu²(1+u)∂_u)(u) ∈ K[[u]],
    U=u+λu²+O(u³).

Within `K(x)((u))`, x is transcendental over `K((u))`: comparing the
u-coefficients of a putative finite polynomial relation in x proves this
directly. Accordingly `E=u²(1+u)/x` is also transcendental over `K((u))`.
If U were algebraic over `K(u,E)`, a nonzero cleared relation could be
written

    Σ_j E^j P_j(u,U)=0,   P_j∈K[u,Y].

Each evaluated coefficient belongs to `K((u))`. Transcendence of E forces
all evaluated coefficients to vanish. At least one `P_j` is a nonzero
polynomial, proving algebraicity over `K(u)`. Thus the descent is justified
without an appeal to an unproved general constant-field principle.

Preservation of the vector field gives

    U′u²(1+u)=U²(1+U).

Set `h=1/U−1/u` and `r=(1+U)u/[U(1+u)]`. Partial fractions then give
`dh=dr/r`, with the sign as written in the draft.

Assume U algebraic and extend constants to an algebraic closure C of K.
In `N=C(u,U)`, if h is nonconstant, `N/C(h)` is finite separable. Extending
`d/dh` and taking the norm of r produces

    (Norm(r))′/Norm(r)=[N:C(h)].

The left side is a rational logarithmic derivative in h, hence is
`O(1/h)` at infinity, while the right side is a nonzero integer. This is
impossible. Equivalently, on the smooth complete curve of N a pole of h
would give a pole of dh of order at least two, whereas `dr/r` has only
simple poles. Therefore h is constant. Then r is constant. Since its
formal expansion has `r(0)=1`, one has `r=1`, which directly forces `U=u`.
This contradicts `λ≠0` in the displayed expansion.

Finally, algebraicity of either source coordinate would imply algebraicity
of both and hence of U over `K(u,E)`. The contradiction therefore proves
individual transcendence of both source coordinates at `k=3`.

## 7. Review disposition

No mathematical correction is required in sections 1–5 of the reviewed
final draft. The explicit injective substitution recommended during review
is present in its section 2. Using C consistently after extending constants
in the section 5 norm proof remains an optional notation improvement, not
a proof gap. The constant extension, irreducibility, nonconstant map,
finite-separability, and coefficientwise-regrouping hypotheses have all
been checked.
