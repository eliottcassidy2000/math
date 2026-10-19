# Independent universal-proof review of the two-exponent classification

2026-10-19. **ACCEPTED.** Reviewed the full producer note
`planar_two_exponent_flow_classification_20261019.md`, SHA-256
`3f0500d2c3321feac096ae0f427f5c3906e6e4d05789d3595f327c8669499b9c`.

The final revision adds artifact links and the producer's control count;
the reviewed mathematical argument is unchanged. This review does not
independently certify that control count.

The algebraicity locus and rational-map locus are correct for
`H=x^k t^m(1+x²t)`, `m≥2`, `k≥0`, characteristic zero, and every nonzero
scalar time. This review checks the universal mathematical proof, including
source-completion issues and specialization. It does not replay the
companion controls or claim additional irreducibility, individual degrees,
or coordinatewise nonrationality in the algebraic cases.

## 1. Source and scaling hypotheses

Put `u=x²t`, `a=k−2m`, and `E=x^a u^m(1+u)`. Direct differentiation gives

    δu=−aEx,
    δx=x^(a+2)u^(m−1)[m+(m+1)u].

The latter expression is a rational-coordinate expression; the actual
source derivation is polynomial for every `k≥0`. It raises t-adic order
by at least `m−1≥1`, so all scalar-time source exponentials exist.

The invariant equation is genuinely a positive-degree equation in either
coordinate when the other is algebraic. Even at `k=0`, its equation in X
has degree two with nonzero leading coefficient `T^(m+1)`. Consequently
exclusion of simultaneous algebraicity proves individual transcendence.

The dilation weight of δ is `a+1`. For `a≠−1`, the fixed-scalar to
formal-time argument from packet 127 therefore applies without alteration
except for this weight. The common coefficient-array ring is

    K[x,τ][[t]] ∩ K[x,t][[τ]].

The substitution `τ=λz^(a+1)` is injective coefficientwise in t for
either sign of `a+1`, and the extension in z has finite degree `|a+1|`.
Cleared algebraic relations therefore descend and regroup legitimately.

## 2. The finite cover and exactness test

For `a≠0`, extend the initial-value constant field to an algebraic closure
C and choose `η^a=E`. In the formal-time flow, set `v=U_τ,w=η/X_τ`.
Then

    w^a=v^m(1+v),   dτ=−w dv/(aηE).

The nonzero derivative `dv/dτ=−aEX_τ` makes v transcendental over C.
The simple valuation at `v=−1` makes the curve polynomial irreducible of
degree `|a|`; for negative a, invert w. The algebraic flow field and the
curve field have transcendence degree one, and the former is finitely
generated over the latter. Thus the relevant cover is finite separable.
Normalized trace of τ gives a rational primitive of `w dv`. The enlarged
constant field causes no descent difficulty, since the obstruction is
proved over C itself.

Projection to the w-character of the cyclic curve extension puts a
primitive in the form `wR(v)`. This is valid also when `|a|=1`, when the
character group is trivial. Its equation is exactly

    a v(1+v)R′+[m+(m+1)v]R=a v(1+v).

At a finite pole away from 0 and −1, differentiation creates an unmatched
higher-order pole. A pole of order r at 0 requires `ar=m`; a pole at −1
requires `ar=1`.

Thus a positive a not dividing m forces R to be polynomial; its leading
degree and then the degree-one coefficient comparison exclude it. For
`a=−b<0`, every rational solution is polynomial. A degree at least two
requires `bd=m+1`; degree at most one is impossible. This requires b to
be a proper divisor of `m+1`. At `b=1`, the residue of
`dv/[v^m(1+v)]` at −1 is `(-1)^m`, so exactness fails. The only remaining
negative candidates are therefore precisely the stated `b≥2` proper
divisors. This is an all-integer proof, not an inference from controls.

## 3. The resonant scalar-time case

At `a=−1`, the restricted one-variable flow is

    U=exp(λu^m(1+u)∂_u)(u)
      =u+λu^m+O(u^(m+1)).

The source embedding `t=u/x²` puts this U in `K[[u]]` and makes
`E=u^m(1+u)/x` transcendental over `K((u))`. Expanding a cleared relation
as a polynomial in E therefore correctly descends algebraicity of U from
`K(u,E)` to `K(u)`.

The draft's partial fractions have the correct sign:

    du/[u^m(1+u)]
      =dB(u)+(-1)^m dlog((1+u)/u).

For example,

    B(u)=−Σ_(j=2)^m (-1)^(m−j)/[(j−1)u^(j−1)].

Preservation of the differential gives `dh=dr/r` with
`h=−[B(U)−B(u)]/(-1)^m` and
`r=(1+U)u/[U(1+u)]`. In an algebraic function field, a nonconstant h
would permit the norm to the rational field in h, giving a rational
logarithmic derivative equal to a positive integer. Its behavior at
infinity makes this impossible. Then h and r are constant, and the source
constant term forces `r=1`, hence `U=u`, contradicting `λ≠0`. The general
resonant exclusion is therefore valid at every nonzero scalar time.

## 4. Constructive algebraicity and scalar specialization

The two claimed birational parameters and time primitives were derived
independently.

For positive `a|m`, set `r=m/a` and `s=1/(xu^r)`. Then

    u=Es^a−1,
    x=1/[s(Es^a−1)^r],
    t=s²(Es^a−1)^(2r+1),
    δs=−1/[s^a(Es^a−1)^r].

Thus the primitive

    P_E(s)=−∫s^a(Es^a−1)^r ds

is polynomial of degree `m+a+1`, with `δP_E(s)=1`.

For `a=−b`, `b≥2`, and `r=(m+1)/b≥2`, set `z=x/u^r`. Then

    u=1/(Ez^b−1),
    x=z/(Ez^b−1)^r,
    t=(Ez^b−1)^(2r−1)/z²,
    δz=−1/[z^(b−2)(Ez^b−1)^(r−2)].

The corresponding polynomial

    P_E(z)=−∫z^(b−2)(Ez^b−1)^(r−2) dz

has degree `m−b` and derivative one along δ. In both cases the inverse
formulas establish the required full rational function field equality.

The source exponential extends to `K(x)((t))`: the same valuation estimate
works for negative Laurent orders and rational x-coefficients. It is an
algebra automorphism there. Applying it to either primitive therefore
gives the exact identity `P_E(C_λ(v))=P_E(v)+λ` on the existing source
branch, at every scalar time. Since the displayed polynomial has nonzero
leading coefficient independent of λ, it proves algebraicity directly.

The draft's general specialization safeguard is also correct. After
clearing denominators, cancel the maximal common power of `τ−λ` from
the polynomial relation. Cancellation of its evaluated relation is valid
in the coefficient-array domain. Evaluation at τ=λ is defined
coefficientwise in the t-adic arrangement, and the remaining specialized
polynomial is nonzero. Thus a formal-time relation cannot fail to prove
scalar-time algebraicity merely because its original coefficients share
a vanishing time factor.

## 5. Rational-map classification

The lemma excluding `P(R)=P+λ` for `deg P≥2` and rational R is sound.
A finite pole of R over an algebraic closure would persist under P, so
R is polynomial. Degree comparison makes it affine, `R=αv+β`, with
`α^d=1`. If `α=1`, the next-to-leading coefficient forces `β=0`. If
`α≠1`, R has finite order, and iteration contradicts nonzero λ in
characteristic zero.

Apply this over `K(E)`. A rational source coordinate pair makes its
birational parameter image rational, so the lemma excludes rational maps
whenever the primitive degree is at least two. The positive primitive
always has degree at least four. In the negative family, degree one
requires

    m−b=1,   m+1=br,   b,r≥2,

or `b(r−1)=2`. The unique solution is `b=2,r=2,m=3,k=4`.

At this exceptional pair, `z=1/(x³t²)` satisfies `δz=−1`. Substituting
`z−λ` in the inverse formulas gives exactly

    D=1−λx³t²,   B=(1+x²t)D²−x²t,
    X=xD/B²,   T=tB³/D².

The draft correctly proves T nonpolynomial using `gcd(D,tB)=1`.
The `a=0` geometric-series formula is also correct and has nonpolynomial
first coordinate for every nonzero λ.

## Disposition and scope

No mathematical correction is required in the reviewed draft. It proves:

- Algebraicity exactly at `a=0`, positive divisors `a|m`, and negatives
  `a=−b` with `b≥2` a proper divisor of `m+1`.
- Individual transcendence outside that locus.
- Rationality of the pair exactly at `a=0` and `(m,k)=(3,4)`.
- Nonpolynomiality of these nonzero-time rational maps.

The primitive degrees are not asserted to be scalar algebraic extension
degrees, and the nonrational-pair argument is not promoted to separate
coordinatewise nonrationality. These limits are necessary and are
respected by the note. All frozen packets reviewed earlier were preserved;
this review performed no Git operation or companion-checker replay.
