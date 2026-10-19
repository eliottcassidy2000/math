# Negative-divisor source flows: cyclic and alternating/symmetric groups

2026-10-19. **PROVED; frozen after independent review.** This completes the exact degree,
individual-coordinate, Galois, and radical-solvability analysis on the
negative-divisor algebraic locus of packet 128. The fixed scalar-time
source equation is treated explicitly, not inferred from a generic
parameter specialization. The general alternating/symmetric argument
uses the classification-dependent theorem of Jones, as in packet 135.
No historical novelty is claimed.

## 1. Parameters and theorem

Let K have characteristic zero. Write the negative algebraic locus as

    a=−b,   b≥2,   m+1=br,   r≥2,
    k=2m−b,   n=m−b=(r−1)b−1≥1,
    H=x^k t^m(1+x²t),   λ∈K*.

Let `(X,T)=exp(λ{−,H})(x,t)` be the source t-adic flow. The birational
parameter of [packet 128](planar_two_exponent_flow_classification_20261019.md)
is

    u=x²t,   E=H=x^(−b)u^m(1+u),   z=x/u^r,
    F=K(x,t)=K(E,z),
    u=1/(Ez^b−1),
    x=z/(Ez^b−1)^r,
    t=(Ez^b−1)^(2r−1)/z².                           (1)

Define

    P_E(v)=−Σ_(j=0)^(r−2) binom(r−2,j)(−1)^(r−2−j)
                      E^j v^(b(j+1)−1)/(b(j+1)−1).

Then

    P'_E(v)=−v^(b−2)(Ev^b−1)^(r−2),
    deg P_E=n,   leading coefficient ℓ=−E^(r−2)/n.   (2)

The actual source parameter `Z=X/(X²T)^r` satisfies

    P_E(Z)=P_E(z)+λ,                                (3)

and X,T are obtained from (1), replacing z by Z. In particular,
`F(X,T)=F(Z)`.

**Theorem.** The polynomial in (3) is irreducible of degree n over F
and over `C(z)`, where `C=overline{K(E)}`. Each coordinate individually
generates its root field:

    F(X)=F(T)=F(Z),   [F(X):F]=[F(T):F]=n.           (4)

The groups and radical cases are as follows.

1. If `r=2`, the geometric splitting-field group is cyclic C_n
   (trivial when n=1). The arithmetic group over F is

       C_n ⋊ Gal(K(ζ_n)/K),                         (5)

   with the usual cyclotomic action. All these coordinates are
   expressible by radicals. The case `b=2` is the rational exception
   `(m,k)=(3,4)` of packet 128.

2. If `r≥3`, the geometric splitting-field group is A_n exactly when
   b and r are both even, and is S_n otherwise. In a geometric A_n
   case put

       δ=(−1)^(b/2−1)n.

   The arithmetic group is A_n if δ is a square in K, and S_n otherwise.
   In a geometric S_n case the arithmetic group is always S_n.

3. Within `r≥3`, the coordinates are expressible by radicals only in
   the cubic case `(b,r,m,k)=(2,3,5,8)`. Every other case has `n≥5`
   and each coordinate is not expressible by radicals.

The geometric convention fixes the algebraic closure of K(E); E and
z have not been numerically specialized.

## 2. The pure-power case r=2

Now `n=b−1` and

    P_E(v)=−v^n/n,
    Z^n=z^n−nλ.                                    (6)

For n≥2 the right side has n simple zeros over an algebraic closure.
The polynomial `V^n−(z^n−nλ)` is Eisenstein at any such zero, viewed
as a valuation of `C(z)`. Thus it is irreducible of degree n there and
over F. The case n=1 is immediate.

Over C all n-th roots of unity are present, so the splitting field is
`C(z,Z)`, with cyclic automorphisms `Z↦ζ_n^j Z`. Over F put
`F_n=K(ζ_n)(E,z)`. The same irreducibility argument gives degree n
over F_n, and the splitting field is `N=F_n(Z)`. Hence

    [N:F]=n [K(ζ_n):K].

For `σ∈Gal(K(ζ_n)/K)`, apply σ to the coefficients in the presentation
`F_n[V]/(V^n−z^n+nλ)` and fix the class of V. This gives an automorphism
of N fixing E,z,Z and supplies a section of the constant-field quotient.
Conjugation acts on `Z↦ζ_n Z` by `ζ_n↦σ(ζ_n)`. This proves exactly the
semidirect product (5), without assuming a splitting of that quotient.
The cyclotomic factor is abelian, so the group is solvable; (6) also
provides an explicit radical expression.

The individual-coordinate claim is not deduced from primitivity here,
since C_n is imprimitive when n is composite. Instead use

    χ(V)=V/(EV^b−1)²,
    θ(V)=(EV^b−1)³/V².

Z is transcendental over C, since z is algebraic over `C(Z)`. Thus any
equality `χ(ζZ)=χ(Z)` or `θ(ζZ)=θ(Z)` is an equality of rational
functions in a transcendental variable. The pole set of χ at
`EV^b=1`, or the zero set of θ at the same points, must then be
preserved by multiplication by ζ. This forces `ζ^b=1`. But also
`ζ^n=1` and `gcd(b,n)=gcd(b,b−1)=1`, so ζ=1. Each coordinate has
trivial stabilizer in C_n and degree n over `C(z)`. Degree bounds
then give (4) over F. This includes all composite n, and n=1 is again
immediate.

## 3. Distinct critical values when r≥3

Choose `ρ^b=1/E`. The b nonzero critical points `ρζ`, `ζ^b=1`, each
have derivative multiplicity `r−2`. Zero is an additional critical
point of multiplicity `b−2` when b≥3; for b=2 it is not critical.
The values at the nonzero critical points are

    P_E(ρζ)=κρ^(b−1)ζ⁻¹,
    κ=(−1)^(r−1)b^(r−2)(r−2)!
                          /∏_(j=1)^(r−1)(bj−1) ≠0. (7)

This follows by substituting into the primitive. The remaining finite
sum is `−∫_0^1 v^(b−2)(v^b−1)^(r−2)dv`, whose elementary beta-integral
value is the displayed κ. The values are pairwise distinct and nonzero;
when zero is critical its value is zero. Consequently the polynomial
has at least two distinct critical points and all distinct critical
values.

In particular P_E is indecomposable, by the elementary argument in
[packet 135, section 3](planar_positive_divisor_galois_classification_20261019.md).
Briefly, in a hypothetical polynomial composition every point over a
critical point of the outer factor would be critical with the same
image. Critical injectivity forces the inner factor to be a translated
pure power. There can be only one outer critical point, which then
makes the whole polynomial a pure power, contradicting its multiple
critical points. The block/intermediate-field correspondence and
Lüroth's theorem therefore make the generic-target group primitive.
These standard equivalences and the critical-injectivity criterion
also appear in [Kreso–Tichy, Lemma 3.4 and Proposition 3.9](https://arxiv.org/pdf/1601.07316).

## 4. Generic-target alternating and symmetric groups

For `P_E(V)−Y` over `C(Y)`, formal lifting at the distinct critical
values gives a single `(r−1)`-cycle at each nonzero critical point.
When b≥3, zero gives a single `(b−1)`-cycle. All other roots are
fixed by the corresponding local inertia. The lifting justification
is the same coefficientwise construction as in packet 135: an index-e
critical point is resolved by the parameter substitution `Y−Y_0=η^e`,
while all unramified roots lift in the original parameter.

If b≥3, the `(b−1)`-cycle fixes

    n−(b−1)=(r−2)b≥3

points. If b=2 and r≥5, an `(r−1)`-cycle fixes `r−2≥3` points. In
both cases the primitive group contains A_n by [Jones, Corollary 1.3](https://www.cambridge.org/core/services/aop-cambridge-core/content/view/BBF3021D8434AD2C0BB0D3D6CE637008/S000497271300049Xa.pdf/primitive-permutation-groups-containing-a-cycle.pdf):
a primitive finite permutation group containing a single cycle with
at least three fixed points contains the alternating group. The
cycle length need not be prime. The general theorem uses the
classification of finite simple groups, which is an explicit dependency
of this argument.

The two remaining pairs are `(b,r)=(2,3)` and `(2,4)`. The former has
degree three and a transposition, hence group S3. The latter has degree
five and a single 3-cycle. A primitive group containing a 3-cycle
contains the alternating group; this classical small-cycle criterion
is stated in [Kreso–Tichy, Remark 3.11](https://arxiv.org/pdf/1601.07316).
Thus it contains A5 as well.

If b or r is odd, one of the local cycles is odd: its length is b−1
or r−1 respectively. Therefore the group is S_n. If both are even,
the monic derivative is a constant times a square polynomial, and
the discriminant is square over C(Y), as computed in section 6.
The group is then A_n. This covers all `r≥3`, including the two
small cases.

## 5. Transfer to every actual fixed nonzero scalar time

Let B be the finite critical-value set from (7), including zero only
if b≥3. Any nonzero difference of two members of B has the form

    κρ^(b−1)(ζ_i−ζ_j),

where each ζ is a b-th root of unity, or zero when applicable. It
cannot equal the fixed `λ∈K*`: taking b-th powers would make
`E^(b−1)` algebraic over K. Therefore

    B ∩ (B+λ)=∅.                                    (8)

Let L and L_λ be the splitting fields over `C(Y)` of `P_E(V)−Y`
and `P_E(U)+λ−Y`. Their common Galois subextension can ramify only
at infinity, since their finite branch sets are B and `B+λ`.
Tame Riemann–Hurwitz forces any such subextension to be trivial:
for degree d and h points above infinity,

    2g−2=−2d+(d−h)=−d−h,

and `g≥0,h≥1` imply `d=h=1`. The precise characteristic-zero formula
is given in the [Stacks Project, Section 53.12](https://stacks.math.columbia.edu/tag/0C1B).
Hence `L∩L_λ=C(Y)` and the two Galois extensions are linearly disjoint.

Embed `C(z)` into L_λ by `Y=P_E(z)+λ`. Then

    Gal(L C(z)/C(z))=Gal(L/C(Y)),

and the left side is the splitting-field group of (3). Thus both the
entire permutation group and irreducibility are preserved. In
particular the actual source branch has degree n over `C(z)` and F.
This is not an appeal to an unspecified generic-specialization theorem.

At λ=0, the source equation has the factor `V−z`; times depending on E
may also make (8) fail. Neither is in the claimed scope.

## 6. Arithmetic square class for r≥3

If the geometric group is S_n, arithmetic equality with S_n follows
immediately. Suppose b and r are both even, so n is odd. With

    f_Y(V)=[P_E(V)−Y]/ℓ,
    h(V)=V^((b−2)/2)(EV^b−1)^((r−2)/2),

we have `f'_Y=(n/E^(r−2))h²`. Thus

    Disc(f_Y)=(−1)^(n(n−1)/2)
                (n/E^(r−2))^n Res(f_Y,h)².

As `r−2` is even, its square class is

    (−1)^((n−1)/2)n=(−1)^(b/2−1)n=δ.              (9)

After substitution `Y=P_E(z)+λ`, the square factor remains nonzero
by separability. The arithmetic group contains the geometric A_n and
is contained in S_n; the discriminant sign character gives A_n exactly
when δ is square in `K(E,z)`. Since this is a purely transcendental
extension of K, that is equivalent to δ being square in K.

For example `(b,r)=(2,4)` gives geometric A5 and δ=5: over Q the
arithmetic group is S5, and over `Q(√5)` it is A5. Already
`(b,r,m,k)=(2,6,11,20)` gives `n=9,δ=9`, so both groups are A9 over Q.

## 7. Individual coordinate generation for r≥3

In the generic root field `C(v)/C(P_E(v))`, consider

    χ(v)=v/(Ev^b−1)^r,
    θ(v)=(Ev^b−1)^(2r−1)/v².

Neither function belongs to `C(P_E(v))`.

For χ, at a nonzero critical point c with `Ec^b=1`, its pole has
order r. If `χ=R∘P_E`, that pole order would be a multiple of the
local ramification index `r−1`. Since `r≥3`, r is not divisible by
r−1, a contradiction.

For θ, zero is a pole but no other finite point is a pole. Yet
`P_E(0)=0`, and the order of zero in P_E is b−1, strictly less than
`deg P_E=n` because `r≥3`. Hence `P_E(v)=0` has some nonzero finite
root as well. If `θ=R∘P_E`, the pole of R at zero required by θ's
pole at v=0 would produce a pole at that other root, a contradiction.

The linear disjointness from section 5 gives `L∩C(z)=C(Y)`. Thus
neither χ(v) nor θ(v) can become rational over the source field
`C(z)` after the fixed-time substitution. The natural A_n/S_n action
is primitive (also for the S3 case), so its root field has no proper
intermediate fields. Each nonrational coordinate generates the whole
degree-n field. Degree bounds descend this conclusion to F, proving
(4). This is separate from the direct cyclic stabilizer argument
needed in section 2.

## 8. Radicals and completion of the algebraic-locus picture

All r=2 cases are radical by (6), and their arithmetic groups are the
solvable semidirect products (5). For `r≥3`, the only degree below five
is `(b,r)=(2,3)`, with `n=3` and group S3. It is radical-solvable.
All other r≥3 cases have degree at least five and arithmetic group
A_n or S_n, which is nonsolvable. Each coordinate generates the full
root field, so its normal closure is that same splitting field;
therefore neither coordinate is expressible by radicals. The standard
criterion is recalled in [Sharifi, *Abstract Algebra*, Section 10.3](https://math.ucla.edu/~sharifi/notes/algebra-ch10.html).

Together with packet 135 and the rational a=0 case of packet 128, this
settles exact individual degrees and radical solvability across that
entire algebraic locus. The radical cases are precisely the rational
a=0 locus, the positive quartic `(m,k)=(2,5)`, every negative r=2
Kummer case, and the negative cubic `(m,k)=(5,8)`.

This does not classify numerical specializations, invariant-dependent
times, or other Hamiltonian families. It does not alter the previously
proved polynomiality obstruction and gives no planar-JC conclusion.

## Verification scope

The [small exact checker](../../04-computation/experiments/planar_negative_divisor_galois_classification_20261019.py)
passed 41 primitive, coordinate, Kummer, and discriminant checks on five
specified witnesses. Its [normal output](planar_negative_divisor_galois_classification_20261019.json)
and [optimized output](planar_negative_divisor_galois_classification_20261019_optimized.json)
are byte-identical. These are identity witnesses, not a Galois-group
computation or a substitute for the universal proof.

The [independent proof review](planar_negative_divisor_galois_classification_20261019_independent_review.md)
accepted the full mathematics without corrections. The final changes
after its pinned draft concern only review and verification metadata.
