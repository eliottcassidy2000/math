# Prime fibre degree forces one leading root on every W_m

**Status: PROVED DERIVED FILTER + FINITE-EXACT + INDEPENDENTLY AUDITED.**
This is a necessary leading-coefficient filter, derived from proved
mechanisms. It gives no full quintic exclusion, no prime-degree root
gluing, and no reduction of arbitrary Keller maps to these surfaces.

## 1. All graph degrees, all prime degrees, unrestricted polynomial mates

For every integer `m>=1` use the actual surface and charts

    W_m=(P1_x x P1_Z)\{Z=x^m}, t=1/(Z-x^m),
    x=1/r, t=-r^m-r^(2m)b,
    omega=dx wedge dt=r^(2m-2) dr wedge db.

Let `p>=2` be any prime. Suppose `F in O(W_m)` has degree exactly
`p` in `t` on the original affine chart and has a **polynomial** mate
`G in C[x,t]` with `J(F,G)=lambda!=0`. There is no degree bound on
`G`, and it is not required to extend to `W_m`. Let `A(x)!=0` be
the leading `t` coefficient of `F`.

**Theorem.** Either `A` is a nonzero constant, or

    A=c(x-h)^(kp), c!=0, h in C, 2<=k<=2m.           (1)

Every finite location `h` is retained; it is not moved by an assumed
surface translation. In the nonconstant case the leading boundary
section has the single finite zero of multiplicity `kp` at `h` and
original infinity multiplicity

    p(2m-k).                                        (2)

For constant `A`, there is no finite zero and the original infinity
multiplicity is `2mp`. In particular, for `W_2` and quintic fibre degree,
the only leading patterns left by this necessary gate are

| k | Leading coefficient | Finite multiplicity | Original infinity multiplicity |
| --- | --- | --- | --- |
| 0 | c | none | 20 |
| 2 | c(x-h)^10 | 10 | 10 |
| 3 | c(x-h)^15 | 15 | 5 |
| 4 | c(x-h)^20 | 20 | none |

These are open leading strata for the next investigation, not asserted
Keller examples. Lower coefficients, global approximate roots and full
mate obstructions remain additional work.

The closest inherited mechanisms are the leading-degree target-shear
argument used in [quartic boundary, Section 2](planar_jc48_sep08_quartic_boundary.md),
the general [radical leading exactness theorem](planar_jc48_sep08_leading_exactness.md),
the all-degree reciprocal-polynomial lemma in
[continuing10_20260908_dg_linear_all_m, Section 2](continuing10_20260908_dg_linear_all_m.md),
and the [complete all-m filtration](planar_jc48_sep08_dg_genus.md), Section 2.
Targeted current-canon/result searches for the combined prime-fibre,
single-leading-root statement did not retrieve an identical W_m theorem.
Its components are inherited; this note presents their derived filter
and makes no external priority claim or new prime-degree Keller-rigidity
claim.

The live concepts are polynomial target shears, coprime leading valuations,
the exact coefficient field, reciprocal-polynomial pole degree, and the
marked infinity section. The map retains the original source bracket but
discards all lower coefficients. Its sharp missing-coordinate control is
the actual rational mate in Section 5: polynomiality of the proposed
mate is needed before the leading radical field collapses to `C(x)`.

## 2. Polynomial shears and the prime valuation step

Scale `G` so that `J(F,G)=1`. If its current positive `t` degree is
`n`, with leading coefficient `g_n(x)`, the coefficient at degree
`p+n-1` in its Jacobian is

    n A' g_n-p A g_n'=0.

Consequently

    g_n^p/A^n is a nonzero constant in C.            (3)

This uses the usual constant field of `C(x)` and has no generic-parameter
premise. If `p` divides `n`, write `n=ps`. Equation (3) gives
`g_n=c A^s` for a complex constant `c!=0`; subtract the polynomial
`cF^s` from `G`. Its leading term cancels, its degree decreases, and
its Jacobian remains one. The current mate stays polynomial throughout.

A degree-zero mate is impossible, either initially or after these
reductions:
for `G=B(x)`, one has `J(F,G)=-F_t B'`. If `B'!=0`, this has
positive `t` degree `p-1`; if `B'=0`, it is zero. Neither is one.
Thus a finite sequence of constant polynomial shears produces a mate
of positive degree `n` with `p` not dividing `n`.

Since `p` is prime, `gcd(p,n)=1`. At each finite root of `A`, (3)
gives

    p*ord(g_n)=n*ord(A).

Hence every root multiplicity of `A` is divisible by `p`. Factoring
over `C` and absorbing a nonzero scalar into a polynomial root gives

    A=V^p, V in C[x], V!=0.                         (4)

There is no assumption that `G` is global, and no bound on its original
or reduced degree. Primality is used exactly in this coprime valuation
step; the proof does not assume that taking arbitrary powers preserves
a Jacobian equation or a polynomial-mate condition.

## 3. Exactness collapses the possible finite roots of V

Apply the proved leading exactness theorem to the original `F` of
degree `p` and any of its polynomial mates. By (4), a chosen `p`th root
of its leading coefficient is the polynomial `V`. Thus its coefficient
field is exactly `C(x)`, and

    dx/V=dB for some B in C(x).                      (5)

This statement does not require a geometric generic-fibre irreducibility
claim. It follows by formal inversion at `t=infinity` with all rational
denominators and derivations paid in the inherited theorem.

For completeness, the reciprocal-polynomial lemma has a short all-degree
proof. If `V` is nonconstant and has a simple root, `dx/V` has a
nonzero residue there, contradicting (5). Otherwise let `d=deg V`
and let its `s` distinct roots have multiplicities `e_i>=2`. A primitive
has a pole of order exactly `e_i-1` at each root and no other finite
poles. It has no infinity pole, since

    B=B_infinity+gamma x^(1-d)+..., gamma!=0.

Hence its total pole degree is `d-s`, but `B-B_infinity` has a zero
of order `d-1` at infinity. Equality of the total zero and pole degrees
forces `d-1<=d-s`, so `s<=1`. Since `V` is nonconstant, `s=1`.
Therefore

    V=a(x-h)^k, a!=0, k>=2,

or `V` is a nonzero constant. Conversely these reciprocal differentials
are exact, with primitives `x/a` in the constant case and
`(x-h)^(1-k)/(a(1-k))` in the pure-power case. This converse is only
about the differential, not a Jacobian mate of the original `F`.

The single-root lemma is precisely the inherited all-degree result from
the source-linear all-m note. We use only its rational integration
statement, not its additional global-submersion classification.

## 4. Globality supplies the uniform exponent bound and the infinity label

The full all-m filtration writes every global function of `t` degree
at most `p` uniquely as

    F=Q(x,Z)/(Z-x^m)^p,
    deg_x Q<=mp, deg_Z Q<=p.                         (6)

In particular its leading coefficient is `A(x)=Q(x,x^m)`, so

    deg A<=2mp.

Together with `A=V^p`, this gives `deg V<=2m`. The pure-power exponent
in (1) is therefore in the entire integer interval `2<=k<=2m`.
This proves the theorem, including `m=1` and `p=2` as formal instances
of the same necessary statement. Stronger inherited low-degree results
may already exclude those instances; no new low-degree closure is claimed.

The restriction of the numerator in (6) to the degree-m graph is a
binary section of degree `2mp`. At original infinity its local expression
is `r^(2mp)A(1/r)`. Its order is `2mp-deg A`, proving (2) and
the constant-leading alternative. This is an actual labelled infinity
calculation, not an interchange of infinity with the unique finite root.

All displayed leading strata can be realized by global polynomials,
without asserting a mate. Indeed the full basis is

    x^i t^(p-j)(1+x^m t)^j,
    0<=i<=mp, 0<=j<=p,

whose leading monomials are `x^(i+mj)` and span every degree from
zero through `2mp`. Linear combination therefore lifts any polynomial
`A` in that range, including every `c(x-h)^(kp)`. The source verifies
all four quintic W_2 lifts with the arbitrary parameter `h` retained.

## 5. Exact hostiles and the strict stopping scope

The polynomial-mate premise is essential. For every `m>=1` and every
prime `p>=2`, the actual function

    F=x t^p

is global: its second-chart expression is

    (-1)^p r^(mp-1)(1+r^m b)^p.

It has the literal rational mate

    G=-1/[(p-1)t^(p-1)], J(F,G)=1.

But its leading coefficient `x` is not a `p`th power. The pole of `G`
on the original affine divisor `t=0` is genuine. Thus leading exactness
for a rational mate alone does not permit the polynomial shear/valuation
conclusion (4). There is no contradiction: its genuine radical field is
`C(x^(1/p))`, where `dx/x^(1/p)` is exact.

Primality cannot be discarded in the valuation inference either. At
composite degree four and reduced degree two, `A=x^2`, `g=x` satisfy
`g^4=A^2`, although `A` is not a fourth power. The literal pair
`F=x^2t^4+t^2`, `G=xt^2` realizes the leading cancellation but has
Jacobian `-2t^3`. It is a top-identity hostile, **not** a Keller pair or
a counterexample to an unproved composite-degree theorem.

Finally, the leading condition is not sufficient: the global `F=t^5`
has constant leading coefficient but its two partial derivatives vanish
on `t=0`, so it has no polynomial mate. The four W_2 quintic rows in
Section 1 retain the lower coefficients as the precise unpaid work.
No approximate-root globality, classification of quintic mates, or
exclusion of all global pairs on W_m is inferred from this leading filter.

## 6. Reproducible finite controls

The proof is uniform in every integer `m>=1`, every prime `p>=2`, and
every polynomial mate degree. The source's explicitly named controls
use primes `2,3,5,7`; they verify the highest Jacobian identity, constant
ratio differentiation, all nonzero residue classes in each named prime,
polynomial shear invariance, and the reciprocal primitive formulas.
They are not an extrapolated prime or mate-degree census.

The four complete W_2 quintic leading lifts retain arbitrary finite `h`
and are checked on the whole second chart. The rational-mate hostile is
checked in the named primes and graph degrees `1,2,3`, while its displayed
formula proves the all-m/all-prime claim analytically. Composite and
necessary-only controls keep the distinct failed implications separate.

    python3 -B 04-computation/planar_jc48_sep08_prime_leading.py
    python3 -B -O 04-computation/planar_jc48_sep08_prime_leading.py

Both complete producer modes pass **166 always-active gates**, with
identical 528-byte output. Frozen SHA256:

    source 8b41eacdd8811321aa0eb76e72fa164638592c60219ec66ed646ada1bd366d8b
    output a7445045213bc308096149a228ba724916dc881d24a3e3a3fac9c169affcc0fb
    semantic 81ccf950f55485c4bb6484ddb7c30b493fe5500591d12d903cf2fbba34bbfe5a

The source and output are frozen. The [complete independent audit](planar_jc48_sep08_prime_leading_audit.md)
accepts the all-prime, all-m proof, incoming reciprocal-polynomial
supplier, all166 gates and the strict polynomial-mate scope. The named
finite controls support the unbounded argument above.
