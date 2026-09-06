# Independent audit of the fixed-anchor Laurent model obstruction

**Status: INDEPENDENT ANALYTIC AND SOURCE AUDIT PASS.** The audited object is
[the anchored model theorem](open_frontier_sep06_laurent.md), with its
[standalone exact producer](../../04-computation/open_frontier_sep06_laurent.py).
This referee accepts uniform model positivity for every real R>=384 and
algebraic same-anchor double cancellation for some R in (256,384). It does
not convert the model into an actual binomial path row.

## 1. Exact parameter and coefficient types

The B roots are `a_i=13R^i/(1+R+...+R^4)`, and the fixed regime is
`q=5,h=4,x=1,g=14,k=9`. The sum of the roots is exactly 13. Applying the
displayed coefficient recurrences gives monic C,D and their next
coefficients -12,-11, respectively. Thus the reversed polynomials have
constant coefficient one and linear anchors `(13,12,11)` without a free
normalizing scale left over.

The continuous shape parameter R is distinct from the Laurent variable.
In Section 3, `t=1/R` is a reciprocal shape parameter; sigma is the
positive magnitude of the negative Laurent phase. All divisions use positive
parameters on the actual domain. The t=0 endpoint is used only for a
polynomial limit, not as a finite-degree admissible model.

For ordered positive roots, `sign B'(a_i)=(-1)^(4-i)=(-1)^i`. Therefore
the ten signed evaluations in the source imply positive C/B' and D/B'
residues. These give strict interlacing with B, including positive simple
C and D roots. At w>0 the signs of `wB^2-2vCD` alternate at 0, each B root,
each intervening C root, and infinity. These ten intervals exhaust its
degree ten and give positive simple roots. At w=0 it is `-2vCD`, with
nonnegative real roots, allowing multiplicities. Reversal gives the stated
all-weight real-rooted pencil, with nonnegative coefficients. No claim about
weights is inferred from a finite weight sample.

## 2. Independent raw-carry and uniform-polynomial check

Write `lambda=13/S_R`. Scaling all base roots by lambda gives

    beta_lambda(z)=lambda beta_base(lambda z),
    W_lambda(z)=lambda^2 W_base(lambda z),
    R_skip,lambda(z)=lambda R_skip,base(lambda z).

The last exponent is one, not two: the crossing factor is the original
`2z`, and both raw contiguous carriers have lower exponent -1. Thus at
`sigma=s_base/lambda`,

    sigma Q_lambda(-sigma)
        = lambda w_base(s_base) + r_base(s_base).

With `s_base=u t^3`, `S_R=Sbar(t)/t^4`, this yields exactly

    t^2 Sbar(t) sigma Q_lambda(-sigma)
      = 13t^6 w_base(u t^3)+t^2 Sbar(t)r_base(u t^3)=T(t,u).

This independently checks the key transport identity with its positive
multiplier. The ordinary coefficient identities in the theorem also retain
the exponent -1 term: the `u^18` hit extraction is divided by s^2, and the
`u^16` skip extraction is multiplied by -2/s.

The direct symbolic construction in the producer was read in full. Its
first polynomial coefficients arise from

    [u^9](1+u)^14 sum_i b_i s^(5-i)u^(2i) / s.

After reciprocal substitution this is precisely the displayed quartic
f(t,u). For the full response, the binomial identity

    [z^j](O^2+z^-1 E^2)=binom(28,2j+2)

justifies the source's direct coefficient formula. Each reciprocal power
has its nonnegative degree checked before reversal. The resulting T and
its substitution `u=1/26+t z` are polynomial identities constructed by
finite convolution and binomial expansion, not fitted identities.

The interval routines perform exact rational Horner evaluation. At each
step they take all four products of interval endpoints, so signed intervals
and intervals containing zero are handled correctly. The bivariate routine
first encloses each coefficient polynomial in t, then applies the same
interval rule in u or z. Correlations are discarded only by enlargement;
this cannot falsely certify a positive lower bound.

The replayed rational certificates establish the following entire-domain
statements, rather than samples:

* all ten normalized interlacer evaluations exceed 1/14 on
  `0<=t<=1/256`;
* the lower and upper f endpoints, after their exact zero constant terms
  are removed and division by t, have opposite strict signs;
* `363<partial_u f<366` on the whole enclosing rectangle;
* `1/40<T(t,1/26+t z)<33/32` for
  `0<=t<=1/384`, `0<=z<=1/100`.

Since f(t,0)<0 and f is strictly increasing through the entire interval
from zero to the upper tube boundary, the trapped root is the unique
smallest positive root. Its derivative is nonzero, so the same branch is
continuous and real analytic. Substitution of this root into the strict T
rectangle gives positive full response for every real R>=384.

The asymptotic constants also agree independently. The first root tends
to u=1/26, hence `sigma/R -> 1/(13*26)=1/338`. The limit
`T(0,1/26)=47439/48334` gives
`Q/R -> 338*47439/48334=47439/143`.

## 3. Algebraic cancellation and exact scope

The negative response at R=256 and the positive response at R=384 occur
on this same admissible first-root branch. Continuity therefore gives
a zero of the full response in the open parameter interval. Clearing
rational R denominators produces two bivariate polynomial equations.

The coprime specialization at R=256 proves their resultant is nonzero.
There is no hidden degree collapse in this inference: the original P has
degree four, whose leading coefficient is a positive binomial coefficient
times the full fifth elementary symmetric function of the positive roots;
the carried polynomial `sigma Q(-sigma)` has degree nine, with a nonzero
leading coefficient from the squared full B term. The skip polynomial has
lower degree. Both degrees remain full for every R>0.

Consequently only finitely many R can give common roots, each algebraic.
The first-root phase is then algebraic over that algebraic R, hence algebraic
over the rationals. The least root-producing R in (256,384) is a valid
exact specification; no numerical zero or guessed resultant factor is used.

The model preserves real-rootedness, the Euler identities, all-weight pencil
compatibility, the original common zero, and all three linear anchors.
It changes the higher coefficients of the genuine composition polynomial
`(1,13,55,84,35,1)`. Thus actual two-rung separation and actual coefficient
transport stay open. Positivity of this model only obstructs an implication
from the retained relaxed predicates, and a proposed equality to a cone of
strictly negative same-zero generators.

The only requested textual repair was a local cross-reference: the R=256
negative control is stated in Section 4 and the source, rather than Section 2.
No mathematical repair was required.

## 4. Reproduction and pins

Both commands were independently rerun:

```
python3 -B 04-computation/open_frontier_sep06_laurent.py
python3 -B -O 04-computation/open_frontier_sep06_laurent.py
```

Normal, optimized, and frozen outputs are byte-identical, with 203 exact
gates. This is a source and proof audit plus independent replay; it is not
advertised as a second independent producer.

* Source SHA256: `37fdc13870f7cf3122f7eebd413cbbfcf8df6a3451744027a8172b1ecdf5f8f0`.
* Output SHA256: `c9d40353e33585d6c1c86d465e09f2862345c8e5af83744595d4357af015e8dd`.
* Semantic digest: `096cdf5d1d82f5d5dd28d2805b365732312d66a6e6ca8b57c288dacdc6cdb247`.
