# Independent audit of the smaller-endpoint-eight trinomial strip

**Verdict: PASS, relative to the inherited trinomial semigroup classification.**
The proof and all thirty symbolic certificates in
[the moments report](overnight_20260906_moments_width8.md) were independently
reviewed. No universal mathematical step is inferred from bounded sampling
of the return parameter.

The normal-form coverage is complete. A collided first return has
`AB<=a`; for a<=8 this leaves exactly the listed thirty coprime (a,A,B)
families. For a>8,c<=8, the two necessary inequalities
`g|c-b` and `B<g` imply
`a<=g(g-1)-c<=(c-1)(c-2)-c`. Exhausting that finite rectangle using the
original charge equations independently recovers exactly the five listed
collided supports.

The content removal is valid at every admissible integer parameter.
`Ay+Bz=ka` implies `y+z<=ka/A<kg`, since `g>=floor(a/A)+1`. Thus every
nonempty row coefficient is positive. A removed polynomial content factor
cannot vanish there, because it divides every coefficient. Removing
`t^ymin` is legitimate on t!=0, and `v=t^B` maps the complex torus onto
itself. Nonzero primitive resultant therefore excludes every common torus
zero with all three original coefficients nonzero. Extreme-coefficient
normalization and reflection preserve the original moment-zero predicates.

The independent companion reconstructs every raw compressed coefficient
from literal falling products, without importing the producer's code or
calling a resultant routine. It builds integer Sylvester matrices and
evaluates their determinants by fraction-free elimination. For degrees m,n
in v and coefficient degrees d_P,d_Q in g, the determinant has degree at
most `n*d_P+m*d_Q`. At that bound plus one distinct integers, it agrees
with the stored resultant up to a fixed nonzero scalar. The degree bound
therefore proves the full polynomial identity, including unbounded g.
This is a second certificate path rather than a statistical parameter test.

All thirty identities pass, using 674 integer determinant evaluations.
An independent binomial shift reconstructs all 352 strictly positive
coefficients of `R(g+gmin)`. The five opposite-endpoint moment pairs are
also reconstructed directly with their multinomial weights, and their
only common factor is a torus monomial. The audit has 63,843 exact checks.
Normal and optimized outputs agree byte for byte.

- [Independent source](../../04-computation/overnight_20260906_smith_moments_audit.py)
- [Matching output](overnight_20260906_smith_moments_audit.out)

```text
source_sha256:
f3c257aa46bd3b5fe9eb933a09779ae77a9c7f955d203a68f7d7ddf012ed1ee8
output_sha256:
9a32da57667d5da5d88c0158a43ce4cf34f6f21d2c1a321ab04efc09ae285154
```

The proved scope is three-term Laurent polynomials with smaller endpoint
degree at most eight and arbitrary opposite endpoint. It is not an
arbitrary-support Laurent theorem, a full all-trinomial result, or a
coefficient-independent general Gaussian-moment detection bound.
