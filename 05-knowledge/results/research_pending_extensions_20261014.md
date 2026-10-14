# Pending extensions after the completed research checkpoints

14 October 2026. Recovery note for unfinished mathematical follow-ups. This is not a published theorem packet. The completed global surrogate-descent theorem is documented separately in `aliquot_surrogate_global_descent_20261014.md`.

## Genuine prime-power aliquot chambers

Write an even integer as `n=x p^k`, where `x=2^a`, `a>=1`, and `p` is an odd prime. An even successor requires odd `k`. For odd `k`,

`v2(sigma(p^k))=v2(p+1)+v2(k+1)-1`.

The exact size difference is

`(p-1)(s(n)-n)=(2x-p)p^k-(2x-1)`.

Consequently the abundance threshold depends on the base prime: deficiency for `p>2x-1`, equality only for `k=1,p=2x-1`, and abundance otherwise. The surrogate's local core/size tradeoff does not extend: `108=4*3^3` maps to `172=4*43`, increasing while retaining the dyadic core.

If an edge preserves the odd prime base, reduction modulo `p` forces `p | 2x-1`. Such an edge is nondecreasing, with equality only at a perfect number. Thus no nonconstant cycle can use a single fixed odd prime base throughout.

For an odd-length cycle in these chambers, a proposed precise restriction is `gcd_i(a_i+1)=1`. The proof route is short: if this gcd is `d>1`, any prime divisor `ell` of `2^d-1` divides every dyadic divisor-sum factor. Hence `n_(i+1)=-n_i mod ell`; odd closure forces every member to be divisible by `ell`. Their sole odd prime bases must all equal `ell`, contradicting the preceding strict-growth restriction for a nonconstant cycle. This argument is ready for a properly scoped successor note, but no frozen publication packet was completed.

Literature scope: the CWI survey by Garcia, Pedersen and te Riele, section 3.3, attributes to Kanold the stronger classical fact that an amicable pair cannot have exactly two distinct prime factors in each member. Source: https://ir.cwi.nl/pub/4143/04143D.pdf . Original reference reported there: Kanold, *Uber befreundete Zahlen II*, Math. Nachr. 10 (1953), 99–111. Do not present the two-cycle obstruction as a new result. The original publication was not yet inspected in this follow-up.

## Uniform Paley-code obstruction: unfinished extension

The published Golay packet treats the prime-cyclic Paley tournament at 11. A suggested generalization uses the Paley adjacency matrix `A` at order `q=3 mod4` and `M=I+A`:

`M M^T=((q+1)/4)(I+J)`.

At an odd prime `ell | (q+1)/4`, the row space of `M` is self-orthogonal and its rows have zero coefficient sum. Adjoining the all-ones word and extending by minus coefficient sum suggests the same extended-code quadratic identity used at 11. The exact identification with the quotient kernel still needs to be written and checked uniformly.

For `ell>=5`, positional coefficients `0,...,q-1` have sum `1` and square-sum `-1` modulo `ell`. If the quadratic identity applies to a normalized charge equal to a scalar `a`, it forces

`-1+a^2+(1-a)^2=0`, hence `a=0` or `a=1`.

This would exclude the other `ell-2` scalar charges before Hamiltonian constraints. Prime order `q=19,ell=5` is a concrete first target. Literal order `q=27,ell=7` must use the additive group `C3^3`, not a cyclic group of order 27; its appropriate group-algebra quotient remains to be specified. These are proposed extensions, not completed claims.
