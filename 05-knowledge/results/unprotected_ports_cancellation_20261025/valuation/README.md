# Exact characteristic spectrum of the regular action cut polynomial

2026-10-25. **PROVED:** the complete coefficient-characteristic spectrum
below, for every odd prime p and integer r>=1. **FINITE-EXACT:** the concrete
50-vertex tournament realizing p=5,r=2. The algebra is independent of the
existence of a tournament with a prescribed full automorphism group.

## The action and its integer polynomial

Put q=p^r and h=p^(r-1). A regular elementary abelian group C_p^r acts
on q gaps. Its identity has q singleton cycles, and each of its q-1 other
elements has h cycles of length p. Its subset-orbit polynomial and positive
shift are therefore

    P(t)=((1+t)^q+(q-1)(1+t^p)^h)/q,
    Q(u)=-P(-1-u)=[u^q+(q-1)((1+u)^p-1)^h]/q.       (V1)

Both are integer polynomials. In particular P counts actual subset orbits;
integrality of Q also follows by its integer change of variable. Q has
nonnegative coefficients. See [the positive-jet proof](../counting/POSITIVE_JET.md)
for the general tournament construction, the distinction between full
automorphism groups and chosen subgroups, and an explicit C_5^2 certificate.
Here the depth is the multiplicity of u=0 in Q over the specified field.

Define the integral polynomial

    C_p(u)=((1+u)^p-1-u^p)/p
          =sum_(i=1)^(p-1) binom(p,i)/p * u^i.       (V2)

Its linear coefficient is one. Expand (V1), combining the j=0 term with
the separate u^q term before dividing:

    Q(u)=u^q+(q-1)*sum_(j=1)^h
                  [binom(h,j)*p^j/q]*u^(q-p*j)*C_p(u)^j. (V3)

Each scalar in square brackets is an integer. Indeed

    v_p binom(h,j)=r-1-v_p(j),
    v_p [binom(h,j)*p^j/q]=j-1-v_p(j).               (V4)

For the first identity, use binom(h,j)=(h/j)binom(h-1,j-1).
Modulo p, the coefficients of (1+x)^(h-1) are (-1)^i, by dividing
(1+x)^h=1+x^h by 1+x. Thus the second binomial factor is a p-unit.
The only possible denominator of the bracket in (V3) is a power of p,
and its valuation is nonnegative. For odd p and j>=2 it is at least one.
The j=1 bracket is exactly one.

## All coefficient characteristics

Over characteristic zero the exact depth is h, and its leading coefficient is

    [u^h] Q=(q-1)*p^(h-r).                           (V5)

This follows directly from (V1): the powered factor begins with p^h*u^h,
and h<q. The exponent h-r is nonnegative for odd p and r>=1.

In the physical characteristic p, reduce the integral identity (V3):

    Q(u) mod p = u^(q-p)*(u^p-C_p(u)).               (V6)

The coefficient of u^(q-p+1) is -1. Hence the exact depth is q-p+1.
The residual polynomial has the finite-logarithm form

    C_p(u) mod p = sum_(i=1)^(p-1) (-1)^(i-1)*u^i/i. (V7)

This is a polynomial identity in F_p[u]; every displayed i is a unit.
It does not require a choice of logarithm or convergence. Crucially, (V6)
is obtained after integer division, not by dividing by q in F_p.

There is a separate exceptional regime. If a prime ell divides q-1, then
q is a unit modulo ell and (V1) gives Q(u)=u^q modulo ell. If ell divides
neither p nor q-1, the coefficient (V5) survives. Thus the complete answer is

| Coefficient field characteristic | Exact depth |
| --- | --- |
| 0, or prime ell not dividing p(q-1) | h |
| p | q-p+1 |
| prime ell dividing q-1 | q |

The regimes are disjoint and cover every field. In particular characteristic
two always lies in the last row, since q is odd. A depth increase is not
exclusive to the characteristic matching the underlying group prime.
At r=1 the first two depths coincide at one. The modular jump in the middle
row is strict for r>=2. The proof requires odd p.

## A genuine tournament and the scope of the bridge

The explicit 25-vertex Cayley tournament in
[POSITIVE_JET.md](../counting/POSITIVE_JET.md#4-a-concrete-50-vertex-certificate)
has full automorphism group C_5^2. Replacing each vertex by a transitive
pair gives a 50-vertex tournament whose canonical gap action is regular
on 25 gaps. Thus its invariant constant-XOR-toggle rank is one, while

    depth over Q = 5,
    depth over F_5 = 21,
    depth over F_2 or F_3 = 25.

Its rational leading coefficient is 3000, and

    Q mod 5 = 4u^21+3u^22+3u^23+4u^24+u^25.

The checker and [independent review](INDEPENDENT_REVIEW.md) retain full
automorphism and canonical-clone control. For an unbounded genuine family
at p=5, the positive-jet note uses the cited tournament regular representation
theorem. The algebra above is valid for p=3,r=2,3, but that cited existence
theorem specifically excludes those two elementary abelian groups. An
action formula alone must not be promoted to a full-automorphism claim.

The connection to the [cyclic wraparound calculation](../bridges/README.md)
is precise but limited. Both require integral data before specialization
in the physical characteristic. The cyclic calculation loses division by p
when extracting its constant mode; the cut polynomial loses division by
p^r in its orbit average. Here (V3) supplies the missing integral carrier.
The resulting residual finite logarithm is a new explicit object for the
connection. There is no identification of the cyclic boundary operator
with a black-line move, nor a consequence for black-core connectivity.

There is also a literal arithmetic carry behind (V2). Define

    C_p(x,y)=((x+y)^p-x^p-y^p)/p.

This is integral and satisfies

    C_p(x,y)+C_p(x+y,z)=C_p(y,z)+C_p(x,y+z),

by cancellation of the four p-th powers. The polynomial in (V2) is
C_p(1,u). Thus the residual finite logarithm is the specialization of
an associative additive carry for p-th powers. This elementary identity
gives an explicit map from the earlier arithmetic-carry theme to the
present orbit polynomial; it supplies no claim about a Collatz orbit.

## Controls and next question

Run `python check.py` and `python -O check.py`. The checker compares literal
binomial substitution of P with an independent polynomial power in (V1),
checks integrality and positivity, all valuations in (V4), (V6)-(V7), and
the complete characteristic spectrum for the listed finite universe.
It includes the r=1 equality boundary and formal p=3 cases, and explicitly
tests the additional characteristic-two jump as a hostile to an exclusive
physical-characteristic explanation. No Python assertions are used as gates.

**OPEN question:** for a general odd-order gap action, can an integral
cycle-index refinement predict the full characteristic spectrum without
expanding every subset orbit? The rational minimum-cycle formula alone
cannot do this: it omits the divisibility of its positive leading coefficient.
The present elementary abelian family supplies an exact test case for any
proposed answer, with both a bad-denominator regime and a vanishing-multiplicity
regime that must remain distinct.

A cheap orthogonal probe replaces the elementary abelian regular action
by the cyclic regular action C_(p^r). Its exact positive polynomial is

    Q_cyc(u)=1/p^r * sum_(a=0)^r phi(p^a)
                            *((1+u)^(p^a)-1)^(p^(r-a)).

**FINITE-EXACT:** for (p,r)=(3,1..3),(5,1..3),(7,1..2), its rational
depth is one and its characteristic-p depth is p^(r-1). Thus at q=25
the two regular actions already have different modular depths,5 versus21.
The checker computes these eight cases directly. **OPEN in this packet:**
the suggested all-r cyclic modular-depth formula. This is a precise next
valuation problem, not an asserted extension of (V6), whose proof uses
the elementary abelian cycle inventory.

One all-order comparison is already exact. Differentiating the integral
cyclic polynomial cancels the q denominator and gives

    Q_cyc'(u)=sum_(a=0)^r phi(p^a)*(1+u)^(p^a-1)
                            *((1+u)^(p^a)-1)^(p^(r-a)-1).

Modulo p only a=0 and a=1 survive. Consequently

    Q_cyc'(u) = u^(q-1)-(1+u)^(p-1)*u^(q-p) mod p
              = Q'(u) mod p.                       (V8)

Both integer polynomials are monic of degree q and vanish at zero.
Since the ordinary derivative kernel in F_p[u] is F_p[u^p],

    Q_cyc-Q mod p lies in span{u^p,u^(2p),...,u^(q-p)}. (V9)

This is an (h-1)-dimensional space, empty at r=1. Thus the ordinary
derivative can be identical while the root depths differ drastically.
The q=25 case has depths5 and21 despite (V8). Coefficients at positive
p-multiple degrees are the missing coordinates; Hasse coefficients can
retain them. This is a second explicit observation-kernel bridge to the
older cyclic work, beyond the common theme of division by the physical prime.
