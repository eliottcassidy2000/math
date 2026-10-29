# Fermat certificates: a constructive projection and a bound on fixed recipes

Write `F_n=2^(2^n)+1`. PROVED here: F1 and F2, with elementary proofs below.
F1 organizes classical finite abelian group facts; F2 is a scoped consequence
of cyclotomic irreducibility and roots of unity. No novelty priority is claimed.
FINITE-EXACT: witnesses at `5<=n<=10` and specified hostile cases. OPEN:
compositeness for **every** `n>=5`. No heuristic about the number of further
Fermat primes is used in either proof.

## F1. A moving certificate that is not blind like source 2

Set `m=2^n`, `N=F_n`, and `G=(Z/NZ)^*`. Define the computable homomorphism

`h_n:G -> G`, `h_n(a)=a^(N-1) mod N = a^(2^m) mod N`.

**Theorem.** Its image is exactly the odd-order subgroup `G_odd`. Its fibers
all have size `|G_2|`, the order of the 2-primary subgroup. Moreover,

`N composite <=> G_odd != {1} <=> exists a in G: h_n(a)!=1`.

Conditional on compositeness, at least `2/3` of all units supply a witness;
the exact fraction is `1-1/|G_odd|`.

*Proof.* If `N=product q^e`, CRT decomposes `G` into cyclic groups of orders
`q^(e-1)*(q-1)`. Their 2-parts have exponent `2^(v_2(q-1))<=2^m`, since
`q<=N`. Raising to `2^m` kills these parts and acts bijectively on every
odd-order part. This proves the image and equal-fiber statements. It is
surjective onto the odd subgroup, but is not necessarily the identity there,
so "projection" here does not assert idempotence.

If `G_odd` is trivial, each exponent `e` must be 1 (otherwise the odd prime
`q` divides the unit-group order), and each `q-1` must be a power of 2. A prime
of the form `2^r+1` has `r` a power of 2: an odd factor of `r` would factor the
binomial. Thus each `q` is a Fermat prime `F_j`. Pairwise coprimality of
Fermat numbers forces `j=n`; hence `N` itself is prime. The converse follows
from `|G|=N-1=2^m`. If the odd subgroup is nontrivial its order is at least 3,
which proves the fraction bound. QED.

This proof assumes neither squarefreeness of Fermat numbers nor their
factorizations. Factorization is used only to prove the group description.
Given any proposed `a`, a nonunit gives a factor via gcd, and a unit with
`h_n(a)!=1` gives a directly checkable Fermat compositeness witness. The
algorithm need not know the group order or factors to verify that witness.

The unit witness `x=h_n(a)` has odd order; conversely any nontrivial odd-order
unit is a certificate, because a Fermat prime's unit group has power-of-two
order. This is the earlier
[odd-order/squaring-digraph bridge](../four_vertex_glue_fermat_20261029/FERMAT.md),
related to [Somer–Křížek (2004)](https://dml.cz/handle/10338.dmlcz/127904).
The constructive step here is to retain the source `a` and the full clock of
`m` squarings, rather than inspect only the universally available source 2.

For `n>=1`, `h_n(2)=1` whether `F_n` is prime or composite: the exact order of
2 is `2^(n+1)`, which divides `2^(2^n)`. Base 2 is therefore a mandatory
hostile control. At `n=0`, `F_0=3` and `h_0(2)=1` as well. In contrast, the
checker verifies `h_n(3)!=1` for `5<=n<=10`, using both binary modular powering
and an independent loop of `2^n` squarings. It gives `h_n(3)=1` for `1<=n<=4`.
**This finite experiment does not assert that base 3 works for every composite
Fermat number.** For a complete primality decision, the classical test is
Pépin's test with the half exponent; it is a different claim.

At `F_5=641*6700417`,

`|G_2|=2^14=16384`, `|G_odd|=5*52347=261735`.

The unit-witness fraction is `261734/261735`. The inherited explicit unit
`3303966555` has order 3 modulo `F_5`; its square is `3041328343`. This gives
a positive certificate even though base 2 remains silent. Checks at
`9,25,45,65,85,91,257` explicitly include repeated prime factors in the
general Hall-fiber formula.

The quantifier gap remains exact: the fraction statement starts **after**
assuming `F_n` composite. It cannot be reversed into a proof that the witness
set is nonempty for every `n>=5`. An exhaustive search terminates exactly
when the corresponding compositeness assertion is true.

## F2. Fixed algebraic charts cannot supply uniform odd-order certificates

The incoming THM-4847
[`fermat-primes-barriers-to-proving-exactly-five-the-function-field-analogue-and-why-tournament-counts-are-blind-while-three-state-counts-are-pepin`](../../../01-canon/theorems/THM-4847-fermat-primes-barriers-to-proving-exactly-five-the-function-field-analogue-and-why-tournament-counts-are-blind-while-three-state-counts-are-pepin.md)
recovers the finite-prime-cover obstruction and the difference between
cyclotomic factorization over finite fields and over the rationals. The next
statement gives a more quantitative restriction on one proposed family of
certificates. It does not rule out all formulas, all algorithms, or all
uniform certificates.

Fix `s>=0`, set `D=2^s` and, for `n>=s`, set

`t_n=2^(2^(n-s))`, so `F_n=t_n^D+1`.

Choose **fixed** `P,Q in Z[X]` and a fixed odd integer `d>1`. At indices where
`gcd(Q(t_n),F_n)=1`, consider the rational recipe

`x_n=P(t_n)*Q(t_n)^(-1) mod F_n`.

**Rigidity theorem.** If `x_n^d=1 mod F_n` for infinitely many admissible
indices, then `x_n=1` at every admissible index. Consequently a finite bank of
such fixed charts, with a finite set of odd certificate exponents, cannot
supply a nontrivial odd-order certificate for every `n>=5`.

*Proof.* Divide `P^d-Q^d` by the monic polynomial `f=X^D+1`, giving an integer
remainder `R` of degree less than `D`. The congruence is exactly
`f(t_n)|R(t_n)`. If it holds infinitely often, growth of `t_n` and
`deg R<D` imply that eventually `R(t_n)=0`; hence `R=0` as a polynomial.

Now `f=Phi_(2D)` is irreducible over Q. At a primitive `2D`-th root `zeta`,
`Q(zeta)!=0`: otherwise `f|Q` and no index would have an invertible denominator.
Thus `P(zeta)/Q(zeta)` is a root of unity of odd order in `Q(zeta)`.
That field's roots of unity all have power-of-two order, so the ratio is 1.
Irreducibility yields `f|(P-Q)`, and therefore `x_n=1` at every admissible
index. The finite-bank conclusion follows by pigeonhole. QED.

The root-of-unity input can also be checked by degrees: if this field contained
a primitive odd `r>1` root, it would contain a primitive `2D*r` root, forcing
degree `phi(2D)*phi(r)>phi(2D)`. The elementary cyclotomic facts are described
in [Ben Lynn's cyclotomic-field notes](https://crypto.stanford.edu/pbc/notes/numberfield/cyclo.html).

**Effective certificate-size bound.** Put `S=||R||_1`, the sum of absolute
coefficients. If a recipe gives a nontrivial odd-order certificate at `t=t_n`,
then `R!=0` and

`t <= 1+S`.

Indeed, when `t>1+S`, the leading coefficient dominates the lower coefficients,
so `R(t)!=0`, while `|R(t)|<=S*t^(D-1)<t^D+1`. Divisibility is impossible.
The `R=0` case is already trivial by the preceding proof. This bound is
deliberately conservative; optimizing its additive constant is irrelevant.

Reduce `P,Q` modulo `f` first. If their coefficients have absolute value at
most `B>=1`, then `||P||_1,||Q||_1<=DB`, and reduction modulo `X^D+1` cannot
increase coefficient 1-norm. Therefore any nontrivial certificate obeys

`2*(D*B)^d >= t_n-1 = 2^(2^(n-s))-1`.

This inequality also applies separately to *varying* recipes. With fixed depth
`s`, bounded coefficient height and bounded odd exponent cannot suffice at
arbitrarily large `n`. At least one must grow. Increasing chart depth is a
third escape. The moving exponent in F1 is outside the fixed-recipe theorem.
This gives a concrete design constraint rather than a claim that creative
uniform constructions are impossible.

## Hostile tests and the next positive construction

Oddness is essential: the fixed recipe `x_n=t_n` has exact even order `2D`
modulo `t_n^D+1` for every `n>=s`. Thus replacing "odd" by "arbitrary" in F2
is false. Invertibility of `Q(t_n)` and constancy of the recipe are essential
typing conditions. A factor-dependent CRT construction is allowed, but then
the factor and relative phase are extra arithmetic inputs, not a free chart.

The checker evaluates 990 admissible samples from lags `s=0,1,2`, certificate
exponents `3,5,7`, six numerators and three denominators, through `n=8`.
It includes the constant `3303966555` as a positive order-3 example at `F_5`
in each of the three charts (three successes), as well as the uniform
even-order hostile family. These finite tests validate the implementation
and the stated boundaries; the proof of F2 is the symbolic argument above.

The useful next search is for a family of **moving** sources and odd clocks
with a demonstrable nonidentity condition, not merely a new picture of the
squaring graph. Record `(n,a,m,h_n(a))`, or `(n,x,d)` with `d` odd,
`x^d=1`, and `x!=1`, as the certificate. Any proposed bounded-degree chart
should first be reduced modulo `X^(2^s)+1` and subjected to the size bound.
For Collatz, retain forcing and numerator under attachment; for Fermat,
retain source and the growing certificate clock. These are analogous
information-retention requirements, not a proved map between the conjectures.
