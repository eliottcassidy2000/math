# Polynomial evaluation gives an exact decreasing horizon quotient

2026-10-25. **PROVED:** a polynomial form of the moment transfer. It is an
equivalent representation of [the proved recurrence](README.md), not an
additional storage or running-time claim.

For a bank with fixed defect profile D and remaining total degree S>0,
encode its signed cycle-count weights by

    H(Y)=sum_m w_m Y^m,
    A_S(Y)=product_(r=1)^S (Y-2^r).

The retained moments are exactly H(2),H(4),...,H(2^S). Over Q these nodes
are distinct. Therefore banks have the same future responses exactly when
their polynomials are congruent modulo A_S. Polynomial remainder, evaluation
vector, and the associated quotient algebra Q[Y]/(A_S) are equivalent
S-dimensional presentations. For a bank restricted to h distinct admissible
exponents m the rank is min(S,h), as proved by the transpose Vandermonde
argument in the main note.

For a nonzero batch (ell,k,e), put S'=S-ell(k+e) and

    c=beta(k,e)/(ell^k*k!) * 2^[A_ell(k)-k*D(ell)].

Its action before quotienting is

    H'(Y)=c*Y^k*H(2^(ell*k)*Y).                      (E1)

If S'>0, then for every destination node Y=2^r,1<=r<=S', the argument
of H is the source node 2^(ell*k+r), with 1<=ell*k+r<=S. Consequently

    A_(S')(Y) divides A_S(2^(ell*k)*Y).               (E2)

Thus (E1) descends to a well-defined linear map of the two evaluation
quotients. The multiplier Y^k and coefficient c do not change this fact.
When S'=0<S, every nonzero completing batch has k>=1; its terminal sum
is H'(1)=c*H(2^(ell*k)), again an available source evaluation. An already
terminal state carries its scalar unchanged. This proves exactness at
the boundary as well as in the interior.

The hostile bank is itself the small polynomial

    H(Y)=8Y-6Y^2+Y^3=Y(Y-2)(Y-4).

It is zero modulo A_2 but H(8)=192. The quotient is therefore safe under
the decreasing horizon (E2), and unsafe when reused at a larger horizon.
This identifies the missing coordinate exactly; it is not evidence that
the underlying signed histories are absent.

## Relation to cyclic Hasse coordinates

The [cyclic difference model](../bridges/README.md) in characteristic p
uses k[epsilon]/(epsilon^p), a local algebra with a nilpotent coordinate.
Here Q[Y]/(A_S) is a product of S copies of Q, since its evaluation nodes
are distinct. Both describe the kernel of a specified observation map,
but their algebras and allowed continuations differ. One uses derivatives
at a repeated root; the other uses values at distinct roots. Their common
requirement is to retain exactly the coordinates that future operations
can read.

In an odd coefficient characteristic ell, the nodes 2^r repeat with
period ord_ell(2). This suggests a smaller modular observation algebra.
It does not yet give an integer-counting algorithm: the transition contains
denominators ell_length^k*k!, which need not be units in that field. The
characteristic-jump calculation in [the valuation note](../valuation/README.md)
is a concrete hostile to reducing an orbit average before its integer
division. A valid modular recurrence must first supply an integral
normalization and enough precision for these divisions. That is the
next question; the rational quotient above is already exact.
