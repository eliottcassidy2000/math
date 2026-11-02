# Tensor compensation for every integer residue modulus

**PROVED, elementary extension of the inherited tensor lemma; same-agent
audit.** This strengthens the arithmetic statement while leaving its
tournament application at the actual parity hypothesis.

Let m>=2, a nonzero modulo m, and delta=min_(t in Z)|a+mt|. Let P be an
integer multiaffine polynomial whose top coefficient is congruent to a
modulo m and all proper coefficients are divisible by m. For |x_i|<=1,

    sum_(S subset[r]) (partial_S P(x))^2
        >=delta^2 product_i(1+x_i^2).                    (C)

This is sharp: choose the top coefficient to be a least absolute residue
of a and all others zero. Divisibility and the full derivative jet are
essential. The one-variable polynomial x-1 at x=1 has jet norm1<2 if the
proper-coefficient restriction is dropped. A single polynomial value can
vanish even with the restriction, as3x-2 at x=2/3 modulo2 demonstrates.

## Inductive proof

Use the coefficient-reversed jet convention of the inherited
[tensor proof](../session139_incidence_moduli_and_stability_20261102/tensor_matching/PROOF.md).
The empty coefficient lies in a+mZ and every other coefficient lies in mZ.
Split at the final coordinate into vectors (u,v+xu). The coefficient
vector producing v is entirely in mZ; consequently u,u+v,u-v all belong
to the same permitted coset in one lower dimension. By induction their
squared norms are at least L=delta^2 product_(i<r)(1+x_i^2).
Put A=||u||^2,C=||v||^2,D=<u,v>. Then A>=L and A+C+-2D>=L, whence

    A+||v+xu||^2
      >=(1-|x|+x^2)A+(1-|x|)C+|x|L
      >=(1+x^2)L.

The base is the distance delta from the integer coset to zero. No field,
primality, or assumption that a is a unit modulo m is required.

## Eisenstein proof without an inert-prime restriction

The inherited [Bernstein proof](../session139_incidence_moduli_and_stability_20261102/stability/EISENSTEIN_BERNSTEIN_TENSOR_PROOF.md)
uses Q0=I,Q1=[[1,1/2],[1/2,1]],Q2=[[1,1],[1,2]] and a unimodular
integral factor of Q2. Put zeta=exp(i*pi/3). Any evaluation at zeta or
its conjugate is congruent modulo m Z[zeta] to a times a power of zeta.
That power is a unit even when Z[zeta]/m has zero divisors (for example
m=3 or6). Multiplying by its inverse preserves the complex norm and puts
the evaluation in the coset a+m Z[zeta].

Every element of this coset is z=(a+mb)+mc zeta. If c=0 its norm squared
is at least delta^2. Otherwise

    |z|^2=(a+mb+mc/2)^2+3m^2c^2/4
           >=3m^2/4>delta^2,

because delta<=m/2. Thus the exact norm floor of this Eisenstein coset is
delta^2. This replaces the older F4 nonvanishing argument with a stronger
norm statement. The remaining Bernstein argument is unchanged: each
type2 coordinate supplies two retained slices in the same top-residue
coset, and the positive Bernstein weights yield(C). Negative coordinates
are handled by sign changes; delta is unchanged on replacing a by-a.

## What actually transfers to tournaments

The perfect-matching parity calculation for principal Pfaffians gives the
required coefficient divisibility modulo2. It does not generally give it
modulo3,5,11 or a composite modulus. For example the proper coefficient-2
of3x-2 is zero modulo2 but not modulo5. Therefore(C) supplies new norm
certificates whenever a stronger residue coset is independently proved;
it does not automatically strengthen every tournament Pfaffian bound.
This preserves the exact boundary between a general tensor theorem and
the arithmetic data available in its application.

The arithmetic commonality with prime wheels is integral residue
information. The invariants differ: coset distance controls jet energy,
whereas collision partitions control boundary deletion. A shared prime
or the same cardinality alone supplies no map between those predicates.
