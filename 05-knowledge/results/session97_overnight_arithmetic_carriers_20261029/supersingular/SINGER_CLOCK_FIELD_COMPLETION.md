# Additive ports plus a named Singer clock and identity recover a field

**Session97, 2026-10-29. Status:** PROVED elementary completion theorem;
FINITE-EXACT verification for all six reconstructed F32 carriers. The finite
field/Singer construction is classical. This note contributes an explicit
source-preserving synthesis with the repo's joint-port reconstruction, not
a new classification of finite fields or an F1 realization.

Only this note and its checker are new. Earlier checkpoints are unchanged.

## 1. Exact hypotheses and theorem

Let `(V,+,0)` be an F2-vector space of dimension m≥1. Let sigma be a **named
linear Singer clock**: an invertible linear map acting as one cycle on the
`2^m-1` nonzero vectors. Choose a named `e≠0` to serve as multiplicative
identity. Every nonzero vector has a unique expression `sigma^i(e)`, with
i taken modulo `2^m-1`.

Define zero to be absorbing and define

`(sigma^i e) * (sigma^j e) = sigma^(i+j) e`.                (1)

Then `(V,+,*)` is a field of order `2^m`. It is the unique field structure
with the given addition and identity e **in which sigma is multiplication
by sigma(e)**. Merely requiring sigma to remain a linear map would not be
a uniqueness condition; the realization as a multiplication operator is
part of the statement.

**Proof.** On the nonzero vectors, (1) is the cyclic group law on exponents:
it is associative and commutative, has identity e, and every element has
an inverse. For `x=sigma^i(e)`, left multiplication by x is exactly the
linear map `sigma^i` on all V, including zero. Hence

`x*(y+z)=sigma^i(y+z)=sigma^i(y)+sigma^i(z)=x*y+x*z`.

This proves distributivity; the case x=0 is immediate. The additive vector
space supplies an abelian group, so all field axioms hold. In any other
field satisfying the stated realization, multiplication by sigma^i(e)
must again be sigma^i, forcing (1). No nontrivial counting argument is used.

In the earlier all-m theorem, four-port observations recover the addition
for m≥4 after marking zero. The present theorem completes that additive
carrier once the linear clock and its identity phase are also supplied.
It does not infer those extra marks from the four-port counts alone.

## 2. What is retained, what is redundant, and what cannot be forgotten

* **Additive ports** determine the affine geometry and, after marking zero,
  its addition. They do not by themselves specify multiplication.
* **The Singer orbit** supplies a cyclic multiplicative skeleton, but without
  a chosen e there is no distinguished multiplicative identity. For m>1,
  different e give genuinely different operations on the fixed named set.
  The m1 field F2 is the trivial exception.
* **Linearity of the clock** is load-bearing: it is exactly what proves
  distributivity. A permutation of the same order and cycle type is not enough.
* **The named generator of the clock subgroup** contains more information
  than the field multiplication alone. Replacing sigma by sigma^k with
  `gcd(k,2^m-1)=1` and keeping e yields the same multiplication, but a different
  named one-step clock. Thus the cyclic subgroup plus e suffices for the
  operation, while later queries about the source's next clock step require
  the generator to remain marked.
* Exponents in (1) are residues modulo the period. They are **not an unbounded
  elapsed-time counter**. A later operation sensitive to integer winding
  needs that separate clock coordinate.

If `e'=sigma^t e`, the new multiplication is

`x *' y = sigma^(-t)(x*y)`.

The linear map sigma^t is a clock-preserving field isomorphism from the old
structure to the new one and transports e to e'. These are isomorphic fields,
not identical operations on the named carrier. This is precisely the source
versus unmarked-isomorphism distinction the construction must preserve.

A simple missing-clock hostile uses the same F2^3 addition (bitwise XOR)
and the same identity1. Polynomial multiplication modulo `x³+x+1` gives
`2*4=3`; modulo `x³+x²+1` it gives `2*4=5`. Both irreducible cubics define F8,
but they define different products in the fixed bit chart. The corresponding
multiplication-by2 Singer clocks differ. Addition and identity alone do not
select between them.

## 3. All six length32 carriers, with the phase conventions explicit

The six-code affine reconstruction used labels0,...,31, with additive zero
at label31. The source clock is

`0→1→...→30→0`, with label31 fixed.

Here **label0**, not label1, is chosen as the multiplicative identity e.
Let `alpha=sigma(e)`; it is label1. The power basis
`e,alpha,alpha²,alpha³,alpha⁴` consists of labels0,1,2,3,4.
The four-port addition then determines the primitive polynomial of alpha:

| Source code generator | Primitive polynomial bitmask |
|---:|---:|
| 73235 | 41 |
| 83903 | 55 |
| 92135 | 47 |
| 102641 | 37 |
| 118733 | 61 |
| 129989 | 59 |

Bit i is the coefficient of x^i. All six additions are distinct on the
named set, although the orbit clock and formula (1) give the same
multiplication table. Each combination is a field. Thus the source's
four-port geometry has not been discarded after naming the clock: it
selects which additive structure is compatible with the marked primitive
power basis. The checker reconstructs these additions again from octad
incidence rather than assuming a field chart, and then verifies every
addition/multiplication pair against exact polynomial arithmetic.

## 4. The two Frobenius meanings must not be conflated

In the reconstructed F32, field Frobenius is `F(x)=x²`, which on nonzero
clock coordinates sends exponent i to2i. It is **not sigma**, which adds1
to the exponent. The exact compatibility is

`F sigma = sigma² F`.

The checker verifies this identity and additivity of F for all six carriers.

In the earlier arithmetic examples, the source31-cycle could arise as
**F5-Frobenius on the finite branch roots** of a curve. The newly reconstructed
F32 is an auxiliary algebra on the port set; its field characteristic is2.
The source branch permutation becomes multiplication by a primitive element
in this algebra, not its squaring Frobenius. No embedding of F32 into the
characteristic-five branch splitting field, and no Jacobian endomorphism
identification, is claimed. The source action is transported as a named
permutation, with its interpretation explicitly changed.

## 5. Hostile clock certificate

For source code73235, retain its reconstructed addition and the identity
label0. Conjugate the source clock by the permutation swapping labels1 and2
and fixing every other label. The resulting clock is still a31-cycle on
nonzero labels and fixes additive zero31, but is not linear.

Use its orbit to define multiplication by (1). For `(a,b,c)=(1,0,1)` the
checker obtains

`a*(b+c)=label16`, but `a*b+a*c=label6`.

Thus the cyclic multiplicative skeleton remains a group, while the additive
branch equation fails. This is an exact reason to retain both structures
and the linear compatibility guard; neither cycle size nor recurrence
period supplies the missing condition.

## Reproduction and stopping boundary

Run `python check_singer_field_completion.py` and also with `python -O`.
Outputs agree. The checker reconstructs all six additions from fourth-order
incidence, verifies the clock's linearity and the field axioms, compares
all products with the six primitive-polynomial charts, tests both Frobenius
interpretations, changes the identity phase, and supplies the nonlinear-clock
and missing-clock hostiles. It uses only Python's standard library.

This closes the specified port-plus-clock completion task. It does not
extend the family search, construct an F1 scheme, or infer new statements
about Sha, integer Collatz or unrestricted tournament composition.
