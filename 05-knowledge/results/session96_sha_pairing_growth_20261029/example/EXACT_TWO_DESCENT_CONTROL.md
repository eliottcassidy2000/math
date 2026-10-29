# A concrete arithmetic control: the congruent-number curve at one

**CITED METHOD + REPRODUCED CLASSICAL EXAMPLE.** This is a finite local-to-global
certificate for a specified 2-primary Sha, not a new result about this familiar
curve and not a proof of the general finiteness conjecture.

Consider E/Q given by

    y²=x(x-1)(x+1).

The example connects the session's abstract stopping criterion to an actual
elliptic curve. It retains the additive equations of its two-coverings; a
multiplicative square-class skeleton alone does not make the certificate work.

## 1. The exact descent interface

We use the standard full rational 2-torsion descent. A primary source for the
Kummer map and its local interpretation is Morgan–Paterson, [On 2-Selmer
groups of twists after quadratic extension](https://eprints.gla.ac.uk/261997/1/261997.pdf),
§7.2, equation(7.5), with roots0,1,-1. In three-coordinate notation the
Kummer map at a non-torsion point is

    delta(P)=(x,x-1,x+1) in (Q*/Q*²)^3,

whose coordinate product is a square. At a root, the vanishing coordinate
is replaced by the product of its differences from the other two roots.
At infinity the image is(1,1,1). The 2-Selmer group is the subgroup of global
square-class triples that lie in these local Kummer images at every place.
Equivalently it parametrizes the everywhere locally soluble two-coverings

    d1 z1²-d2 z2²=t²,
    d3 z3²-d1 z1²=t²,                                  (1)

in projective3-space, with d1d2d3 a rational square. We import this descent
identification, not a rank computation or a tabulated Sha value.

## 2. There are only sixteen source triples

Every Selmer coordinate has even valuation at every odd prime ell. Here is
a direct check of the local Kummer condition. If v_ell(x)<0, all three of
x,x-1,x+1 have that valuation; their product is a square, so the valuation
is even. If x is integral, at most one of the three is divisible by ell,
and its valuation must be even. The special torsion images obey the same
rule. Thus each global coordinate can be represented by one of1,2,-1,-2.
The third coordinate is determined as a square class by the first two,
leaving sixteen candidates.

This bound keeps all local primes through a proof. It does not prune a
candidate just because a finite prime sample failed to find an obstruction.

## 3. Real signs and the primitive 2-adic guard

Over the real numbers a finite point has x>=1 or -1<=x<=0. The signs of
its Kummer triple therefore have one of the forms

    (+,+,+),  (-,-,+).

The special torsion values and infinity obey these forms too. Eight of the
sixteen candidates remain.

Any Q2-point of(1) can be scaled to a Z2-vector with at least one odd
coordinate. Therefore a necessary condition is a **primitive** solution
modulo8. The companion verifier exhausts all8^4 vectors for each candidate
and checks both congruences. All four surviving-real candidates with
|d1|=2 have no primitive solution. The four remaining triples are

    (1,1,1), (-1,-1,1), (1,2,2), (-1,-2,2).             (2)

Modulo8 survival alone would not prove local solubility. Here a separate
positive certificate does: every triple in(2) is the Kummer image of an
actual rational point, respectively

    O, (0,0), (1,0), (-1,0).

Equivalently primitive rational points on their coverings are respectively

    (z1,z2,z3,t)=(1,1,1,0), (0,1,1,1),
                   (1,0,1,1), (1,1,0,1).

Thus there are at most four Selmer elements by the local exclusions and at
least four by distinct rational Kummer images. The bound is exact.

## 4. Arithmetic conclusion and all-height stopping

The standard Kummer exact sequence gives

    0 -> E(Q)/2E(Q) -> Sel_2(E/Q) -> Sha(E/Q)[2] -> 0.

The four distinct rational images exhaust the four-element Selmer group,
so Sha(E/Q)[2]=0. Any nonzero element of a 2-primary torsion group has a
nonzero multiple of order2. Consequently

    Sha(E/Q)[2^infinity]=0.

Also E(Q) is finitely generated, and
|E(Q)/2E(Q)|=2^rank(E(Q)) |E(Q)[2]|. The latter torsion group has four
elements, so the exact count gives rank zero. No assumption that the entire
Sha is finite was used, and no claim about its odd-primary parts follows
from this calculation alone.

This is an actual finite certificate, stronger than noticing that every
finite pairing model has square cardinality. It succeeds because the exact
local covering equations leave no 2-primary obstruction at the first layer.
A general algorithm producing such a stopping certificate for every curve
and every relevant prime is not supplied.

## 5. Why the additive structure is essential

All sixteen candidates satisfy the multiplicative product-square rule.
Eight even survive the real sign test. The relations involving1 in(1)
remove the remaining four false candidates at2. Forgetting those relations
keeps a strictly larger class set; forgetting primitivity makes the zero
vector pass every congruence. Neither shortcut is a valid F1 descent.
This is a concrete instance of the repository's earlier warning that a
multiplicative root construction still needs its additive branch equations.

The checker verifies the complete finite universes, all four rational
covering witnesses, and the exact image group law. Its finite arithmetic
is combined with the stated descent theorem and valuation proof. It is
not a standalone implementation of Selmer groups for arbitrary curves.
