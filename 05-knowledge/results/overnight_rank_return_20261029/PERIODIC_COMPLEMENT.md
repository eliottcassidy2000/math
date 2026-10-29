# An actual complement recurrence: 35 boundary states, order 18, and an F243 clock

2026-10-29, 12:00 UTC continuation. **PROVED scoped all-size recurrence and
observer quotient, with FINITE-EXACT transfer/Smith certificates.** Same-agent
review. This is a particular periodic family of literal complements, not
the general strong-atom or re-marking carrier.

## E1. The source and its finite transfer

Use [DEFECT_TRANSFER D2--D4](DEFECT_TRANSFER.md). For n=2m, complement the
interval tournament S_n along the path with m runs of length two. Let b_r,
where n=8+2r, be its native three-run count with starts (0,1,2) and terminals
(n-3,n-2,n-1). The defect pairs are exactly

    (2i+1,2i+2), (2i,2i+3),       0<=i<m-1.             (E1)

Introduce the first four vertices with the three named starts, repeat a
two-vertex bulk rule r times, and introduce the last four vertices with
the three named terminals. At each bulk boundary the only active vertices
are the last even and odd vertices. Rename them 0,1. The next vertices 2,3
have backward neighbours {1},{0}; after their introduction rename 2,3 to
0,1. Component pairings and anonymous tails are retained exactly as in D2.

Closing the reachable state set under this bulk rule gives 35 states. The
receipt contains every state, every integer transition, the initial row a
and terminal column b. It proves

    b_r = a B^r b,          B in Z^(35 by 35).           (E2)

This is an all-size finite-transfer construction, not a fitted numerical
sequence. The explicit prefix, suffix, all outgoing transitions and closure
are rebuilt from component joins. Twenty initial sizes agree with the
independent general frontier evaluator; the values at n=20,100,1000 also
agree. The ancestor is the March run family, through its October literal
cut/hinge rule; no new priority is claimed for that decomposition.

## E2. Minimal scalar order eighteen over Q

Let

    f6(X) = X^6-X^5-6X^4+4X^3-7X^2+5X-4,
    f11(X)= X^11-13X^10-16X^9-268X^8-382X^7-1170X^6
            -1048X^5+2072X^4+3197X^3+459X^2-1232X-576,
    p(X)=(X+1) f6(X) f11(X).                            (E3)

**PROVED with an exact finite certificate:** p(B)b=0. Thus p is an
all-size homogeneous recurrence for (b_r), starting at r=0. The receipt
checks the entire 35-coordinate zero vector, not only an initial stretch
of scalar agreement. The first eighteen values are recorded; they begin

    6,38,684,10712,163650,2502194,38336592,587179288,...

The leading 18-by-18 Hankel determinant is nonzero:

    47433022093166182715564496945519191748739613768199188883059268875523574484092755413006479858071804837888.

Therefore the scalar recurrence order is exactly eighteen over Q. This
also bounds below every rational linear realization of this scalar series.
It is not a lower bound on every nonlinear or realization-constrained
combinatorial representation. The displayed factor degrees 1,6,11 are
from this actual transfer, not a universal rule for all complements.

## E3. A smallest closed quotient for a declared boundary language

Here is a precise, limited carrier problem. Supply arbitrary boundary
weights w in R^35. Allow only repeated bulk insertion B followed by the
terminal observation b. Two weights are equivalent when every future
observation agrees. Retain

    Psi(w) = O w^T = (w.b, w.Bb, ..., w.B^17b)^T,        (E4)

where O is the 18-by-35 integer observation matrix (the formula uses w as
a row). The monic identity p(B)b=0 supplies an integer companion update
on these eighteen observations, over every commutative ring R.

**PROVED.** Equality in (E4) is equivalent to equality of every future
bulk observation. Hence image(O) is the smallest carrier for this declared
weighted observer language: any exact representation must distinguish its
different elements, and the companion rule closes its evolution.

The exact Smith factors of O are

    1,1,1,1,2,4,4,4,4,4,8,8,72,288,576,576,2304,L,
    L=6616316393038183526006648240572416.                (E5)

Independent prime-power elimination checks every nonunit prime layer.
At modulus m its exact state count is

    |image(O mod m)| = m^18 / product_i gcd(m,d_i).       (E6)

In particular the minimum sizes are 16 at modulus two and 3^12 at modulus
three. The source is arbitrary weights on these 35 specified boundary
states, not all tournament attachments. Actual prefixes lie in a smaller
realization set; the particular initial trajectory below has only one
binary state and 726 ternary states. Nor does this operation language
include another complement or arbitrary re-marking.

This separates three previously easy-to-confuse quantities: 35 combinatorial
boundary types, 18 rational observable dimensions, and precision-dependent
finite carrier cardinalities. The full target relation remains supplied.

## E4. The actual ternary recurrence has an F243 chart

Put, over F3,

    g(X)=X^5+X^4-X^3-X^2-1.                             (E7)

The integral polynomial satisfies

    p mod 2 = X^3 (X+1)^15,
    p mod 3 = X (X+1)^2 g(X)^3.                         (E8)

But the actual scalar series modulo three has the smaller MINIMAL polynomial

    q(X)=(X+1)^2 g(X)^2,              degree 12.         (E9)

The verifier checks the projected recurrence through 35 consecutive
Krylov observations. Since the series already obeys p, Cayley-Hamilton
extends that check to all r. Its leading 12-by-12 Hankel determinant is
one modulo three, proving minimality. Modulo two all b_r vanish, by the
eighteen even initial values and the integral recurrence.

The checks X^243=X mod g and gcd(g,X^3-X)=1 prove g irreducible, since its
degree is prime. Thus F3[X]/g is literally F243. Its element alpha=X mod g
has exact multiplicative order 121: alpha^121=1 but alpha^11!=1. The degree
gives Frobenius order five, while multiplication by alpha gives order 121.
These are different operations. The arithmetic source is

    3^5-1 = 242 = 2*11^2.

Therefore alpha generates the entire norm-one subgroup of F243^* over F3.
This is a concrete occurrence of 5 and 11 in a native tournament response,
with a field map, not an identification by equal cardinalities.

The actual scalar sequence has exact period 726=2*3*121 modulo three,
with no transient. E9 explains the factors: the -1 sector has order six,
and the doubled g sector has order 363. A direct recurrence orbit returns
to the initial eighteen-value vector for the first time at 726. Modulo
six the period is also 726; modulo two it is one. These are all-r statements
certified by a repeated full recurrence state, not guesses from scalar hits.

## E5. The nilpotent coordinate is observable even when the period agrees

The minimal ternary response module is

    F3[X]/(X+1)^2  +  F3[X]/g^2.                        (E10)

The second summand is F243[epsilon]/epsilon^2. Explicitly, in F3[X]/g^2,

    epsilon=g(X)/g'(X),        alpha=X-epsilon.          (E11)

The derivative is a unit; Taylor expansion gives g(alpha)=0. The receipt
checks epsilon!=0, epsilon^2=0, alpha^121=1, X^121!=1 and X^363=1.
This is the actual transport decomposition X=alpha+epsilon.

Cubing kills epsilon. It is consequently NOT an automorphism of this local
algebra, although cubing is Frobenius on the field subalgebra. A field-only
summary discards five observable ternary coordinates. Its loss is not
detected by the total clock period: quotienting only this epsilon layer
leaves dimension seven and still has clock order lcm(6,121)=726, because
the other doubled sector already contributes the factor three. The minimal
degree-twelve response proves that this quotient loses future observations.

The standard modular-group-algebra nilpotent mechanism is inherited; compare
[THM-2839, prime-power unit-mass spectrum and provenance](../../../01-canon/theorems/THM-2839-prime-power-unit-mass-full-spectrum-and-q11-response-provenance.md).
That theorem's physical-source restrictions remain; no LRC current is
transported here. The new operation is the explicit recurrence (E1)--(E11).

The recovered [66/70 lattice-content packet](../aliquot_66_five_fifteen_20261028/README.md)
gives sigma(66)=144 from {2,3,11}; the September
[Lucas/circulant note](../collatz_circulant_20260930_circulants_lucas_cubic_monotile.md)
records L5=11 and its separate finite clock constructions. Those are
different sources and operations. The present norm-one order 121 neither
proves a Collatz return nor identifies an aliquot orbit with this module.
It gives a new concrete object against which such proposed maps can be tested.

## E6. Reproduction and next targets

```text
python -X utf8 -B 05-knowledge/results/overnight_rank_return_20261029/periodic_defect.py
python -X utf8 -B -O 05-knowledge/results/overnight_rank_return_20261029/periodic_defect.py
```

The normal and optimized checks rebuild the finite transfer, its polynomial
identity, nonzero Hankel minors, Smith layers, independent boundary queries,
field quotient and exact periodic orbit. No gate relies on Python assertions.
The Smith factors concern the named observation lattice, not a full native
three-run Smith form. No global minimal attachment carrier is claimed.

Open next targets: find a path-level involution explaining the even counts;
extend E3 to additional periodic run words and allowed complement operations;
and determine when a realization-constrained family can discard nilpotent
response coordinates that arbitrary weighted states demonstrably need.
