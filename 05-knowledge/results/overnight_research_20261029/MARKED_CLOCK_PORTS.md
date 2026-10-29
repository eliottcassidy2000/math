# The fully named paired carrier for matched ternary cross joins

2026-10-29, second hourly continuation. **PROVED for the explicit integral
matrix family; FINITE-EXACT supporting audit.** This closes the named-port
question in [Eisenstein clock joins](EISENSTEIN_CLOCK_JOINS.md), without
asserting an arithmetic realization as a Selmer or Sha group.

## M1. Full observations distinguish every source tuple

Use exactly K(d,e,a,b), B(a,b), Q and the known order-three clock
T=diag(A,A) from E1. Assume d,e are nonzero integers and
t=de-Q(a,b) is nonzero. The source domain includes both signs, so global
sign reversal is a closed operation. Cross attachments add arbitrary
integer(u,v) to(a,b); the same clock-unit re-markings remain allowed.

At every nonsingular query observe the finite paired group

    D_K=Z^4/K Z^4,
    lambda_K([x],[y])=x^T K^(-1)y modulo Z,

together with ALL FOUR named generator classes [e1],...,[e4]. An
isomorphism of observations must preserve their names and the pairing.
At a singular query the observation is a distinct terminal symbol.
The comparison is equality of all future experiment traces from the
original source. Separate branch queries can start from the same source;
no decoder surviving a terminal query in a single supplied run is claimed.

**The complete minimal future equivalence is equality of(d,e,a,b).**
Thus these coordinates are a closed sufficient carrier and distinct
tuples can always be distinguished. This is a statement about which
states can be identified, not about the number of machine words needed
to encode an integer tuple.

### Constructive proof

The group order at a nonsingular probe is

    |D_K|=[de-Q(a,b)]^2.

Interpret the terminal symbol as determinant0 solely for the numerical
decoder. The nine axis probes in E4 recover z=de,a,b: the cubic
coefficients of the two quartic response polynomials are2(2a+b) and
2(a+2b), and the quadratic coefficient recovers t. Every coefficient
comes from five integer values along its axis, with their origin shared.

Now use the attachment(-a,-b). The resulting matrix is diag(dJ,eJ),
which is nonsingular, and the orders of its first and third named
classes give D=|d| and E=|e| respectively. This already recovers d,e
if the domain is restricted to positive sources, but does not justify
forgetting signs in a domain closed under global sign reversal.

Choose the integer L=2D+2E+|z|+1 and, from the original state, attach
(L-a,-b). This produces B(L,0) and t_L=z-L^2, with

    |t_L|>2 max(D,E),
    (K_L^(-1))12=-e/t_L,
    (K_L^(-1))34=-d/t_L.                              (M1)

The two displayed pairing residues distinguish +E from -E, and +D
from -D: their differences have absolute value strictly between0 and1.
Hence both signed factors are recovered. This uses at most eleven
branch observations, counting the nine shared-axis observations once.
All later chosen attachments are legal integer cross couplings and
remain nonsingular by the bounds just given.

Conversely(d,e,a,b) constructs K and all named data, and propagates
cross attachments, known clock-unit re-markings, and global sign.
This proves sufficiency and necessity. It does not allow arbitrary
internal block operations; the recovery genuinely uses only cross probes.

## M2. The induced clock is contragredient

The primal lattice clock satisfies T^T K T=K. On D_K the induced clock
is therefore

    tau([x])=[T^(-T)x],                                (M2)

since T^(-T)K=KT. It preserves the pairing because
T^(-1)K^(-1)T^(-T)=K^(-1). Its action on the named classes is completely
specified by the known integer matrix T^(-T); in this model it requires
no independent extra invariant once the four class names are retained.

Using T instead of T^(-T) on the dual presentation generally fails.
This is the same type distinction as the earlier functional-port
condition T^T P=PC, not an identification of those two actions.
An abstract Smith decomposition without the maps of the named classes
would again discard information needed by future attachments.

## Hostiles, bridge and next boundary

The old(d,e)=(1,36) versus(4,9), B=0 example has equal z and gcd,
so its entire UNMARKED paired-response histories agree. Its first
named port has orders1 and4. This identifies the exact lost coordinate.

A smaller sign hostile is(1,1,0,0) versus(-1,-1,0,0): both initially
have the trivial named paired group. Under the common attachment(5,0)
they have t=-24 and opposite first named pairing residues. Initial
pairings alone are therefore insufficient, even after including names.

The arithmetic analogy is precise but limited: a generator bank and
its transport map can carry data lost by an abstract subgroup order.
Incoming [actual two-saturation certificates](../session97_overnight_arithmetic_carriers_20261029/halving/phase2/README.md)
retain those maps for specific elliptic curves. We have NOT constructed
a map from these K matrices to those descent groups. Neither(M1) nor
the square-order classification bounds the prime support or exponents
of Sha. Those remain separate arithmetic obligations.

The next decisive test is native multi-run transport: can its allowed
operations implement a separating family of cross probes while retaining
original port names and clocks? The81/27 control, its single-run modulo9
separator and incoming two-run modulo3 ranks24/22 must remain in that
test. A complete arbitrary-integer probe family does not establish that
a native cut language has the same power. The newly supplied two-run
integral gluing makes three-run transport the next concrete native test.

## Reproduction and audit

`python -B marked_clock_ports.py` and the same with `-O` decode every
nonsingular source with d,e in{-6,...,-1,1,...,6}, a,b in{-2,...,2}:
3,512 states. They use exact Fraction Gaussian elimination for every
branch determinant and named pairing, including376 singular branch
values. The inverse-column denominator lcm independently gives named
port orders at the cancellation probe. The primal and dual clocks,
relation-lattice preservation and global sign operation are checked
literally. There are2,472 sources where the incorrectly reused primal
clock fails the induced pairing test. The all-integer proof is above;
the finite bounds are supporting controls, not its quantifiers.
