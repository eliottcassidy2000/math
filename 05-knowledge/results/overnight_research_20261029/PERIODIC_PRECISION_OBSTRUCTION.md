# An exact precision recurrence and the boundary of the Frobenius lift

2026-10-29. **PROVED all-height formulas using finite integer/p-adic
certificates; FINITE-EXACT realization tests.** Same-agent alternative-path
verification. This concerns the specified periodic native response and
its weighted observer language, not arbitrary tournament composition.

The source,35 boundary types, bulk matrix B, single terminal selector,
degree18 recurrence P, and integer observer O are inherited from
[PERIODIC_COMPLEMENT](../overnight_rank_return_20261029/PERIODIC_COMPLEMENT.md).
The successful residue-field lift is
[PERIODIC_FROBENIUS_LIFT](PERIODIC_FROBENIUS_LIFT.md): on the twelve-dimensional
ternary response it gives sigma T sigma^-1=T^487, of order five, with
the observation selector transported.

The inheritance portfolio is native completion / precision dynamics /
Frobenius realization. The live board is original ports; source weights;
clock image sizes; contracting coordinates; selectors; root carries.
The niche predecessor is
[THM-2727, fixed-rail-target-deck-frobenius-realization-no-go](../../../01-canon/theorems/THM-2727-fixed-rail-target-deck-frobenius-realization-no-go.md):
a formal coefficient map need not act on the corresponding physical bank.
Its LRC source is different. The common research move is to test the
actual operation and source diagram, not just formal clock identities.

## P1. All precision levels and all bulk powers

For h>=1 let

    W_h = image(O modulo 3^h).

This is the smallest carrier for arbitrary weights on the35 boundary
states under bulk insertion and the ONE supplied terminal observation.
Let T act by bulk insertion. For every r>=0,

    log_3 |T^r W_h|
       =12h+5(h-2)_+ +(h-2-2r)_+.                    (P1)

Here x_+=max(x,0). This counts actual image states at the given source
precision, not the quotient of the image lattice by its own3^h multiple.
It gives exact recurrences in both h and r. At r=0 it recovers
12h+6(h-2)_+, the previously known observer Smith count.

**Proof.** Let C be the18-by18 integral recurrence companion. The literal
source equation is C O=O B^T. Work over Z3 with ambient S=Z3^18 and
source lattice L=image(O). Its elementary-divisor exponents are twelve
zeros and six twos, by integer Smith computation.

P(0)=2304 has valuation2 and P'(0)=4352 is a unit. Successive digit
lifting therefore gives a unique root rho=0 mod3 in Z3, and v3(rho)=2.
Write P=(X-rho)G. Both G(0) and G(rho)=P'(rho) are units. Hence

    E=G(C)/G(rho)

is an integral projector onto a free contracting line; C acts there by
rho. Its complementary summand is free of rank17, and C and C^-1 are
integral polynomials on it. The source lattice L is stable under both
projectors, and its unit part is stable under C and C^-1.

The explicit E O certificate has every entry divisible by9 and an entry
not divisible by27. Thus the source lattice on the primitive contracting
line has exponent exactly2. Removing it from the known lattice profile
leaves twelve zeros and five twos on the unit summand. Applying C^r
changes only the contracting exponent, to2+2r. Reducing these18 lattice
directions modulo3^h gives P1. This proves all h,r; a finite grid is not
being extrapolated. Independent integer Smith checks at r=0,1,2,3 and
local checks through h=12, including r=487, support the symbolic argument.

## P2. No full precision-carrier conjugacy at modulus243

At h=5, P1 gives

    |T W_5|=3^76,        |T^487 W_5|=3^75.             (P2)

Any bijection conjugating T to T^487 would preserve these image sizes.
Thus **no such bijection exists**, even if it is nonlinear and unrelated
to the old formula for sigma. In fact every positive exponent a>=2 has
the same strict first-image discrepancy at h=5.

This is the first precision at which this particular image-cardinality
test separates T and T^487; it does not assert that a full compatible
lift exists at every earlier level. T already becomes noninvertible on
W_3, at modulus27. Consequently the successful F243 field chart and the
integer precision ring Z/243 have sharply different transport laws.

## P3. The actual initial trajectory has a different, sharp boundary

Let y_r=(b_r,...,b_(r+17)) be the actual response window. Its initial
boundary weights are retained and satisfy y_0=O a^T exactly. Its
contracting coordinate E y_0 has valuation4, rather than the lattice's
valuation2. Therefore its exact preperiod modulo3^h is

    tau_h = ceil((h-4)_+/2).                           (P3)

The unit component has exact period

    pi_h = 726*3^(h-1).                               (P4)

To prove P4, the full unit operator C^726 is I modulo3. The certificate
shows (C^726-I)y_unit has valuation exactly1, and the actual residue
trajectory has period726. If A=I+D with D divisible by3, the binomial
identity A^3-I=3D+3D^2+D^3 increases the valuation of (A-I)v by exactly
one. A power prime to3 leaves it unchanged. Thus the least additional
factor at precision h is3^(h-1). This proves P4 at every height.
Direct orbit enumeration verifies the initial four precision levels.
The nonunit coordinate vanishes precisely at P3; the complementary
operator is invertible, so no earlier state can lie on its eventual cycle.

On the actual orbit, the source-fixed map forced by a clock conjugacy is

    y_r -> y_(487r).                                  (P5)

Since487 is coprime to every pi_h, P5 is a bijection whenever tau_h is
zero or one. Thus the residue lift restricted to this actual trajectory
has compatible source-fixed extensions through h=6, even though the
full arbitrary-weight carrier already fails the image test at h=5.
These orbit bijections are not asserted linear on the full carrier.

At h=7 the preperiod is two. P5 then maps a transient point into the
cycle and collides with a later point. The literal scalar witness is

    b_1=38,     b_529255=1496 mod2187,
    y_(487*1)=y_(487*529255) mod2187.                   (P6)

Here529254=pi_7. The two entire source windows are distinct; their
reindexed images coincide. Any source-fixed conjugacy must obey P5,
so none exists. More generally tau_h>=2 for every h>=7, and acceleration
reduces its preperiod, excluding a bijection at every such level.

## P4. Boundary re-markings and pure bulk words also fail

The35 single-boundary-state columns of the old twelve-coordinate observer
have22 distinct values. Applying the lifted Frobenius takes27 of the35
columns outside this set. In particular it cannot be induced by ANY
permutation of these35 boundary types. The receipt supplies column0 as
a literal witness. This does not exclude arbitrary weighted sums or
operations outside the specified state language.

Nor is sigma a pure bulk word: every T^r commutes with T, whereas
sigma T sigma^-1=T^487 is different from T. A legal operation language
that centralizes a supplied clock cannot realize a normalizing map
with a different clock exponent. This criterion can be applied once the
actual generators and their commutation identities have been supplied.

## Consequences and limits

The source is an actual periodic tournament-complement family. The map
tested is a clock conjugacy preserving its original observation/source
diagram. The full-weight obstruction reads image cardinalities; the
physical-prefix obstruction reads a retained contracting amplitude.
The different valuations2 and4 explain the different boundaries243
and2187. Counting only stable periods misses both.

The same discipline is required for the native symplectic completion
and any future arithmetic Selmer bridge: specify whether the source is
an arbitrary weighted carrier, one physical orbit, or an arithmetic
subgroup, and retain its transition maps. These exact obstructions do
not prove Sha finiteness or exclude exceptional Collatz integers.

Reproduce with `python -B periodic_precision_obstruction.py` and -O.
All checks use exceptions. The receipt stores the original source hash,
compatible Hensel residues, projector and source amplitudes, integer Smith
checks, modular image checks, literal low-precision orbits and P6.
