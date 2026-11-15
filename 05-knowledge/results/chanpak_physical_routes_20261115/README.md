# Physical residue crossings, and the return theorem an exact proof needs

This session constructs literal paths between tournament witnesses,
instead of merely placing unrelated witnesses on opposite sides of a
target. The constructions work for infinite families. They do not
yet attach to both witnesses for an arbitrary prescribed odd N, so
Chan–Pak remains open.

## Actual routes

- **The earlier dyadic modulus, two single-arc reversals.** A seven-site
  skeleton, inflated at an untouched port, has Hamilton counts
  9H+9Q, 9H+10Q, 9H+11Q. All graphs are strong, and each step changes
  one actual arc outside the port. A triangle followed by singleton
  padding makes Q=3·2^k and H=3, giving 9M+27 → 10M+27 → 11M+27
  in one residue for every k≥2. A different skeleton gives fixed
  step 2H and reads 4Q₃+9Q₂ plus a multiple of H.
  [Port routes](PORT_ROUTES.md).
- **A step fixed while the height grows.** Equal-Q module cores,
  changed by at most two literal arc reversals at a time, give
  52704·2^t+645 → 52704·2^t+651 → 52704·2^t+657 modulo 6 for every
  t≥0. A separate transitive reader parameter makes the modulus
  6(2^s−3), with the baseline independently made arbitrarily high.
  It can match a sublinear power scale but does not have the old
  dyadic shape. Equal-H cores furnish actual high-dyadic routes with
  the opposite limitation: modulus comparable to height.
  [Equal-Q and equal-H routes](EQUAL_Q.md).
- **Between the earlier 721 and 727 witnesses.** One exact route
  visits 721 → 723 → 725 → 727 with changes to respectively two,
  three, and one marked arcs. Moves of at most two marks and count
  step at most two cannot connect these endpoints at that width.
  [Physical local route](LOCAL_ROUTE.md).

## Why the bridges do not automatically place themselves

The high-dyadic seven-site reader cannot absorb a second independent
wall at any of its available ports: a nontrivial second module
makes the second step *strictly larger* than the first. This is an
exact two-slot response defect, not a search cutoff. Conversely the
fixed-step reader admits a complete target-first inverse, but the
literal ×2 and ×3 clocks are synchronized. At fixed H its possible
middles are a finite union of 4C·3^t+9B·2^t−H rays, with only
O_H(log X) counts up to X. [Synchronized port](SYNCHRONIZED_PORT.md).

The old interleaved rectangles illustrate a different obstruction.
All six reduced +2/+4 rectangles through spine width 20 were
independently recovered and counted. Their differing support
midpoints obstruct appending an exterior-only +8 control that keeps
both original steps. Allowing a moving modulus does give exact
four-count roads, but a full classification of their exterior
calibrations shows that none has step 3·2^k. At height N, routes
restricted to one apex-arc change of size o(sqrt(N)) have at most
two effective words even with neutral end padding.
[Exterior calibration](EXTERIOR_CALIBRATION.md) ·
[One-apex component bound](ONE_APEX.md).

## The number-theory synthesis: certify integer returns, not orbit mixing

Hasse and Artin constrain finite multiplicative orbits; the exact
affine model shows why that is different from an integer return.
It even gives a semigroup that reaches all residues at each prime
power individually but fails modulo 6. Ramanujan's 7 and 21 and
the class-number example motivate distinguishing local eligibility
from the actual component, without identifying the mechanisms.
Tunnell/BSD, Friedlander–Iwaniec, Hadamard orthogonality and
Furstenberg's ×2×3 problem each suggest a different precise test:
a signed certificate, correlations at the scale of the full
integer target, exact mixed pairing, and independence of the
two clocks. [Sourced comparisons](NUMBER_THEORY.md) ·
[CRT firewall](CRT_FIREWALL.md).

The resulting proof target is a finite graph of marked physical
interfaces with **integer return shifts**. If all directed moves
lift from every sufficiently high realization, a positive return
and gcd one of the cycle shifts give an entire upper interval from
one physical anchor. No reverse move and no individual unit step
are necessary. For the special moving modulus 2^s−3 there is also a
linear-width guard programming every **difference** residue while
neutral on both modular address readings. It does not turn a
contrast into an absolute count. [Directed return theorem](DIRECTED_RETURNS.md)
· [Odd contrast](ODD_CONTRAST.md).

## Replay

Run the `verify_*.py` files from this directory. To replay the
width-20 inverse census, compile `rectangles20.cpp` (C++17) and run
the resulting executable. The infinite claims use exact transfer
formulas and the stated lemmas, not extrapolation of those controls.
