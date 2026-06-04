# HYP-2197: LRC and Unit-Distance Frontiers Need Channel-Complete Ledgers, Not Raw Scalar Beams

**Status:** OPEN, supported as an S624 bridge program.

**Claim.** The productive common object between the LRC `n=14` program and the
unit-distance `n=22` frontier is a channel-complete frontier ledger. A raw
frontier scalar is too lossy in both problems: LRC needs survival masks plus
owner/carry/pinch/CRT side channels, while unit distance needs edge gains plus
direction support, compactness, canonical orbit budget, deletion resilience,
and geometry/obstruction labels.

This hypothesis is deliberately a bridge rather than a new proof. It asks for a
shared impairment protocol:

- remove or coarsen one LRC witness-residue/shell channel and measure false
  survivor creation;
- remove or coarsen one unit-distance direction/gain/orbit channel and measure
  edge loss or spurious extension;
- compare the damage ledgers by Tournament Analysis over channels and proof
  obligations, not over runners or point sets.

## Known Parent Signals

LRC parent signal:

- HYP-2187/S599w-x makes the LRC discrete witness state-local: a speed `v`
  contributes a survival bitmask over `Z/(2n-1)`, and adding speeds is bitwise
  intersection.
- HYP-2164 through HYP-2175 show that this scalar quotient is not enough for
  `n=14`; owner routes, carry cocycles, pinch certificates, CRT/lift data, and
  Cprime windows are load-bearing side channels.

Unit-distance parent signal:

- HYP-2188/S617 makes the Moser-carrier edge count state-local by frontier
  gains.
- HYP-2194/S622 and HYP-2196/S623 show that small unit-distance carriers expose
  load-bearing side channels under impairment: direction support, high-gain
  packets, width, compactness, and canonical orbit budget.

## S624 Evidence Summary

The S624 computation should supply the actual bridge table. The intended
minimum evidence is:

1. an LRC witness-orbit jackknife over `C=2n-1` survival masks;
2. a unit-distance impairment summary imported from the S622/S623 Moser lane or
   recomputed on a small carrier;
3. a single Tournament Analysis comparing channel types under at least two
   gauges, with edge flips and Hamiltonian-path counts reported;
4. an assumption challenge recording which predicate survives the quotient and
   what information is destroyed.

## Assumption Challenge

Candidate Tournament Analysis vertices include LRC witness residues, witness
orbits, shell/gcd strata, owner routes, carry cocycles, unit-distance
directions, gain packets, deletion cores, canonical orbit classes, obstruction
filters, and proof obligations. S624 will not assume vertices are runners,
arcs, or point sets.

The quotient should preserve the predicate "does this channel impairment change
the frontier decision?" It destroys exact continuous LRC time geometry, exact
unit-distance embedding data outside the tested carrier, and any full proof
certificate not represented by the retained side channels. The challenged
assumption is that frontier-gain beams scale by widening alone.
