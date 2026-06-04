# HYP-2197: LRC and Unit-Distance Frontiers Need Channel-Complete Ledgers, Not Raw Scalar Beams

**Status:** OPEN, supported by S624 channel-jackknife computation.

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

S624 adds `04-computation/lrc_unit_distance_channel_ledger_s624.py` and stores
`05-knowledge/results/lrc_unit_distance_channel_ledger_s624.out`.

### LRC Witness-Orbit Jackknife

The LRC computation uses the `C=2n-1` survival bitmask and then deletes one
`<2,-1>` witness orbit from the test. The damage is the number of residues that
become false survivors in the impaired quotient.

At `n=14`, `C=27`, the AP row, repo `Vstar=(1,...,11,13,24)`, and doubled AP
are all blocked at the `C` tick. Their shell-channel damage is identical:

| Orbit gcd | Orbit size | Released residues | Kill-depth histogram |
|-----------|------------|-------------------|----------------------|
| `1` | `18` | `18` | `{1: 18}` |
| `3` | `6` | `6` | `{1: 6}` |
| `9` | `2` | `2` | `{4: 2}` |

Thus the gcd-9 channel has tiny mass but high redundancy. Forgetting it releases
only two residues, yet those residues are covered four times in the intact
ledger, giving redundancy price `2*4=8`, larger than the gcd-3 channel's
`6*1=6`. This is the LRC analogue of a small but load-bearing unit-direction
packet.

Controls sharpen the point. At `C=15` and `C=21`, AP/Vstar show the same
wide-shallow plus narrow-redundant split. At prime `C=23`, AP has one witness
orbit, so the raw bitmask has no shell-channel choice to damage; side-channel
scarcity moves to lift/carry data. The `n=12` Vstar control has two survivors
and depth histogram `{0:2,1:18,2:2}`, showing that the jackknife also detects
already-surviving channels rather than only blocked rows.

### Unit-Distance Fingerprint

S624 recomputes a lightweight Moser-carrier jackknife at target `9`, width `30`.
The full shell reaches `18` edges. Dropping one antipodal Moser direction gives
loss histogram `{0:3, 2:6}`: three directions are replaceable at this size and
six are edge-critical. The replaceable directions have usage `0` in the full
cluster, while each critical direction has usage `3` and loss `2`, giving
usage-loss price `6`. Gain ceilings at the same target give:

| Gain ceiling | Best true edges |
|--------------|-----------------|
| `1` | `9` |
| `2` | `15` |
| `3` | `16` |
| `4` | `18` |

This matches the stored S622/S623 pattern at a smaller target: direction
support and high-gain packets are load-bearing side channels, not just beam
implementation details.

### Channel Tournament

S624 compares nine channel/proof-obligation vertices under two switches:
proof-relevance and scaling-value. Both tournaments are transitive with score
histogram `{0:1,...,8:1}`, zero directed `3`-cycles, singleton SCCs, and one
Hamiltonian path. But the gauges disagree on `19/36` edges.

The proof gauge ranks `channel-complete ledger`, `LRC owner/carry/pinch fiber`,
and `UD gain-threshold solver` at the top. The scaling gauge ranks `raw
survival bitmask`, `UD canonical orbit budget`, and `UD gain-threshold solver`
at the top. The edge flips are the result: proof relevance and scaling value
are different quotients, and a raw fast frontier can be excellent for scaling
while still too lossy for proof.

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
