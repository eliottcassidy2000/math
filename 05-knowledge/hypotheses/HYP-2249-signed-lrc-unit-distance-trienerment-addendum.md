---
id: HYP-2249
status: OPEN method hypothesis with exact finite gauge evidence
source: codex-2026-06-05-S674b
related:
  - HYP-2248
  - HYP-2241
  - HYP-2240
  - HYP-2202
  - HYP-2201
  - HYP-2132
  - THM-401
---

# HYP-2249 Addendum: Signed LRC Gauges and Unit-Distance Trienerment Edges

## Addendum Claim

The LRC and unit-distance tournament maps should keep one extra equality/gauge
coordinate before taking a binary tournament quotient.

For LRC:

```text
v_i -> +/- v_i
```

is an observer-gauge symmetry.  It leaves every loneliness value
`||v_i t||` unchanged, but it changes the runner-runner pair-clock address.
In particular, alternating signs turn opposite-sign pair clocks from
differences into sums.

For unit distance:

```text
distance < 1, distance = 1, distance > 1
```

should be treated as a three-state trienerment object.  The equality layer
`distance = 1` is the wall/tie layer that may carry a unit Hamiltonian spine.
Collapsing it immediately into either short or long erases exactly the datum
needed to ask whether the mandatory tournament Hamiltonian path is geometric
or an artifact of the chosen order.

## S674b Evidence

S674b adds
`04-computation/signed_lrc_unit_distance_trienerment_s674b.py` and stores output
in `05-knowledge/results/signed_lrc_unit_distance_trienerment_s674b.out`.

### Signed LRC Gauge

The script samples `756` times for total `n=14`, so the moving-runner count is
`13` and the shell modulus is:

```text
C = 2n - 1 = 27.
```

Families tested:

```text
AP13   = (1,2,...,13)
Vstar  = (1,2,3,4,5,6,7,8,9,10,11,13,24)
2AP13  = (2,4,...,26)
```

For all three families and all four sign patterns
`all_plus`, `all_minus`, `alternating_plus`, `alternating_minus`, the observer
clearance sequence is exactly invariant:

```text
observer_invariant=True
```

But the pair tournaments change sharply.  For AP13:

| sign pattern | unique tournaments | unique fingerprints |
|---|---:|---:|
| all plus | `100` | `18` |
| alternating plus | `355` | `29` |

The all-plus versus alternating-plus fingerprints differ at `588/756` sampled
times.

The key arithmetic readout is the pair clock split:

```text
same-sign differences = {2:11,4:9,6:7,8:5,10:3,12:1}
opposite-sign sums    = {3:1,5:2,7:3,9:4,11:5,13:6,
                         15:6,17:5,19:4,21:3,23:2,25:1}
```

So alternating AP13 exposes the odd interior pair sums

```text
3,5,7,...,25
```

inside the exact LRC shell `C=27`.  The observer projection cannot see this
address.  The runner-runner tournament can.

This is the constructive version of the earlier pair-sum modulus idea: the
multi-sieve additive face is not a new observer predicate.  It is a signed
pair-clock side channel living over the same observer shadow.

### Unit-Distance Trienerment

The unit-distance probe orients three distance states:

```text
S = short, distance < 1, follows base order
U = unit,  distance = 1, tie/trienerment layer
L = long,  distance > 1, reverses base order
```

Then it resolves `U` either positive or negative.

Small exact point sets:

| point set | n | pair states | base order |
|---|---:|---|---|
| square + center | `5` | `S=4,U=4,L=2` | no unit spine |
| Eisenstein hex patch | `7` | `U=12,L=9` | unit spine |
| Eisenstein radius-2 patch | `19` | `U=42,L=129` | unit spine |
| Eisenstein `19+2` fringe toy | `21` | `U=44,L=166` | unit spine |

The square-center toy is the impairment: the center is short from all corners
and isolated from the unit graph, so a unit Hamiltonian spine is impossible.
Any tournament Hamiltonian path must pay a short/long address.

The Eisenstein rows are the opposite signal: the equality layer carries a unit
spine through the `n=21` fringe toy.  With the unit-positive resolution the
tournaments are strongly connected on the larger rows; with the unit-negative
resolution and no short edges, they collapse to the reverse transitive order.

For the `n=21` fringe toy:

```text
pair states: U=44, L=166
base path states: U=20
unit_positive: c3=157, SCC=(21,)
unit_negative: c3=0, SCC=(1,1,...,1)
```

This does not classify true optimal unit-distance sets.  It gives the finite
diagnostic: before asking whether the Hamiltonian path has flopped from unit
to nonunit, keep `S/U/L` and ask which layer the path's edges actually occupy.

## Tournament Analysis

Vertices are proof lenses:

```text
signed_observer_gauge
alternating_pair_sum_clock
short_unit_long_trienerment
unit_spine_order
owner_address_derivative
raw_phase_or_edge_count
```

Pairwise observable:

```text
predicate preservation,
side-channel exposure,
exact finite evidence,
transfer to LRC14 / unit-distance n=21,
overclaim safety
```

The lens tournament is transitive:

```text
score_hist = {0:1,1:1,2:1,3:1,4:1,5:1}
c3 = 0
H = 1
```

Outscore order:

```text
signed_observer_gauge
> short_unit_long_trienerment
> alternating_pair_sum_clock
> owner_address_derivative
> unit_spine_order
> raw_phase_or_edge_count
```

## Transfer To LRC14

The signed-speed lesson is not "negative speeds make a new LRC problem."
Independent sign choices are invisible to the observer.  The lesson is:

```text
same observer shadow + different signed pair-clock address
```

This is exactly the address-coordinate theme from HYP-2240/HYP-2241.  A future
LRC14 proof can use alternating signs as a stress test for a proposed quotient:
if the quotient claims to retain pair-sum/shell information but cannot
distinguish all-plus AP13 from alternating AP13, it has probably discarded the
wrong address.

The AP13 alternating split also gives a clean target for THM-401:

```text
C=27 is the shell modulus;
opposite-sign pair clocks enumerate the odd interior sums;
the missing boundary shells are proof-address debt.
```

## Transfer To Unit Distance

HYP-2202 separated graph-level unit Hamiltonian paths from canonical
tiling-order flops.  HYP-2249 adds a stricter diagnostic:

```text
record the whole S/U/L word along each chosen Hamiltonian path.
```

For small Eisenstein rows, the unit word stays pure along a constructed spine.
For the square-center impairment, short edges destroy the unit-spine premise at
`n=5`.  This suggests a concrete next search:

1. Generate candidate optimal or near-optimal point sets.
2. Classify every pair as `S/U/L`.
3. Search for unit-spine Hamiltonian paths.
4. If none exist, record the minimum number and placement of nonunit
   impurities along any Hamiltonian path.

The "flop" should be measured by impurity profiles, not by a binary tournament
edge flip alone.

## Open Questions

- Can the alternating AP13 signed pair-clock split be turned into a formal
  pair-sum-sieve lemma over `C=2n-1`?
- In LRC14 owner/carry fibers, which sign gauges preserve the observer shadow
  while changing private-owner or D/U/N obligation ledgers?
- For unit-distance optimum candidates at `n=21` and beyond, is there always a
  unit Hamiltonian spine, or does the first short-edge impairment force a
  mixed `S/U/L` Hamiltonian word?
- Is the useful tournament vertex set points, distance states, path
  obligations, or deletion-owner records?  S674b says path obligations and
  distance states are better than raw points for this question.
