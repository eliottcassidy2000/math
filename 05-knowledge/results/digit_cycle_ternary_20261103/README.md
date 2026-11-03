# Decimal carries, ternary symmetry and a return at 2107

2026-11-03. **PROVED scoped / FINITE-EXACT / CITED background.** Same-agent
proof review, with separate exact algorithms and literal source controls.

**Continuation:** [the6325 packet](../return6325_completion_20261103/README.md)
now supplies both actual retained-name returns and next children18977/18979.
It also gives exact attachment-depth ranks and the marked clone obstruction.
The completed6325 target search should not be repeated.

The cycle `4,8,3,7,2,6,10,5,9,4` now has several explicit extensions.

* [DIGITS_AND_TOWERS](DIGITS_AND_TOWERS.md) proves a nearest-neighbour
  formula for s10(4n), the smallest decimal-append observer, and
  `sum_(n<10^L) z^(s10(4n))=(1+z^2+...+z^10)(1+z+...+z^9)^(L-1)`.
  Every fixed dyadic multiplier has a finite-range decimal carry formula.
  Adding halving legality makes every positive integer distinguishable.
* [MINIMAL_CLOCK](MINIMAL_CLOCK.md) gives the sharp linear dimensions
  `2,4,10,28,...`, with `r_next=3r-2`, for the ternary affine clock.
  Thus the 4-to-10 jump costs exactly its order-six doubling phase.
  The previous 27-coordinate radical contains the literal nine-clock.
* [SYMMETRY_AND_PORTS](SYMMETRY_AND_PORTS.md) embeds its 54-element group
  into the actual common cyclic-triple relation of the 81/27 pair.
  Oriented intersections have sizes 27 and 9. An infinite family on
  `4*3^(k-1)` vertices has identical labelled triples but automorphism
  orders differing by `3^(k-1)`. Exact shared-port orbit and regularity
  recurrences classify all attainable stabilizer orders in this family.
* That family also admits a native normalization into the old intrinsic
  halving interface: one arc repair per bottom triangle, followed by two
  forced pair contractions to the odd ternary core. Old names and fibres
  remain explicit; this concerns vertex order, not path-count division.
* [RETURN2107](RETURN2107.md) constructs the next marked return with
  **16 native moves, no old deletions and peak order 12**. Its old root
  and spare survive, its profile is `(2107,151,153,2109)`, and its next
  reset produces **6325=25*11*23** with full Q4 phases retained.

The even Collatz chart is explicit: on multiples of 2^a, use x/2 above
height a and 3x+2^a at height a. It is conjugate to the ordinary map.
At a=1 the route `6->20->10->32->16->8->4->2` retains every native step.
Neither this conjugacy nor the new finite return proves convergence or
uniform recursive tournament closure.

[RECOVERY](RECOVERY.md) records the concept board and prior sources;
[AUDIT](AUDIT.md) gives proof scope and finite universes;
[HANDOFF](HANDOFF.md) specifies the next decisive source/operation tests.

Reproduce from the repository root:

```text
python -X utf8 -B 05-knowledge/results/digit_cycle_ternary_20261103/reproduce.py
python -X utf8 -B -O 05-knowledge/results/digit_cycle_ternary_20261103/reproduce.py
```
