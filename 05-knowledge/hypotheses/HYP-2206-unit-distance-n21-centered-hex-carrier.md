# HYP-2206: The Unit-Distance n=21 Row Avoids the H=21 Trap by Splitting Into a Spine Plus Centered-Hex Bulk

**Status:** OPEN, supported by S630 synthesis of existing S625/S626 carrier
rows and the incoming THM-409/HYP-2205 SC perspective theorem.

## Claim

The relevant `n=21` unit-distance statement is not "there is a tournament with
`H=21`."  It is a carrier split:

```text
u(21) = 57 = 20 spine edges + 37 tile/bulk edges,
```

where `37 = 3*3*(3+1)+1` is the next centered-hex shell number after
`1, 7, 19`.

So the exact `n=21` Moser-carrier row does the opposite of a forbidden
Hamiltonian-path collapse.  It keeps the unit Hamiltonian spine as one side
channel and places the surplus in a shell-shaped bulk channel.  The integer
`21` is the vertex count, not a legal tournament `H` value.

## Evidence

S630 reuses the carrier facts established in HYP-2202 and HYP-2203:

- triangular carrier at `n=21`: `47 = 20` spine `+ 27` bulk;
- exact/Moser carrier at `n=21`: `57 = 20` spine `+ 37` bulk;
- every checked Moser witness through the exact `n=21` row has a unit
  Hamiltonian spine.

The `37` bulk term is a centered-hex number:

```text
C_hex(3) = 3*3*(3+1)+1 = 37.
```

This gives a sharper recursive pattern than the raw edge count: the non-lattice
Moser lane beats the Eisenstein triangular carrier by adding shell/bulk
capacity while preserving the traceable unit spine.

## S634 Exact-Core Update

S634 completed the first next test below.  A Hamiltonian-path DP on the five
stored exact `n=21`, `57`-edge graph6 cores from S614 found all five traceable.
Therefore the split

```text
57 = 20 unit-spine edges + 37 bulk edges = 20 + C_hex(3)
```

is graph-real for those exact cores, not only visible in the displayed Moser
slab family.

The core profiles are not identical.  Core 1 has degree histogram `{5:18,
8:3}`, while cores 2-5 each have a degree-3 vertex and more varied
deletion-edge decks.  This suggests the next `n=22` test should enumerate
endpoint-compatible Hamiltonian spines and one-vertex ears from the exact
cores, rather than treating all `57`-edge cores as one raw scalar class.

## Relation To n=7

The `n=7` row is where two different phenomena first separate:

- the compact centered hexagon appears;
- fixed lexicographic unit-flip tournaments lose an all-unit directed
  Hamiltonian path;
- the unit graph itself remains traceable.

Thus `n=7` is a gauge flop and a shell event, not an intrinsic loss of unit
spines.  The `n=21` row then becomes the larger test: after the initial
Eisenstein/triangular simplicity has ended, the exact Moser carrier still keeps
the unit spine, and its bulk lands on the next centered-hex scale.

## Relation To H=7 And H=21

HYP-2204 already records why unit-distance appearances of `7` and `21` are
scalar-collapse warnings rather than literal tournament H-values.  HYP-2206
adds the `n=21` carrier reading:

```text
H=21 forbidden: a single OCF/Hamiltonian-path evaluation cannot equal 21.
unit n=21 legal: the geometric carrier has 21 vertices and 57 unit edges,
                 split into 20 spine + 37 shell/bulk.
```

The side channels are exactly what prevent the false identification.  The
unit-distance construction lives in an additive geometric semiring; tournament
`H` lives in the OCF/strong-component product semiring.

## Relation To SC Perspective Flips

THM-409/HYP-2205 is the primary SC theorem: `Anti(T)` is a coset over `Aut(T)`
and edge reversal induces a canonical involution on rooted perspectives.  S630
keeps a complementary exact audit through `n=7` and uses it as a guardrail:
self-complementarity supplies a conjugation/perspective carrier, not a loophole
for `H=7` or `H=21`.

## Testable Subhypotheses

1. Dense Moser-carrier extremizers preserve a unit Hamiltonian spine through
   the first non-lattice exact rows.
2. The surplus over the spine is governed by centered-hex shell capacities more
   often than the raw triangular-lattice lower bound suggests.
3. A true unit-spine flop, if it exists, must first appear in a non-lattice
   exact optimum whose bulk side channel can no longer be linearly threaded.
4. Any quotient that maps `n=21` directly to tournament `H=21` is too lossy:
   it has forgotten spine/bulk decomposition.

## Next Tests

1. **Done in S634:** audit the five exact `n=21`, `57`-edge graph6 cores
   mentioned in HYP-2203 and compare their bulk shells.  All five are
   traceable and split as `20+C_hex(3)`.
2. Add `bulk = u(n)-(n-1)` as a first-class statistic in the Moser beam ledger.
3. Test whether `bulk` hits centered-hex values at other exact or beam-optimal
   rows, especially near `n=19`, `21`, and `22`.
4. Enumerate endpoint-compatible ears for the five exact `n=21` cores as a
   focused route into the `n=22` `60/61` frontier.
