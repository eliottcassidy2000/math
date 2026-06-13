---
id: HYP-2187
status: OPEN synthesis from S617; small-n audits verified, constructive Royle bridge open
source: user-2026-06-03; codex-2026-06-03-S617
tags: [equinumerosity, equidecomposability, tournaments, even-graphs, Royle-even, Dehn-invariant, path-homology, quotient-fibers]
---

# HYP-2187: equinumerosity is cardinal shadow; equidecomposability is retained fiber data

Namespace note: while S617 was closing, incoming opus-S599v claimed HYP-2186
for the strong-component/scissors-volume synthesis. This file is HYP-2187 and
is complementary: it audits the Royle/even-graph count cautions and the
`(H,beta1,c3-packet)` tournament fiber refinement.

S617 separates two ideas that the repo had been letting touch too loosely.

Equinumerosity says two quotient worlds have the same number of objects.
Equidecomposability says two objects remain equivalent after retaining the
invariants needed by the target predicate.

For tournaments, the bridge should not be a bare bijection. It should be a
fiber functor that records which side-channel survives:

```text
tournament class
  -> cardinal shadow / A000568 count
  -> H-volume fiber
  -> (H, beta1) scissors fiber
  -> (H, beta1, odd-cycle packet) refined scissors fiber.
```

## Verified small-n audit

`04-computation/equinum_equidecomp_bridge_s617.py` verifies:

1. The tournament Burnside count `A000568` is not reproduced by the
   degree-even/Euler graph quotient except at `n=3` in the audited range.
   Counts for degree-even graph iso classes are:

```text
n=3: 2, n=4: 3, n=5: 7, n=6: 16
```

where tournament iso classes are:

```text
n=3: 2, n=4: 4, n=5: 12, n=6: 56.
```

2. The naive even-order Burnside complement is also not Royle-even:

```text
n=4: tournament=4, naive-even=7
n=5: tournament=12, naive-even=22
n=6: tournament=56, naive-even=100
```

Thus the Royle-even equinumerosity needs the intrinsic automorphism-orientation
graph property, not the cycle-space projection and not the Burnside parity
split.

3. Equidecomposability refines `H`. At `n=5`, `H=9` splits into two `beta1`
   classes. At `n=6`, five `H` values split by `beta1`:

```text
17, 23, 33, 37, 45.
```

4. Adding the disjoint directed-3-cycle packet polynomial refines the `n=6`
   quotient further:

```text
iso classes          56
H classes            19
(H,beta1) classes    24
(H,beta1,c3) classes 29
```

## Interpretation

`H` is volume-like: it is a scalar evaluation. `beta1` is the first
Dehn/scissors obstruction. The disjoint 3-cycle packet polynomial is an
early odd-cycle packet side-channel. Equal `H` therefore gives a cardinal or
volume fiber, not equidecomposition.

The useful future bridge is:

```text
Royle-even cardinal theorem + retained tournament fibers.
```

The Royle count may identify the right-sized graph universe, but a proof route
must still say what happens to `H`, `beta1`, score/cut data, and odd-cycle
packets under any proposed correspondence.

## Tournament Analysis

S617's Tournament Analysis uses quotient lenses as vertices rather than
runners, arcs, or graph vertices. The pairwise observable is a majority over:

```text
cardinal evidence,
constructive map,
H preservation,
scissors refinement,
cross-domain reach,
low risk.
```

The tournament is transitive. It ranks:

```text
full tournament isomorphism
> (H,beta1) scissors quotient
> H-volume quotient
> (H,beta1,c3) packet quotient
> degree-even cycle projection
> Royle-even cardinal theorem
> naive Burnside parity split.
```

The reading is not that full isomorphism is the proof target; it is too fine.
The working quotient is `(H,beta1)` plus packet side-channels when `n=6`-style
splitting matters.

## Assumption challenge

The tournament vertices are not automatically runners or arcs. S617 considered
Burnside cycle types, tournament iso classes, degree-even graph classes,
Royle-even graph classes, `H` fibers, `beta1` chains, disjoint `c3` packets,
OCF proof obligations, and quotient lenses.

The quotient must preserve the predicate being used. For equinumerosity, the
predicate is only "same count." For equidecomposability, the predicate is a
scissors class, so forgetting `beta1` or packet data destroys real information.

## Artifacts

- `04-computation/equinum_equidecomp_bridge_s617.py`
- `05-knowledge/results/equinum_equidecomp_bridge_s617.out`
- `07-reflections/equinumerosity-equidecomposability-fiber-bridge-s617.md`
