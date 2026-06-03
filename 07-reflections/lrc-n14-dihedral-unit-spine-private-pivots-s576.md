---
source: codex-2026-06-03-S576
status: proof-route sharpening + bounded finite audit
tags:
  - lonely-runner
  - n14
  - dihedral
  - unit-shells
  - private-pivots
  - composite-slack
  - endpoint-owners
  - tournament-analysis
---

# n=14 Dihedral Unit-Spine Private Pivots

The new useful move is to stop treating the thirteen moving runners at `n=14`
as a flat set. HYP-2091 and the user's simplex/polygon framing say the geometry
is already clean: `m=n-1=13` is odd, so the runner tournament lies on the
regular polygon side, with reflection/converse symmetry rather than antipodal
tie walls.

The arithmetic is where the obstruction lives.

For `n=14`,

```text
C=2n-1=27=3^3.
```

The unit-shell clock has nine reflected shells:

```text
1, 2, 4, 5, 7, 8, 10, 11, 13.
```

Any strict sub-edge row must hit all nine. Otherwise a missing antipodal unit
shell gives the old second-gap witness at `2/27`. That leaves only four
runners for everything else. This is a much more rigid picture than "search
thirteen speeds".

## The Spine

The canonical spine is:

```text
P=(1,2,4,5,7,8,10,11,13).
```

The script `04-computation/lrc_n14_dihedral_spine_pivots_s576.py` fixes this
spine and scans four slack runners among multiples of three through `42`.
Those slack runners are the composite layer: they can cover gcd-3/gcd-9 holes,
D denominators, reflected n-clock gates, endpoint blockers, or pair-sum
shields.

The quotient ledger has:

```text
D obligations: 12
U obligations: 9
N obligations: 7
```

The scan found:

```text
rows scanned: 1001
full D/U/N quotient covers: 531
below 1/14: 0
floor rows: 2
open-gap rows: 1
```

The two floor rows are exactly AP and `V*` in this canonical form:

```text
AP: slack (3,6,9,12), M=1/14 at t=1/14
V*: slack (3,6,9,24), M=1/14 at t=1/14
```

The first open-gap row is:

```text
(1,2,3,4,5,6,7,8,9,10,11,13,36)
M=3/41 at t=17/41
```

So in the normalized unit-spine model, the slack layer does not seem to be the
place where a subfloor counterexample can hide.

## What This Does Not Prove

The bounded spine scan does not prove `n=14`.

The missing part is not numerical; it is structural. A real speed set can use
large representatives of the unit shells, not just the smallest representatives
in `P`. Replacing a large unit representative by its shell minimum preserves
the U obligation but can change D obligations, N obligations, and pair-sum
pinches. So the open lemma is an exchange lemma:

```text
normalize each occupied unit shell to P, or show that the attempted
normalization creates a floor/pair-sum witness before it can hurt the proof.
```

That is the next proof object.

## How It Talks To HYP-2094

HYP-2094 gives a finite-looking topological shadow: `n=14` reduces to `64`
self-converse round classes if open=>round=>lonely and realization-independence
can be nailed down.

This session gives a complementary labelled shadow, now recorded as HYP-2096:

```text
9 unit-spine owners + 4 composite slack owners.
```

The two should be glued. A self-converse round class is too unlabelled for the
obligation proof. A D/U/N owner table without the round class is too local for
the Burnside/dihedral proof. The target is a fibre table:

```text
self-converse round class
  -> unit-shell owner table
  -> four-slack composite debt certificate
  -> endpoint-owner/pair-sum witness
```

That looks small enough to attack.

Rebase note: incoming HYP-2095 computed new all-zero staircase values through
`k=8` and records a near-regular `n=14` score sequence. That is useful signal
for the same polygon/simplex shadow, but it also reinforces the warning here:
unlabelled `H`, `c3`, or score data does not replace D/U/N owners, endpoint
owners, and unit-shell private pivots.

## Cascade Reading

The user's cascade language fits this exactly. Think of the proof as a product
of conditional clearances.

First factors:

```text
the nine unit-shell clearances U_a.
```

Residual product:

```text
the four slack runners carrying D/N/composite/endpoint debt.
```

The unit-shell factors are not optional in the strict sub-edge regime. A
counterexample has to clear all nine before the residual product even starts.
That is why the canonical spine is the right first normal form.

Level-4 hyper-lacunary speeds can still exist, but after this reduction their
freedom is constrained: lacunarity can only live inside unit-shell
representatives or the four-runner slack layer. That converts a vague speed
growth question into a finite-looking exchange problem plus a four-runner
composite ledger.

## Tournament Analysis

This pass deliberately did not put runners as tournament vertices. It used
slack rows over the fixed spine.

Pairwise observable:

```text
M, full D/U/N cover status, private quotient pivot count, slack-shell profile,
and max speed.
```

Switch:

```text
harder row wins: smaller M, then full cover, then fewer private pivots,
then smaller maximum speed.
```

The tournament fingerprint on the 24 hardest sampled full-cover rows is
transitive: no directed 3-cycles, singleton SCCs, one Hamiltonian path. That is
not disappointing. It says this quotient is a sorted proof ledger. The cyclic
behavior, if it exists, must be in the next layer where endpoint owners,
arbitrary unit representatives, or the `64` self-converse classes are attached.

## Assumption Challenge

I considered runners, gaps, fixed sections, section boundaries, residues, cover
arcs, Fourier modes, D/U/N obligations, unit shells, slack rows, pair-sum
blockers, and self-converse round classes as possible tournament vertices.

For this session the right vertices were slack rows over a fixed unit spine.
That preserves the predicate "does the normalized strict-sub-edge candidate
cover every quotient obligation, and what is its exact maximin?" It destroys
realization-dependence of large unit representatives and endpoint ownership.

The challenged assumption is that one global tournament quotient should solve
the proof. Here the proof really has two faces: the polygon outside forces the
unit spine, while the simplex mesh returns as labelled ownership inside the
slack layer.

## Next Lemma

The concrete next lemma is:

```text
Unit-spine exchange lemma for n=14.

If S has one runner in each unit shell modulo 27 and M(S)<2/27, then either
S has M(S)>=1/14 by an n-clock or pair-sum witness, or S can be transformed
to the canonical unit spine P plus four nonunit slack runners without
decreasing M below 1/14 and without losing the private U pivots needed by the
D/U/N ledger.
```

Proving that would leave only the four-runner slack layer, where this audit
found no subfloor full-cover row through slack `<=42` and only AP/V* at the
floor.
