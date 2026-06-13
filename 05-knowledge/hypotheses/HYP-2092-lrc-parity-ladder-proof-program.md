---
id: HYP-2092
status: SUPPORTED by S576 proof-program audit; lemma proofs open
source: codex-2026-06-03-S576
related:
  - HYP-2091
  - HYP-2090
  - HYP-2089
  - HYP-2088
  - HYP-2087
  - HYP-2086
  - HYP-2083
  - HYP-2082
  - THM-397
---

# HYP-2092: HYP-2091 turns the LRC proof into a parity-ladder proof program

HYP-2091 should not be used only as a geometric slogan. It gives a proof
routing table.

For LRC parameter `n`, put `m=n-1` and `C=2n-1`.

```text
even n -> m odd  -> clean polygon tournament, no antipodal tie burden
odd n  -> m even -> wall/mesh tournament, antipodal tie resolutions required
```

Overlaying the clock work from HYP-2083, HYP-2088, and HYP-2090 gives four
cases:

```text
clean unit lane      : no tie burden, no nonunit C-shell burden
clean composite lane : no tie burden, but gcd-stratum descent is needed
wall unit lane       : tie-discharge is needed, no nonunit C-shell burden
wall composite lane  : tie-discharge plus gcd-stratum descent
```

The conjectural proof program is:

```text
L1. clean seam:
    attach D/U/N labels to converse-merged round nodes.

L2. clean exchange:
    every full D/U/N cover on the clean lane has a floor witness or a
    pair-sum witness, unless it descends by a private-pivot exchange.

L3. nonunit descent:
    every gcd-stratum defect for composite C either lifts to a second clock or
    descends to a smaller/neighboring labelled cover.

L4. wall discharge:
    each antipodal tie resolution on the odd-LRC wall lane either creates an
    observer source gap or reduces to a clean-ladder labelled cover.

L5. endpoint compatibility:
    THM-397 endpoint blockers are exactly the owner labels that stop pair-pinch
    escape clocks, so pair-sum covers cannot be treated as unlabelled covers.
```

## Evidence

`04-computation/lrc_parity_ladder_proof_program_s576.py` audits `n=4..18`.
It records the HYP-2091 parity lane, `C=2n-1` factorization, antipodal tie
pairs, labelled tie choices, unit/nonunit summand shells, D/U/N obligation
counts, converse-merged round nodes, and a proof-route label.

The n=14 line is the payoff:

```text
n=14, m=13, C=27=3^3
lane = clean_polygon
tie_pairs = 0
unit/nonunit shells = 9/4
D/U/N obligations = 12/9/13, total 34
round/converse seam = 316 round, 64 self-converse, 190 merged nodes
```

So the fourteen-runner hard row is not a wall-resolution problem. HYP-2091
says its geometry is already on the clean even-`n` ladder. The remaining burden
is composite-clock arithmetic: the gcd-3/gcd-9 nonunit shells of `C=27`, plus
owner-labelled D/U/N private pivots and pair-sum endpoint owners.

This is a sharper target than "classify all speed sets":

```text
clean n=14 geometry
-> attach D/U/N labels to 190 converse-merged round nodes
-> prove nonunit gcd descent or full-cover exchange
-> extract floor or pair-sum witness
```

## Which clocks matter

The proof should keep:

```text
pair-sum clocks     : exact 1D maximin candidates and THM-397 endpoint owners
n-clock             : floor/tight equality and observer-source escape attempts
2n-1 unit clock     : antipodal summand witnesses visible to multiplication
D clocks            : small denominator divisibility witnesses
labelled events     : runner owner, endpoint owner, pair-sum denominator
```

The proof can usually ignore, until labels return:

```text
primitive reset length              : normalized integer rows all close at period 1
binary cyclic lonely-word stabilizer : usually trivial in S571/S573
unlabelled round class alone         : loses D/U/N and endpoint ownership
runner-vertex tournament alone       : loses proof-obligation incidence
```

## Proof-shape consequence

The open LRC body should be viewed as:

```text
round labelled words
-> converse-merged round nodes
-> boundary/self-converse seam
-> D/U/N owner-labelled fibres
-> pair-sum endpoint owner fibres
```

A counterexample would have to survive every stage while preserving the labels
that the previous quotients forget. On the clean lane this becomes a
private-pivot exchange theorem. On the wall lane it becomes a tie-discharge
theorem.

## Tournament Analysis

Vertices in S576 are LRC `n`-ladder rows, not runners.

Pairwise observable:

```text
(wall tie pairs, nonunit shell count, D/U/N obligations,
 converse-merged round nodes, extended dihedral size)
```

Switch:

```text
orient toward the larger remaining proof burden.
```

Fingerprint:

```text
score_hist={0:1,...,14:1}
directed_3_cycles=0
sccs=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1]
```

The tournament is transitive because this pass is a routing table, not a chaos
detector. The useful cyclic structure is expected only after adding owner
labels inside a fixed row.

## Assumption Challenge

Possible tournament vertices considered:

```text
runners, LRC n-ladder rows, antipodal shell strata, D/U/N obligations,
converse-merged round nodes, endpoint owners, pair-sum cells, tie diagonals.
```

Chosen quotient:

```text
LRC n-ladder rows, with proof-burden vector labels.
```

Predicate preserved:

```text
which proof lemma family must be applied before the LRC witness can be forced.
```

Information destroyed:

```text
exact runner ownership, exact active pair-sum geometry, and exact endpoint
protection cores.
```

Challenged assumption:

```text
that a single runner-vertex tournament is the right object for the proof.
```

For this stage it is not. Runner vertices return only after the row has been
sent to the correct proof lane and labelled by obligations.

## Files

- `04-computation/lrc_parity_ladder_proof_program_s576.py`
- `05-knowledge/results/lrc_parity_ladder_proof_program_s576.out`
- `07-reflections/lrc-parity-ladder-proof-program-s576.md`
