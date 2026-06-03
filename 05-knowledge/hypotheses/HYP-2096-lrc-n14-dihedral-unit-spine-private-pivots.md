---
id: HYP-2096
status: SUPPORTED by bounded canonical n=14 spine audit; normalization/exchange lemma open
source: codex-2026-06-03-S576
related:
  - HYP-2101
  - HYP-2100
  - HYP-2099
  - HYP-2095
  - HYP-2094
  - HYP-2093
  - HYP-2092
  - HYP-2091
  - HYP-2090
  - HYP-2088
  - HYP-2087
  - HYP-2084
  - THM-397
---

# HYP-2096: n=14 reduces to a dihedral unit spine plus four composite slack runners

The clean even-`n` polygon ladder should be sharpened at `n=14` by separating
the mandatory unit-shell spine from the residual composite slack layer.

For `n=14`, put `C=2n-1=27`. The strict sub-edge regime inherits the
unit-shell gate from HYP-2084/HYP-2088:

```text
U_a,  a in {1,2,4,5,7,8,10,11,13}.
```

So a putative counterexample with `M(S)<2/27` must spend at least nine of the
thirteen moving runners on the nine unit antipodal shells modulo `27`. The
remaining four runners are the only place where the nonunit `gcd 3/gcd 9`
holes, D/N denominator gates, THM-397 endpoint blockers, and pair-sum escape
debt can hide.

## Canonical Spine Probe

The S576 audit fixes the smallest representative of each unit shell:

```text
P = (1,2,4,5,7,8,10,11,13)
```

and scans four slack runners among multiples of `3` through `42`, using the
dihedral quotient of the D/U/N obligation ledger:

```text
D: 12 fixed denominator gates
U:  9 unit shells
N:  7 reflected n-clock gates
total quotient obligations: 28
```

Among the `1001` canonical spine rows in this bounded scan:

```text
full D/U/N quotient covers: 531
below floor 1/14: 0
floor rows: 2
open-gap rows between 1/14 and 2/27: 1
edge-or-above rows: 528
```

The only floor rows are the known AP and `V*` floors:

```text
AP slack  = (3,6,9,12), shells (3,6,9,12), M=1/14
V* slack  = (3,6,9,24), shells (3,3,6,9), M=1/14
```

The first full-cover open-gap row in the stored scan is:

```text
slack=(3,6,9,36), shells=(3,6,9,9), M=3/41, witness=17/41.
```

This says that the four-runner slack layer is not producing subfloor rows once
the unit spine is normalized to its least representatives.

## Proof Program

This is not yet an `n=14` proof. The scan does not cover arbitrary large
representatives of the nine unit shells. The missing lemma is a normalization
or exchange theorem:

```text
Given one runner in each unit shell modulo 27, lower or exchange the unit-shell
representatives to the canonical spine P while preserving the private U pivots
or creating an explicit floor/pair-sum witness.
```

If that lemma is proved, the remaining theorem becomes a four-runner composite
slack problem:

```text
P plus four nonunit/zero-residue slack runners
-> D/U/N full cover or immediate quotient-clock witness
-> THM-397 endpoint owner ledger for pair-sum blockers
-> gcd-3/gcd-9 descent or floor/open-gap witness
```

This complements HYP-2094. HYP-2094 collapses the even-ladder worry set to the
`64` self-converse round classes at `n=14`; HYP-2096 keeps labelled D/U/N
owners and private pivots attached to the clean polygon face. The two lenses
should meet when each self-converse class is assigned a unit-spine owner table
and a four-slack composite debt certificate.

S578/HYP-2100 sharpens the first step: before lowering a large unit-shell
representative, apply the HYP-2095 unblocked-small-pair sieve. In the bounded
one-unit-lift and named local/extreme stress audits, every full lifted unit row
already has such a cheap witness, leaving no no-cheap exchange residual.

S579/HYP-2101 recasts this as an apex-lift certificate sheaf: denominator-14
cheap witnesses are the central chart, side-denominator witnesses are side
charts, and ledger-failure restrictions are legal local sections.  The bounded
site has no residual defects under one-shell restriction.

## Why Doubled Prime n=14 Is Special

The doubled prime `14=2*7` puts the LRC runner tournament on the clean even-`n`
polygon side of HYP-2091: `m=n-1=13` is odd, so there are no antipodal
tie-resolution walls in the runner tournament. But the summand clock is not
prime:

```text
C=27=3^3.
```

That is the point of the doubled-prime wall here. Geometry is clean, but the
clock has nonunit shells. The proof is forced to be a labelled product of
conditional clearances: nine unit-shell factors first, then a four-runner
composite slack product that must carry all remaining debt.

## Tournament Analysis

Vertices in the S576 script are slack rows over the fixed unit spine, not
runners.

Pair observable:

```text
(maximin M, full-cover flag, private quotient obligations,
 slack-shell multiplicity, maximum speed)
```

Switch:

```text
harder = smaller maximin, then full cover, then fewer private pivots,
then smaller maximum speed; ties use lexicographic slack order.
```

The hardest 24 sampled full-cover rows form a transitive tournament:

```text
directed_3_cycles=0
sccs=24 singleton SCCs
hamiltonian_paths=1
```

That transitivity is a warning: this quotient is a ledger and ranking tool. It
does not yet contain the cyclic contradiction. Cycles should be sought after
the arbitrary unit representatives and endpoint-owner fibres are reattached.

## Assumption Challenge

Possible vertices considered:

```text
runners, gaps, fixed circle sections, section boundaries, wall-crossing events,
residues, cover arcs, Fourier modes, D/U/N obligations, unit shells, slack
rows, pair-sum endpoint blockers, and self-converse round classes.
```

Chosen quotient:

```text
slack rows over the mandatory unit-shell polygon spine.
```

Predicate preserved:

```text
the n=14 strict-sub-edge need for all nine unit shells, exact maximin for the
canonical row, and full/failed D/U/N quotient-cover status.
```

Information destroyed:

```text
arbitrary large unit-shell representatives, realization-dependence inside the
64 self-converse round classes, endpoint-owner identities, and pair-sum active
set histories.
```

Challenged assumption:

```text
that runner vertices or unlabelled round classes alone are enough for the
n=14 proof.
```

They are not enough at this layer. The useful object is the labelled
unit-spine-plus-slack fibre: polygon outside first, simplex mesh labels second.

## Files

- `04-computation/lrc_n14_dihedral_spine_pivots_s576.py`
- `05-knowledge/results/lrc_n14_dihedral_spine_pivots_s576.out`
- `07-reflections/lrc-n14-dihedral-unit-spine-private-pivots-s576.md`
