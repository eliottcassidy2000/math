---
id: HYP-2924
title: LRC14 tournament-realizability summit atlas
status: PROOF-INTERFACE / exact finite atlas; not a proof of LRC14
source: codex-2026-06-22-S127
related:
  - HYP-2922
  - HYP-2923
  - HYP-2920
  - HYP-2921
  - HYP-2913
  - HYP-2914
  - HYP-2917
  - HYP-2907
  - HYP-2908
  - THM-201
  - THM-343
---

# HYP-2924: tournament-realizability summit atlas

This hypothesis turns the prompt "make any pairwise relation into a
tournament, then ask which isomorphism classes are achievable" into an exact
LRC14 interface.

## Exact carriers

The executable atlas is
`04-computation/lrc14_tournament_realizability_summit_codex_s127.py`, with
stored output
`05-knowledge/results/lrc14_tournament_realizability_summit_codex_s127.out`.

Carrier A is the **apex-clock tournament**:

- Vertices are seven selected points on `Z/14Z`.
- Pairwise observable is clockwise displacement.
- Switch/gauge: `x -> y` iff `0 < y-x < 7` on the `14`-clock.
- Diameter/collision ties are either excluded or broken by the Hamiltonian
  path of increasing residue.
- The quotient preserves the half-clock order class and diameter-tie count; it
  destroys actual speed divisibility, runner labels, and off-apex `M`.

Carrier B is the **terminal runner-phase tournament**:

- Vertices are the thirteen speeds of an AP/GW terminal row, ordered by
  increasing speed.
- Pairwise observable is the denominator-14 residue `s mod 14`.
- Switch/gauge is the same half-clock cutoff; collisions and diameters are
  tie-broken by increasing speed.
- The quotient preserves apex winding/tie class; it destroys q-threshold
  divisibility and off-apex escape size.

## Exact finite readout

For the seven-vertex `Z/14` apex clock:

- Strict no-diameter policy: `128` labelled rows, `10` tournament isomorphism
  classes, diameter-pair histogram `{0: 128}`.
- Tie-broken policy over all `7`-subsets of `Z/14`: `3432` labelled rows, still
  `10` isomorphism classes, with diameter-pair histogram
  `{0: 128, 1: 1344, 2: 1680, 3: 280}`.
- The regular carousel class appears in the tie-broken atlas with score
  histogram `((3,7),)`, `c3=14`, `scc=(7,)`, and `hp=175`.

For the terminal runner-phase carrier:

- AP and the residue-liar `{1,...,11,13,26}` are the same apex-winding class
  `T0`.  They have the same score histogram, cycle count, SCC, and Hamiltonian
  path count in this quotient, but differ arithmetically: AP has `q=14`,
  `M=1/14`, while the residue-liar has `q=12`, `M=1/12`.
- GW `12->24` is class `T1`, tight with `q=14`, `M=1/14`.
- Loose `8->16` is class `T2`, with `M=2/23`.
- Loose `10->20` is class `T3`, with `M=2/27`.  It shares the same coarse
  score histogram and directed-triangle count as GW, but differs by exact
  isomorphism/Hamiltonian-path count.
- The Farey near-miss `12->36` is class `T4`, with `M=3/41`.

## Proof readout

Tournament Analysis is useful here only after three choices are explicit:

1. Vertex set.
2. Pairwise observable and cutoff.
3. Tie Hamiltonian path.

The right proof shape is a state lift:

```text
bad LRC14 row
  -> exact achieved tournament/OCF class, with q/off-apex data retained
  -> forbidden tournament class or forbidden OCF component is absent.
```

The atlas also gives the main guardrail.  Tournament isomorphism classes are
not enough by themselves: AP and the residue-liar collapse to the same
apex-winding class.  Therefore any summit proof using tournament classes must
carry at least the q-threshold/divisibility data outside the tournament
quotient, and must carry off-apex witness data for near-miss classes such as
`12->36`.

Post-rebase integration: KPS S40, HYP-2923, and mac-mini S57 sharpen the same
warning.  KPS S40 identifies the AP apex winding tournament as a regular
rotational `13`-vertex tournament, while HYP-2923 retracts the stronger
H-max/H-extremal reading: at `13 == 1 mod 4` there is no Paley tournament, and
global Hamiltonian-path maximality for this class is unverified.  The surviving
bridge is necessary-only and magnitude-blind: higher or loose lifts can share
the regular apex class while escaping at q-threshold or Farey-neighbor scales,
and sink-measure is lonely measure rather than tightness.  Mac-mini S57 defines
the full winding realization over time and verifies that achievable classes
over all `t` do not distinguish tight from loose rows; the metric lives at the
optimum isomorphism class, which recovers the residue/three-gap census
HYP-2913.  Thus HYP-2924 should be used as the finite exact atlas for the
apex/terminal branch, while S57 supplies the broader over-time realization and
HYP-2923 supplies the corrected regular-apex / no-Farey-sink guardrail.

## Assumption challenge

The session considered runners, residues, gaps, fixed circle sections, section
boundaries, wall-crossing events, cover arcs, Fourier modes, matroid/cycle
circuits, and proof obligations as possible tournament vertices.  The selected
vertices are `Z/14` apex points and denominator-14 runner phases because they
preserve the exact half-clock cutoff where AP/GW bind.  The challenged
assumption is that raw runner tournaments suffice; they do not.  They are an
interface to a state-lift theorem, not the whole arithmetic proof.
