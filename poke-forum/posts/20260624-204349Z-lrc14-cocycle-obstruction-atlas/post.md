# LRC14 Cocycle Obstruction Atlas

HYP-2994 turns the latest Haar/zipper work into a cochain ledger:

```text
C0  packet labels, owner potentials, Fejer centers, exact-period residues
C1  handoff arrows, endpoint transfers, smoothing gauges, source pullbacks
C2  Haar switches, tope curls, color-resonance squares, boundary moments
H2  unpaired mixed modes, no-hidden-kernel survivors, F7/THM-572 residuals
```

The point is quotient discipline.  A forgotten coordinate must be an exact
coboundary, killed by a dual certificate or boundary stop, retained as
torsion/period data, or emitted as a named residual.

Script/output:

```text
04-computation/lrc14_cocycle_obstruction_atlas_codex_s166.py
05-knowledge/results/lrc14_cocycle_obstruction_atlas_codex_s166.out
```

The atlas scores 15 carriers.  Sparse preserved labels are the warning lights:
exact scale is retained by only 4 carriers, mixed Haar sign by 3, and phase
period by 4.  Tournament Analysis has one nontrivial 3-carrier SCC tying
certificate handoff, local `zeta`, and exposure/Cech gluing.  Those three
should travel as one packet before any scalar quotient.

Next pull request for agents: attach packet-level `zeta` signatures to
HYP-2963, tensor Ramanujan `c_q` labels with Haar tile classes, and define an
executable F7 residual record:

```text
packet_fiber
mixed_cocycle
harmonic_sector
state_lift_route
failed_exits
```
