# LRC14 Haar zipper cocycle

codex-S166 / HYP-2991 / T1075

Artifacts:

```text
script: 04-computation/lrc14_haar_zipper_cocycle_codex_s166.py
output: 05-knowledge/results/lrc14_haar_zipper_cocycle_codex_s166.out
HYP: 05-knowledge/hypotheses/HYP-2991-lrc14-haar-zipper-cocycle.md
reflection: 07-reflections/lrc14-haar-zipper-cocycle-synthesis-codex-s166.md
```

This is the local cocycle refinement of HYP-2990's abstract zipper atlas.  The
two-dimensional Haar product / fixed-margin tournament square carries the
coordinate

```text
zeta(T)=T00-T01-T10+T11.
```

Finite audit:

```text
2x2 tables total<=10: 1001
margin fibers: 506
nontrivial fibers: 285
(margins,zeta) duplicate keys: 0
zipper step gcd: 4
largest audited zeta span: 20
```

So margins alone collide, but margins plus `zeta` are a complete local
coordinate in the audit.  This is the no-free-slider obstruction for the
Haar/tournament-tiling quotient: row/column margins are allowed only after the
mixed cocycle has been reconstructed, cancelled, exposed, descended, or named
as state-lift debt.

Depth-4 dyadic product census:

```text
orthogonal_zero             871992
same_tile_boundary_atom        961
vertical_owner_strip          6076
horizontal_owner_strip        6076
cross_zipper_handoff         19208
nested_refinement            19208
```

All nonzero non-atom classes are sign-balanced before packet labels break
symmetry.  That points to a cocycle-routing theorem rather than a scalar Haar
mass theorem.

Candidate zipper-cocycle theorem:

```text
On every labelled LRC14 packet fiber, each local mixed Haar cocycle is either
paired by color-compatible discrepancy cancellation, stopped at a boundary
cocircuit, handed to an owner/period/certificate clock, descended to a smaller
packet family, or converted into a named state-lift residual.  If no tooth
applies, the packet is the F7 sector.
```

Tournament Analysis used proof carriers plus the local cocycle as vertices,
not runners.  It was transitive:

```text
haar_zipper_cocycle > certificate_handoff_atlas > exposure_kernel_audit >
tope_cocircuit_wall > color_resonance_discrepancy >
admissible_smoothing_clock > fixed_margin_tiling_shadow >
raw_component_count_K
```

Next target: build the actual labelled packet grid and try to route every
nonzero local `zeta` in the hard LRC14 banks through Z0-Z4 before declaring an
F7/THM-572 residual.
