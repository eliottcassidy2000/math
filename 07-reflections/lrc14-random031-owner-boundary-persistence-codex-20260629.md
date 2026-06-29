# LRC14 random031 owner-boundary persistence

HYP-3520 now has two layers in one executable certificate.

The seam-local layer says:

```text
phase-flow side:      u=2t, 282 witnesses, 12 lower-delta bypass hits
owner-boundary side:  seven-owner forbidden seam, four seam-only labels
```

The pure bypass is local and isolated in the legal horizontal+mirror graph:
`R01`, `12` cells, branch split `6/6`, hit components `(43,54)`, endpoint rank
`(2,)`, owners `(23,93,113)`, and `external_horizontal_ports=0`.  That last
zero matters.  The missing labels are not waiting in a neighboring ordinary
packet.  They are boundary charge attached to the deleted seam.

The owner-current matrix is the clean finite object:

```text
45:+1,147:+1,169:+1,173:+1
```

The component-persistence layer says the owner boundary is globally closed at
the current scaffold:

```text
79 legal seam-sheaf components
= 64 rank2_owner_persistent
 + 14 free_hole_bracket_persistent
 + 1 pure_bypass_owner_boundary
unresolved_owner_boundary_components=[]
```

Full seam debt is no longer a loose end.  The `14` full-debt components are
exactly HYP-3511's bracketed free-hole packets.  The two half-open doublets are
still the visible puncture seams:

```text
('R77','R76') -> (4,7,90,93)
('R53','R05') -> (35,38,59,62)
```

After those caps are separated, the owner-support graph on the remaining
components is connected: `65` vertices in one overlap carrier.  That is the
geometric signal from this run.  The complement carries one phase-flow
receiver, while the owner boundary carries one connected support word plus a
single pure bypass debtor.

The quotient failures agree across both layers.  Global owner flow is too
coarse because it introduces owner `55`.  Dead-island owners are too small
because they recover only `45` after subtracting bypass owners.  Raw `12`, raw
`282`, component pair `(43,54)`, delta-2 phase block, owner count, endpoint
rank, component size, and mirror closure all miss or mix terminal owner fields.

So the owner boundary is a cochain word.  The moving phase flow is the
`n*2`/two-adic zipper in the seam complement, while the stationary owner charge
is the `n+2` seam insertion.  That is a span, not a choice between recursions.

HYP-3513 sharpens the route-side warning: the private-firewall bit can be read
through existing colored axes, but full route purity still needs a sidecar.
For random031, HYP-3520 identifies that sidecar as the owner-current/seam word.

The HYP-3493 seam-sheaf table adds a second, more local canary.  Over its
`79` legal components the pure bypass still has word `PDPPOOO`, and the
free-hole components are owner-empty.  But owner count, endpoint rank, branch
balance, component size, and mirror closure each merge `R01` with discharged
context.  The legal small quotients are `flow_class`, `allowed_exit`,
`owner_union`, and `sheet_pgf_bucket`.  This is the same lesson in a smaller
language: compression is allowed only when the owner-boundary predicate stays
constant on the quotient fibers.

HYP-3522 then slices the same boundary word more finely: transport
`(23,93,113)`, branch bracket lift `(147,169)`, residual `(45,173)`.  I read
that as a refinement of HYP-3520, not a contradiction: persistence says which
coordinates cannot be forgotten, while filtration says which subcharges should
be proved by separate lemmas.

Assumption challenge: I considered runners, raw gates, arcs, u-fibers, fixed
circle sections, section boundaries, wall-crossing events, residues, cover
arcs, free-hole packets, owner rows, hard-component owner maps, relative-H1
boundary classes, persistence classes, and proof obligations.  The chosen
tournament vertices are quotient sidecars and owner-boundary proof carriers.
That preserves pure-bypass discharge and the three-exit right boundary while
intentionally forgetting raw runner order and scalar witness counts.

Next formal hook: make the random031 right boundary a three-exit theorem:

```text
rank-2 route
or free-hole bracket
or pure-bypass owner-boundary
```

Then enforce the quotient rule: owner count, seam-owner count, endpoint rank,
component size, and mirror closure are illegal terminal quotients unless the
seam owner word or the named bracket/bypass sidecar survives.
