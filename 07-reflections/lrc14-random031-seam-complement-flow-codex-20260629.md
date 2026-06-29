# LRC14 random031 seam-complement flow

HYP-3510 made one coarse version of the HYP-3482 experiment literal after
HYP-3484's seam-surgery packet and the upstream HYP-3486 fiber graph.  Delete
the two max-delta seam gates, keep the `282` q=`14V` witnesses, attach them to
branch-compatible survivor gates, and add branch-order edges.

The result is stronger than the initial metaphor but coarser than HYP-3486:
the pure witnesses are two branch cycles of `141` each, and the
branch-ordered incidence graph becomes one component after survivor-gate ports
are included:

```text
seam_complement_components_after_branch_order=1
seam_complement_summaries=[{'C':80,'G':122,'W':282}]
```

So the seam is not required for phase transport.  It is boundary debt.  The
six bypass rungs on components `54` and `43` are not the whole carrier; they
are the visible stitch points in a larger one-component complement graph.

The two reframes I would keep:

1. **Two-sheet zipper:** branch sheets first, survivor-gate stitching second,
   mirror completion third.  This is more exact than saying "the mirror zips
   the sheets"; the gate ports already fuse the complement before mirror edges
   are added.
2. **Puncture exact sequence:** the `40` no-gate witnesses are balanced
   `20/20` by branch and live on the branch order between gate ports.  They
   are not separate phase obstruction.  The owner-boundary punctures/seam are
   the residual debt.

This connects cleanly to the old `n+2` versus `n*2` split.  The `n*2` side is
now not just the ordered six-hit bypass blocks from HYP-3483; at coarse
incidence level it is the one-component seam-complement flow carrier, while
HYP-3486 refines that carrier into rank-2 routed cells, free-hole packets, and
one pure bypass component.  The `n+2` side remains the owner boundary
insertion: seam-only owners `(45,147,169,173)` plus the apex layer.

The incoming monad-explorer HYP-3486 note makes the free-hole side less
mysterious: the `14` packets should split into `10` individually bracketed
packets and two same-branch doublets, with both doublets still bracketed by
ordinary rank-2 packets.  That sounds like a finite bead lemma rather than a
new topological wall.

After rebasing over HYP-3490, the graph carrier looks less optional.  HYP-3490
shows the random031 pair-current lane is blocked by private E/branch-touched
labels, so the seam-complement graph is the remaining transport object rather
than another decorative quotient; HYP-3486 gives the proof-safe fiber
decomposition inside that object.

Next experiment: compare HYP-3510's connected incidence collapse with
HYP-3486's legal mirror-run components and free-hole doublet refinement on the
other HYP-3477 hard orbits.
Prediction: random031 is the unique "deleted seam, connected complement"
packet; the others discharge earlier by projection current and will not need
the two-sheet zipper sidecar.
