# Arbitrarily many inserted vertices still need only two cut ports

2026-10-21. **PROVED** for consecutive balanced centre-flanking growth.
This extends Section4 of `BALANCED_FIBRE_GROUP_AND_PORT_DEFECT.md`; it does
not concern arbitrary placement of deleted vertices in arbitrary tournaments.

Let k>=2 and d>=1. Embed T_k in T_(k+d) by leaving labels below k unchanged,
sending the centre k to k+d, and sending labels above k to v+2d. Delete

    {k,...,k+d-1} union {k+d+1,...,k+2d}.

These are the d vertices on either side of the retained centre. The inherited
pointed-Boolean decoder proves this is an actual HP surjection.

For each finite parent coordinate, let E or H be the subset of the d new
bits in its lift. An infinity coordinate instead forces all d bits to1 for
its physical placement. All new left-shore vertices occur at the ends of
the pre- or post-centre A runs; all new right-shore vertices occur at the
beginnings of the pre- or post-centre B runs. Therefore all deleted vertices
occupy only the same two seams (7a) in the main note. Within an occupied seam,
the new A block is followed by the new B block, possibly with either empty.

**Theorem4.** The two-seam correction has size at most2, independently of d.
The pre seam is unoccupied exactly when E=empty and H=full; the post seam is
unoccupied exactly when E=full and H=empty. An occupied seam touching START
or END contributes no correction edge. In every other case the seam's
actual vertex pair contributes one edge. This completely reconstructs the
path-reversal deletion defect from two port pairs and the d-bit lift data.

**Proof.** The pre seam receives the vertices indexed by E on the left
and by the complement of H on the right. Thus it is empty exactly in the
first stated case. The post seam receives the complementary left subset
and H on the right, giving the second case. Since d>=1, these two exceptional
bit patterns are different. No other old adjacency can be cut by the new
vertices. Apply the exact defect identity(7) from the main note. QED.

Thus the vertical insertion group F2^(2d) grows with d, but its effect on
the deleted move target passes through a constant-size incidence sidecar.
This is an exact theorem for this family, not a uniform bounded-state model
for arbitrary tournaments.

## Exact all-depth defect counts

Let q=2^(k-1)+1 and D=2^d. The source has Q^2 paths, where Q=D(q-1)+1.
The numbers with correction size0,1,2 are respectively

    N0 = D^2 + 2(q-2),
    N1 = 2(q-2)D^2 + 2D + 2(q-2)^2,
    N2 = (q-2)^2(D^2-2) + 2(q-2)(D-1) + 1.             (A1)

These sum to Q^2. At d=1 they reduce to the one-step formulas in the main
note.

**Proof.** The pre seam has two actual vertex endpoints exactly when
I is not the empty finite subset; the post seam has two actual vertex
endpoints exactly when K is not the empty finite subset. There are q-2
nonempty finite states and one infinity state in each coordinate.

- Empty/empty: all D^2 lifts have defect0.
- Empty/nonempty-finite, or its reverse: each such parent has one defect0
  lift and D^2-1 defect1 lifts. There are2(q-2) such parents.
- Empty/infinity, or its reverse: each has D defect1 lifts.
- Two nonempty finite states: each of the (q-2)^2 parents has two defect1
  lifts and D^2-2 defect2 lifts.
- Nonempty-finite/infinity, or its reverse: each of the2(q-2) parents has
  one defect1 lift and D-1 defect2 lifts.
- Infinity/infinity: its unique lift has defect2.

The forced full bits at infinity are essential in the third and fifth
rows. Summing the disjoint rows gives(A1). QED.

## Check

    python3 check_all_depth_ports.py > all_depth.out
    python3 -O check_all_depth_ports.py > all_depth.optimized.out
    cmp all_depth.out all_depth.optimized.out

This checker compares literal deleted-target tournaments with the two-seam
formula for every path state at k=2,...,4 and d=1,...,3, then checks every
stratum count against(A1). It imports only the main packet's local elementary
tournament/path constructors, not repository code. Its coverage is reported
in `all_depth_ports.json`. All-order status rests on the proof, not on this
finite control. This addendum has no separate literature-priority claim.
