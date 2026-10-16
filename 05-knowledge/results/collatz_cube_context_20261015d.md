# An exact terminal-cube obstruction, repaired by one native edge

2026-10-15. **PROVED:** the signed-cycle packing inequality below.
**FINITE-EXACT:** an exact free-terminal cube optimum and a one-edge repair
for one declared annular controller. **OPEN:** a uniform overlap rate and
universal ROOT entry. This is a companion to
[odd-ray overlap](collatz_odd_ray_overlap_20261015d.md), not an all-scale
certificate. Self-audited; no independent referee is claimed.

[Program](../../04-computation/experiments/collatz_cube_context_20261015d.py),
[complete certificate](collatz_cube_context_20261015d.json).

## 1. The old failure was not solely the entrywise relaxation

The [annular packet](collatz_annular_overlap_20261015c.md) has a deliberately
retained failed case: X=4096, bin length b=64, every original source
1<=n<=32768, initially 16 shortcut steps per source, followed by 16 steps
from every compiled boundary terminal. Both stages stop at ROOT and retain
the union of observed edges. There are 298 unresolved original sources.
After zero columns are removed the annular Gram has 51 active terminals,
including ROOT. The row multiplicities, source intervals and graph are unchanged
throughout this comparison.

Write G_int=sum_I B_I B_I^T, so the actual coarse quantity is
z^T G_int z/64, with z_ROOT=1 and every free terminal in [-1,1]. Previously

    ||G_int||_entry1 /64 = 9767515472/64
                        = 9.0967076569... X^2,

which misses the chosen 9X^2 budget. A sharper terminal-label calculation is
possible, but it still misses:

    max_(z_ROOT=1, |z_e|<=1) z^T G_int z /64
        = 9716780760/64 = 9.0494572743... X^2 > 9X^2.          (1)

This is an **exact optimum**, certified by matching upper and lower witnesses;
no numerical optimizer or unproved terminal label is used in verification.
It closes a precise unproductive route: every upper bound valid for all labels
in this fixed terminal cube must exceed the budget. This is not a lower bound
on the actual Collatz energy, because the maximizing labels are artificial.
It is not an impossibility theorem for richer native constraints.

## 2. Keep sign compatibility as a context

For a symmetric positive-semidefinite Gram matrix G, define a signed graph on
its terminal columns. Nonzero G_ef gives the desired sign z_e z_f=sign(G_ef),
with capacity |G_ef|. Zero entries are absent edges, not arbitrarily oriented
ties. A cycle is frustrated if the product of its edge signs is -1.

Give each selected frustrated cycle C a nonnegative weight p_C, subject to

    sum_(C containing edge e f) p_C <= |G_ef|.                (2)

Then

    max_(z in [-1,1]^m) z^T G z
          <= ||G||_entry1 - 4 sum_C p_C.                     (3)

Proof: a convex quadratic attains its maximum over a cube at a vertex
(successively choose an endpoint in each coordinate). For a sign vertex, every
frustrated cycle has at least one unsatisfied edge. Its weighted capacities
imply that the total absolute mass of unsatisfied edges is at least sum_C p_C.
An unsatisfied unordered edge loses 4|G_ef| relative to its contribution to
the entrywise norm. This proves (3). The bound remains valid when the ROOT
coordinate is fixed to +1. Cycles may share vertices and edges; only the
capacity constraints (2) matter. This is a finite signed-graph inequality,
with a proof here and no novelty claim.

For (1), the program finds 1092 frustrated triangles and exports an integer
packing of 52 of them. Every triangle sign and every edge capacity is verified
again from the exported packing. Its total weight is 12683678, giving the
upper bound 9767515472-4(12683678)=9716780760. A literal 51-entry sign vector
with ROOT=1 attains the same value. The lower witness is evaluated twice:
by the Gram quadratic and independently by summing the squares of the actual
weighted source rows. This pair of certificates proves optimality without
searching all 2^50 sign assignments.

The first odd-ray cone has a fixed sign per terminal across physical bins,
so its Gram graph is balanced and this improvement is unavailable there.
Later sign changes can create frustration. Thus the failure of that cone
points to another useful invariant rather than merely invalidating a proof:
retain compatibility of terminal signs, in addition to their pairwise absolute
masses. Frustration and native overlap pay different invoices.

## 3. A single new edge pays the remaining budget

The fixed graph has no outgoing edge at terminal 70924. Add just

    70924 -> 35462.                                          (4)

Every edge of the 112-step shortcut path from 35462 to ROOT is already in the
old observed-edge union. The certificate lists that entire suffix. Therefore
(4) compiles the 113-step receipt

    70924 -> ROOT,    t(70924)=-1.                           (5)

The maximizing artificial assignment in (1) gives this terminal +1. The new
native relation excludes that witness, and direct signed transport gives

    ||G_repaired||_entry1 /64 = 9209359132/64
                             = 8.5768840573... X^2 < 9X^2.  (6)

Rebuilding the rows from their original sources gives exactly the same Gram
as transporting the old one. The exact raw weighted payment is 558156340;
there is no zero-label deletion. This is cancellation against the marked
ROOT column with its correct clock parity. The repaired graph still has
285 unresolved original sources, down from 298, and maximum compiled native
clock 180. It is the weighted invoice, not the number 13 of newly grounded
sources, that proves (6).

Within this fixed graph plus free-terminal-cube interface, no additional
relation can be omitted altogether: (1) proves that zero new information
cannot certify 9X^2. One added native edge suffices. This is not a claim of a
globally optimal observation strategy. Discovery inspected up to 16 native
steps at each of 50 active nonROOT terminals; **one edge is the additional
exported certificate, not the total discovery cost**. The 112 previously
checked suffix edges and the original graph are essential dependencies.

## 4. Transfer to the continuing research target

The useful move is concrete: first certify that a bound is limited by the
actual free-label model, then use an extremizing context to identify a native
relation that removes it. Here the relation was a one-edge bridge into an
already-paid suffix. This strengthens the older complementarity example:
a cheap added edge can carry a long integer-clock certificate only when the
entire observed graph is retained.

A candidate all-scale structure should therefore carry

    original weighted rows + observed native graph + Gram sign graph
       + extremal label contexts + source/bin profiles.

Its map to the old Gram preserves the actual observable and all source
multiplicities; discarding the sign graph loses simultaneous compatibility,
and discarding the native graph loses the 112-step reusable suffix. Neither
an extremal label vector nor a frustrated cycle proves a particular integer
reaches ROOT. The next arithmetic question is whether every excessive actual
annular context admits a *quantitatively priced native refutation*, uniformly
in X. The present example supplies an exact witness and repair, not that
uniform theorem.

This can be combined with the cone-with-defect invoice (16) in the odd-ray
note: use native overlap where a contracting packet is available, cycle
packing where physical-clock sign changes create frustration, and an explicit
native extension when the exact cube optimum remains too high. All three
steps preserve the original source and integer clock.

## Verification

```sh
OPENBLAS_NUM_THREADS=1 python3 04-computation/experiments/collatz_cube_context_20261015d.py
OPENBLAS_NUM_THREADS=1 python3 -O 04-computation/experiments/collatz_cube_context_20261015d.py
```

The output contains the full integer Gram, terminal names, packing, extremizing
sign vector, original graph parameters, new edge, old suffix, transport invoice,
and source counts. Tiny balanced/frustrated Gram controls verify both the
positive and hostile boundaries of (3); no sampled label is treated as true.

Normal and optimized Python outputs match byte for byte. SHA256:
`7fe6c5d0df795c23ec3bd363564e7eeb8dad02cfe262d41553d66da4cdb25b28`.
