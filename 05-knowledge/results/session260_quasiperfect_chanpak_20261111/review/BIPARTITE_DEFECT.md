# Focused review: the all-size bipartite defect theorem

**Verdict: PASS. No material correction required.** Reviewed only the new
all-size proof and exact equality criterion in
`../bridges/BIPARTITE_DEFECT.md`; inherited classical results were accepted
at their explicitly named scopes. No large census was repeated.

## Endpoint existence and the injection

The cut/set map is a bijection. Its two domination conditions say precisely
that each unselected A vertex has a selected B out-neighbor on the left and
each unselected B vertex has a selected A in-neighbor on the right.

For a nonempty left side, a terminal SCC wholly inside A would have a forced
outgoing edge to a B vertex outside that SCC, contradiction. Thus its terminal
SCC meets B. A singleton component is immediate; otherwise Camion's cycle
can be cut to end at a chosen B vertex. Concatenating paths through the earlier
SCCs gives the required spanning endpoint path. The initial/right argument
is exactly dual. This uses endpoint freedom only inside a strong component,
not the false claim that an arbitrary tournament has a path starting or ending
at every vertex. Empty sides are correctly assigned one empty path.

Canonical selection of those two paths gives a literal injection. Splitting
at z recovers the cut and hence the dominating set, so different sets cannot
collide. The displayed defect formula accounts for all Hamilton paths and
all dominating sets; each summand is nonnegative for the stated reason.

## Exact equality locus

The all-selected dominating fibre is H(T[B])H(T[A]). Under equality it is one,
forcing both shores transitive. This includes empty shores.

If the root-deleted tournament remains cyclic, a directed triangle must now
be mixed. For two A vertices and one b in B, deleting only b from the selected
set gives a dominating cut with right side A union {b}. Its triangle SCC is
the unique possible nontrivial SCC, since every cycle must use b. If that SCC
is first, cutting a Hamilton cycle at two distinct A vertices gives two
A-starting paths. If it is not first, every earlier SCC is an A singleton,
so any of the at least three cycle cuts through the nontrivial SCC gives a
globally A-starting path. The opposite shore is transitive and contributes
exactly one left path. Thus this single dominating fibre already has positive
excess. The one-A/two-B case is genuinely dual with ending vertices.

These cases exhaust nontransitivity. Conversely a transitive root-deleted
host is exactly the earlier chain-word adapter. The iff is therefore proved,
not inferred from finite equality counts. Parity then upgrades a positive
gap to at least two.

For a fixed bipartite input, allowing arbitrary shore orientations does not
weaken the obstruction: zero defect forces a transitive full host, whose total
order realizes precisely nested cross-neighborhoods. Conversely a chain graph
has such an order. Hence the fixed-input zero-defect criterion is exactly
“chain graph,” with the declared apex/cross-incidence interface retained.

## Tiny independent control

`bipartite_defect_tiny.py` imports no author code. It uses literal permutation
paths and direct dominating-set enumeration on all64 labelled four-vertex
tournaments, all four roots, and all eight cuts per root:256 root cases and
2048 cut-fibre cases. It also checks all four shore orientations of a fixed
2K2 input. All pass; the latter have d=9 and H-values13,15,15,13.

See `bipartite_defect_tiny_receipt.json`. This is only a small hostile/control
check; the all-order conclusions rest on the SCC and fibre proof above.

No universal zero-defect transfer from the solved graph-domination spectrum,
no all-odd Hamilton-count coverage, and no external novelty claim is inferred
by this review.
