# A source switch at 105, and why a neutral birth cannot do it

**PROVED:** the inherited neutrality theorem freezes the ordered nontrivial
SCC word under count-neutral vertex additions and deletions. **VERIFIED
explicit construction:** two vertex replacements switch a 317-ready H105
source into a 319-ready H105 source, with intermediate counts19 and9.
**FINITE-EXACT:** bottleneck9 in the complete fixed-order replacement graphs
at orders8 and9. **OPEN:** an unlimited reset compiler with joint ports.

## R1. The precise obstruction to a neutral reset

By [THM-1900(B), insertion-response calculus](../../../../01-canon/theorems/THM-1900-the-insertion-response-calculus.md),
adding one vertex preserves H if and only if the new vertex is a singleton
SCC. Equivalently H(T-v)=H(T) if and only if v is a singleton SCC of T.
This is an inherited all-order theorem with a repaired proof, not a newly
discovered neutrality criterion.

It follows that the ordered isomorphism-class word of nonsingleton SCCs
is invariant under **every** sequence of count-neutral births, count-neutral
deletions, relabellings and changes of the active root. Internal singleton
gap lengths can change; the nonsingleton components and their order cannot.
Conversely two unmarked tournaments with the same ordered nonsingleton word
are connected by such neutral moves: remove their singleton components to
the common word and insert the desired singleton gaps. This converse has no
promise of retaining ports that were on the removed singletons.

The [all-order target classification P3](../targets/PROOF.md) says every
105->317 source and every105->319 source has exactly one nonsingleton SCC,
and their possible strong cores are disjoint. Therefore the neutral moves
above cannot switch between these two branch families. Root switching or
interior *insertion* cannot evade the conclusion: the target audit already
permits every cut of the new vertex. A macro that changes old arcs, replaces
old vertices, or temporarily changes H is a genuinely different operation.

This is a direct new application of the inherited theorem to the105 hostile.
It does not exclude all possible count-preserving macros or arbitrary graph
rewrites. In particular the following count-dipping macro succeeds.

## R2. An actual two-replacement reset

Rows are bitmasks; bit j of row i means i->j. Set

    A = (254,204,168,16,38,74,156,56),
    B = (46,60,24,16,33,76,159,63),
    C = (158,140,24,16,34,79,159,60).

All three are genuine tournaments and have H105. A is a source singleton
followed by the seven-vertex strong core with Hamilton-order code1346931.
B and C are strong. A has a one-birth child of count317 at cut186; C has
children of count319 at cuts26 and31. A cut is the set beaten by the new
vertex, encoded in the displayed source labels.

Delete vertex5 from A. The surviving graph, in increasing surviving-label
order, is

    (126,108,72,16,6,92,24), H=19.

Adjoin a new vertex with cut123. The resulting graph is B after the exact
surviving-name map

    0->7, 1->0, 2->5, 3->3, 4->4, 6->1, 7->2;
    the newly added vertex corresponds to B's vertex6.

Now delete vertex4 from B. The surviving graph is

    (30,28,8,0,44,79,31), H=9.

Adjoin a new vertex with cut34. The resulting graph is C after the map

    0->0, 1->1, 2->2, 3->3, 5->7, 6->5, 7->6;
    the newly added vertex corresponds to C's vertex4.

Thus the physical reset has count trace

    105 -> 19 -> 105 -> 9 -> 105.

All arcs between surviving vertices are retained at each operation. The
maps compose to preserve six original named vertices:

    original 0,1,2,3,6,7 -> final 6,0,7,3,1,2.

Any of the15 pairs among these six names can be retained as named ports.
This does **not** preserve every initially prescribed pair: original
vertices4 and5 are removed. Nor does it preserve the numerical joint
observer at the retained ports. The checker carries the original pair(1,2)
through both maps, recounts its full two-root deletion cube and stores its
residues modulo49 and343. The initial and final cubes differ. Complete
rows and maps permit further endpoint/attachment coordinates to be computed;
the four-count deletion cube is not advertised as the 17-coordinate observer.

The construction is a finite reusable template at count105. Reversing the
steps gives a reset back to A. It neither proves source-switch connectivity
at arbitrary h nor supplies a monotone one-birth induction. It demonstrates
that the no-neutral-reset obstruction identifies a repairable operation
boundary, rather than a new forbidden count.

## R3. A small exhaustive test of reset depth

The inherited audited strong inventory has four H105 classes at order7,
eleven at order8, and three at order9. The only possible product of two
nontrivial realizable counts equal to105 is3*35: factors7 and21 are
forbidden. An H3 strong component requires3 vertices; the two H35 strong
components require7. Hence no H105 tournament of order8 or9 has two
nonsingleton components. Every such tournament is an exterior transitive
padding of one of those18 strong cores. There are exactly19 classes at
order8 and37 at order9.

For each fixed order, connect two classes if deleting one vertex from each
gives isomorphic tournaments. The edge records that residual H. Traversing
it means delete once and insert once, with exact surviving-name transport.
The checker constructs every such edge, then maximizes the minimum residual
H along a route from a317-ready class to a319-ready class. At both orders
the optimum is9. The displayed route attains it. Thus within either of
these two finite fixed-order universes one must at some point dip to at
most9. This is **not an all-order bottleneck theorem**, and it does not
cover larger replacement blocks or old-edge reversals.

`check.py` uses recursive subset counting; all displayed sources, cards
and branch children also have independent literal-permutation counts.
Canonicalization minimizes over every actual Hamilton order, an
isomorphism-invariant nonempty set, without an orbit heuristic. Complete
bank status is inherited, not regenerated here. The bank hash, target
receipt hash, graph universes and surviving maps are in `check.json`.

## R4. The induction obligation after this test

A constructive controller now needs three distinct mechanisms:

1. A target-aware choice of an actual source, with named ports and a cut,
   for the next3h+2 or3h+4 branch.
2. A proved local reset template when the available source is in the wrong
   branch family, including any permitted count excursion and port deletion.
3. An exact update of the full retained attachment state after that reset,
   including lifted phases. Equality of H before and after is insufficient.

The H5 compiler identifies the analogous incoming difficulty at173: the
parent-count bank1,3,5 is insufficient, and parent9 is necessary and
sufficient there. The golden contraction theorem supplies actual mixed
words with controlled7-adic valuation, but it changes H and retains the
normalized phase. These are complementary pieces of a controller; neither
is already a source-reset theorem in the unrestricted two-port language.

A global reset rule must be given by actual templates and a decreasing
inverse or an invariant guaranteeing termination. Allowing a reset to
delete everything and rebuild an unproved target would merely hide the
original realization problem inside an oracle. The explicit surviving
six-vertex interface here avoids that issue at105, without asserting a
universal analogue.
