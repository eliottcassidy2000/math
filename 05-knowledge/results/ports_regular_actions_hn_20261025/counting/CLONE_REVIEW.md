# Independent review of the clone-cut cancellation

2026-10-25. **PASS**, following separate proof inspections by the root and
the port-review agent. The producer supplied the proof and literal checker;
the independent review concentrated on the inflation fibre, clone mergers,
and equivariance. No complexity or black-path transport claim is inferred.

The exact statement is [CLONE_INVOLUTION.md](CLONE_INVOLUTION.md), equations
(I1)--(I6). The coefficient meaning was independently repaired earlier in
[INDEPENDENT_REVIEW.md](INDEPENDENT_REVIEW.md). The canonical path decomposition
is the proved [THM-4463, section2](../../../../01-canon/theorems/THM-4463-homogeneous-tournament-contacts-and-pinned-scale-floor.md),
whose theorem status and proof were checked directly.

1. An inflated vertex block is a transitive module. It lies in one component
   of the homogeneous-pair path graph. It is contiguous there: a skipped
   interior vertex would distinguish its earlier and later members. Conversely,
   every such interval is a module. Thus every possible source quotient is
   a cut partition of these canonical paths, even if its original vertices
   cease to be maximal clone blocks after inflation.
2. An isomorphism of weighted quotients preserves block sizes and all cross
   arcs. The unique order isomorphisms inside the transitive blocks therefore
   lift it to an automorphism of the fixed inflated tournament. This proves
   the full iff between weighted-isomorphism classes and cut-set orbits;
   merely having the same block-size multiset would not suffice.
3. The chosen orbit must depend on the fixed unweighted tournament alone.
   The least canonical code of a distinguished nontrivial module selects
   precisely an automorphism orbit and meets that condition. This selection
   is finite, but no efficient canonical-labelling algorithm is claimed.
4. Tournament automorphisms have no involution; the elementary inverse-pair
   argument makes the whole finite group odd. Its module orbit consequently
   has odd size. Intrinsic orders make the set of first gaps invariant, so
   symmetric difference with that set commutes with every automorphism.
5. The number of cuts changes by an odd integer. Thus the sign flips and
   no orbit can be fixed, since an automorphism preserves the number of cuts.
   The inflated vertices, arcs and total N are unchanged. Pair-free fibres
   have no gaps and exactly one positive source object, including N=0.

The two hostile controls exercise distinct implications. In C3[2,2,2],
one labelled gap does not descend to cut orbits; all three first gaps do.
In C3[2,2,1], toggling both nontrivial modules preserves sign, so one must
select one odd automorphism orbit rather than their union. These expose the
actual role of odd-order symmetry instead of treating it as a count coincidence.

The root replay invokes the producer's no-import literal checker in normal
and optimized modes. It independently compares weighted quotient canonical
codes with cut orbits for all tournament classes through order5 and both
hostiles. This is computational corroboration; the five steps above prove
the unbounded identity. It gives a direct combinatorial proof of clone
removal, while the computation of unrestricted tournament numbers still
uses the separate degree-response recurrence.
