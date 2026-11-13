# A finite retained-star handoff at29 and its failed uniform enlargement

**FINITE-EXACT physical repair**, with an elementary lower bound on its
arc-edit radius. This is a marked handoff certificate, not a new occurrence
of counts29 or89 and not an all-target ready-return theorem.

The repaired [H9 core](history/ONE_ARC_H9.md), rooted at0, is
`(22,8,2,21,6)`. Apply N R. The graph at the exact target29 is

    T=(54,8,2,21,6,94,31),  root5,
    (H,s,q,F)=(29,9,9,29).

Its only incoming root arc is0→5, but q−s=0, so it is not unit-ready.
Reverse the three old arcs on the unordered pairs02,13,46, obtaining

    U=(50,0,3,23,70,94,15), root5,
    (H,s,q,F)=(29,8,9,30).

All changes avoid root5. Its entire root star and the physical predecessor0
are retained. The new defect q−s is1, so this is a unit-ready rebase at
the same exact count29. Vertex1 has become a global sink; deleting it
preserves all four root observations and leaves a strong six-vertex
tournament. The physical root/predecessor are not deleted. Applying N R
therefore gives the next exact target89. Its returned root again has
defect0; no ready return at89 is claimed.

The three-arc count is minimal among edits on this fixed seven-vertex
source that preserve its root star: there are15 eligible nonroot pairs,
and every one- or two-pair reversal (15+105=120 candidates) was counted
exactly. None has both H29 and readiness defect1. This is a finite local
minimum; changing order or using other source graphs is outside it.

The six-vertex ready source was recovered from the inherited small bank
`actual_hamilton_multiplicities_20261021/session24-multiplicities/extended_census.json`
at `c022b73c3cb5`. Aligning its root and unique predecessor with5 and0,
then adjoining a neutral sink, supplies the displayed witness. The script
recounts it directly instead of treating the bank's source status as a
proof of the new marked incidence.

A tempting extension fails at its first new size. Inflate vertex2 of the
H9 core into TT_b, keeping all other vertices singletons. After N R,
reverse every arc from vertex0 to that block, reverse1→3, and reverse the
arc from the last newborn to vertex4. The named root star still survives
and vertex1 is again a sink. At b=1 this is the repair above. At b=2 the
original next target is53, whereas the modified source has count59 and
is unit-ready. Thus the repair is **not** an exact-target family under
this inflation. Controls through b=6 retain the same distinction.

The later [all-module law](constructors/H9_INFLATION.md) resolves this
failure exactly. For any nonempty substituted module M, the modified
target exceeds the original by 3 C_2(M), where C_2 counts ordered
nonempty two-path covers. Every nonsingleton tournament has C_2(M)>0
by cutting a Hamilton path, so this common inflation preserves the
exact target only for a singleton. Its new readiness defect is H(M);
unit readiness survives precisely for transitive modules.

This failure preserves a useful part of the construction: readiness can
survive an enlargement even while its absolute target drifts. Both target
and readiness must be tested together. The all-parameter exact-prefix
repair and retained-chain renewal elsewhere in this packet satisfy their
own proved target/interface laws; this finite handoff does not replace
those laws or bypass their remaining target-supply boundary.

Run `python3 -B 05-knowledge/results/h17_exact_repairs_20261112/ready_handoff_probe.py`.
