# An actual alternative H141 source with the425 branch

**PROVED by explicit construction and exact counting.** This is a
same-count change of realizer followed by two births. It is not a
count-preserving native path from either old source, and it is not a
one-birth certificate or a general source-return theorem.

The alternative lane's all-order frontier excludes425 from every induced
supertournament of each of its two specified neutral C4 sources of H141.
That obstruction concerns those actual sources. It does not survive
forgetting the source and retaining only its count.

The current bank's **singleton** chart, O=0,A=1,B=61, gives the strong
9-vertex tournament

    T=(254,508,504,240,480,448,384,256,9), H(T)=141.

Adjoin a new vertex beating only old vertex7, then one beating only the
preceding newborn9. In outgoing-mask convention the two cuts are128 and512.
The exact intermediate and output are

    (766,1020,1016,752,992,960,896,256,521,128), H=215,
    (1790,2044,2040,1776,2016,1984,1920,1280,1545,128,512), H=425.

Every old arc is preserved through these two births. All three graphs
are strong. Independent subset-DP counts verify141,215,425. On the
converse source, vertex7 has the exact ready profile

    (H,s,q,F,F-H,delta)=(141,68,69,142,1,0).

Thus the two sole-exit births are literally the converse of the standard
1-ready reset, giving3*141+2. The positive certificate and the negative
frontier jointly demonstrate why an actual source identity is essential.

The exploratory current-bank inverse found42 raw normalized witnesses
at H141 (orders8:14,9:20,10:8), including19 distinct strong classes.
None of those19 strong representatives had a one-birth425 cut; the
displayed two-birth hit was found after237 distinct intermediate strong
classes. These search statistics do not prove all-order absence for
other H141 sources or for omitted weak representatives. In particular,
the inherited complete small-base J3 certificate mentioned by the parent
session already supplies one-birth141->425 existence in a larger census;
this note does not replace that fact with a negative claim.

## A genuine infinite family, and its limited intersection with C4

The same singleton chart with A=1 and B=2^m-3, m>=4, has

    H=9*2^(m-2)-3.

Indeed g(2^m-3)=2^(m-2)+1 for m>=3, and the chart gives
H=g(B)+2B+2. Its last two transitive-base vertices both beat the apex,
so the converse at the last base vertex is1-ready. Appending the same
two sole-exit births yields3H+2 for every m>=4. At m3 the last-two
orientation fails and H15 is not ready: the boundary cannot be dropped.

Here is the structural readiness proof, which avoids extrapolating from
the counts. In any one-apex tournament whose base is transitive, assume
the first base vertex is beaten by the apex and the last two base vertices
beat it. The whole tournament is strong. Its last base vertex r has the
apex as its only out-neighbor. Delete r. Every Hamilton path of the
deletion except those ending at the apex extends to a path ending at r.
Exactly one deletion path ends at the apex, since deleting the apex leaves
the unique transitive path and its last vertex points to the apex.
Thus on the converse q-s=1 and the root has indegree1, giving delta0
and F-H1 by the inherited rooted defect identity.

This specific dyadic source family equals the old neutral C4 value
34*2^k+5 only when k2,m6. The equation is
`9*2^(m-2)=34*2^k+8`. For k>=3 the right side has valuation3 at2,
forcing m-2=3 and an incompatible size. The k1 case fails; k2 gives144.
Therefore the explicit141 rebase is not an all-k repair of the C4 family.

A larger scoped family follows from the same structural argument. If the
apex cut word on a transitive base is the ordinary binary word for4u,
u>=1, its final two bits are00, and its converse has a1-ready root.
Its count is

    g(4u)=g(u)+3u+(u mod2).

For a target H, u< H/3 bounds the complete inverse of this subfamily.
Among the C4 targets34*2^k+5 for k1..8, it supplies parameters

    k1: u=21,23; k2: u=36; k7: u=1140;
    k3,4,5,6,8: no parameter in this subfamily.

These latter rejections do not say the targets have no ready presentation
in the complete old ready grammar. They identify the limit of this
particularly simple source rebase.

## Reproduction

```text
python -B 05-knowledge/results/chanpak_five_core_bridges_20261110/history/replay141.py
```

[alternate141.json](alternate141.json) records the literal source,
extension cuts, full ready profile, finite family controls and the complete
scoped k1..8 inverse. The exploration receipt is
[alternate141_search.json](alternate141_search.json); it is provenance for
the discovery, not a no-other-source theorem.
