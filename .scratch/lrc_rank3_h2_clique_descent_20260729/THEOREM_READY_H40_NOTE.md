# Ranked-r3 `H2<=40` ordered-Hunter closure

**Status: FINITE-EXACT, VERIFIED, SCOPED.**  This closes 62 marked suffix
branches in the locked `K>=20` scalar-hard universe.  It is not a root
closure, does not close every strongest-r3 branch, and does not prove
LRC(14).

## Exact statement

Use the hash-pinned all-root scalar ledger and restrict to rows satisfying

```text
K>=20,
the scalar five-singleton test is hard,
the strongest strict ranked-complement margin has rank r=3,
|H2|<=40.
```

There are exactly 62 such rows.  On every row, let `C` be the literal
post-apex carrier, `h=mu(C)`, and

```text
R3=q1+q2+q3 < 5h/7,
theta=h-R3,
H2={v:c_C(v)>=theta/2}.
```

No five allowed labels cover `C`.

## Structural reduction

THM-2893 and THM-2897 force at least four labels of a hypothetical
five-cover into `H2`.  Every pair among those four labels has

```text
U_C({x,y})>=theta.
```

Thus the four labels form a `K4` in the symmetric theta-heavy graph on
`H2`; after subtracting the corresponding four danger combs, only a
singleton cover obligation remains.

The core is global, not a truncated sample.  If `r(C)` is the component
count and

```text
epsilon=theta/2-h/7=(5h/7-R3)/2>0,
gamma=(99/70)r(C)/7,
```

then the discrepancy estimate seals the tail at

```text
N=ceil(gamma/epsilon)-1.
```

All allowed speeds through `N` are evaluated exactly and equality at the
`theta/2` core threshold is retained.

## Ordered proof partition

The heavy graph and its pair-union weights are undirected and symmetric.
The deterministic elimination order is only an acyclic proof-branch
gauge.  At a current pivot `x`, every `K4` assigned to `x` has its other
three vertices among the remaining heavy neighbors of `x`.

Write

```text
d_x(y)=c_(C minus D_x)(y)=U_C({x,y})-c_C(x).
```

If `d1>=d2>=d3` are the three largest such child coverages, then

```text
c_C(x)+q1+d1+d2+d3<h
```

closes every assigned `K4` at once.  Fewer than three remaining heavy
neighbors is a vacuous closure.  This is the relative Hunter star
identity from THM-2897, with the parent `q1` retained as a valid cap for
the final label.

If the star invoice does not close a pivot, each forward-neighborhood
triangle is enumerated once.  For its `K4`, put

```text
psi4(Q)=sum_(v in Q)c_C(v)
        - max_(spanning trees T on Q) sum_(xy in T)i_C(x,y).
```

The leaf closes whenever

```text
q1+psi4(Q)<h.
```

The maximum tree credit is computed both by brute enumeration of all 16
labelled spanning trees of `K4` and by deterministic Kruskal.  Remaining
tree-hard leaves are subtracted literally, with simultaneous and
sequential subtraction required to agree, and are closed by an exact
longest-component singleton horizon.  Every strict inequality is used in
its proved direction; equality is retained as hard.

## Locked census

```text
branches                                      62
closed / surviving                          62 / 0
H2 vertices                                  1662
heavy edges                                 12090
heavy-threshold equality edges                  0
initial K4s                                 93754

star steps / star-assigned K4s          1454 / 81429
paid steps / paid K4s                    208 / 12325
tree-closed / tree-hard                12286 / 39
tree-margin equalities                          0
literal singleton-closed leaves                 39
empty residuals / witnesses                  0 / 0
singleton scans / independent controls    295 / 67
star-only branches                              29
```

The ordered partition assigns all `93754` original `K4`s exactly once.
As an independent, less compressed positive control, literal enumeration
on all 42 rows with `|H2|<=30` closes all `22631` `K4`s after `213410`
singleton scans and finds no empty residual or witness.

## Locked artifacts

```text
.scratch/lrc_rank3_h2_clique_descent_20260729/ordered_hunter_pivot_H40.py
SHA-256 ee092030eaa5f69e064148c42418851581d717a80752bc73346e22d637acdb27

.scratch/lrc_rank3_h2_clique_descent_20260729/ordered_H40_summary.locked.out
SHA-256 7c8fbbb82282d0f2f2e5c687b5d56a90c1a684b9d5a7487e491f66410f023043

.scratch/lrc_rank3_h2_clique_descent_20260729/ordered_H40_ledger.locked.out
SHA-256 9320782e37c01028cbcdfd6d1076db167b766628a9fde70b23ba0d61b0f377b1

canonical semantic ledger
SHA-256 2088b2016179d14e1b5e76d6c9eddfe88e8352c943f49970bfb201b46583cd47
```

Ordinary and `python3 -O` summary outputs are byte-identical.  Ordinary
and optimized full-ledger outputs are also byte-identical.

Pinned inputs:

```text
pilot.py
SHA-256 9ec1b8d8c697f54fdbf2836638fa6c7dc5284f4ea2644849d98d71398c1f3520

core.discovery.out
SHA-256 5569d8d34b59e5eed2bfa82148648dcb5e63515146ff55b95234a002d0c87ba2

graph_H40.discovery.out
SHA-256 5683dc23f70d78a6c351008fd262a6d3214905650733446feac12155ccb58ad7
```

## Canon route recommendation

Do **not** place this in THM-2902.  The current remote THM-2902 is already
a promoted, sharply scoped omission-six ranked-H1 root closure, whereas
this result concerns 62 ranked-r3 marked suffixes spread across the
`K>=20` universe and yields no root closure by itself.

The later remote scan at
`origin/main=08bbf52695` found THM-2903 proved and THM-2904/2905 reserved;
no new ID should be reserved.  THM-2905 is the closest live namespace,
but its stated intended scope is the global piecewise-linear `G5`
envelope and whole-root consequences.  This finite-core ordered-star/tree
certificate is therefore best kept as a clearly scoped sibling control,
or folded into THM-2905 only if that theorem is deliberately expanded to
separate the global envelope from the literal finite-core refinement.

A stronger route is to postpone promotion until the graph-free
continuation closes a materially larger portion of the 143-row
strongest-r3 universe, then canonize the common ordered-Hunter algorithm
and its complete scoped census together under an ID allocated by the root
session.
