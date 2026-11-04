# Two actual nonzero rearmings of the prepared H=1429 output

**PROVED — [independent analytic and exact replay audit](AUDIT.md) passed.** The new positive statement is about one actual prepared
tournament, not a reusable two-branch grammar. The all-clock statements below
are exact response formulas with an explicit, unique neutrality clock.

## 1. Inherited interface and the missing operation

The input is the literal clock-three output of the audited
[PREPARED_FOUR](../../multiport_grammar_20261103/component/PREPARED_FOUR.md).
It has eleven physical names `0,...,10`, physical root **8**, and rows

```
(1790,1668,1672,1778,1670,1942,1718,1280,1631,128,512).
```

Its actual observers are `H=1429`, `F_8=1563`, and root defect134. Rows
encode outgoing adjacency by bit masks. The newest auxiliary9 has zero
start gap; the inherited C2 converse excludes gap1 or2 at **every** vertex.
Thus neither of the preparations here is re-marking, or another induced C2
step. Both change old arcs and retain every old physical name.

The inherited exact joint interface is
[construction/RESULT](../../multiport_defect_grammar_20261103/construction/RESULT.md).
Write `P` for the path polynomial including the empty path, `W_r` for
nonempty paths ending at the retained root, `A_a` for paths starting at an
auxiliary, and `Q_ar` for paths from a to r. Full square-zero coefficients are

```
H=[P], F=[W_r P], S=[P A_a],
L=[W_r P A_a], K=[P^2 Q_ar], G=L+K.
```

The two permanent C2 births give `(H,F) -> (H+2S,F+2G)`. All physical
root readings in this note refer to old label8, including when a temporary
helper is the final array entry.

The closest mechanism is the C2 joint response, the canonical hostile is
its every-auxiliary nonreadiness theorem, and the failed near miss is equality
of H without a joint root response. The least-used sidecar that succeeds
here is the complete one/two-ready funnel classification in
[READY §§1–3](../../two_entry_funnels_20261103/READY.md). The live comparison
coordinates are H, physical F, start gap, joint deficit, actual strongness,
and whether a declared class accepts the output again.

## 2. A useful exact joint identity

Suppose a has r as its **sole out-neighbour**. Put B=T−a and C=T−{a,r}.
Then

```
L=F_r(B),       K=P_C(2),       F_r(T)=F_r(B)+P_C(2),
G=F_r(T).                                                   (1)
```

For L, the A_a path must avoid r, since W_r already uses r, so that path is
the singleton a. For K, Q_ar is exactly the two-vertex path ar. To split F,
if a is in the path ending at r, that path ends ar and its earlier part is
an arbitrary path of C; the other path supplies the second ordered bin.
If a is in the other path, that path ends at a and deleting a leaves the
ordinary F_r(B) cover. This proves all identities, including empty-bin
boundaries, without an assumption that fixed or eligible paths are retained.

Consequently a target with S−H=sigma has exact C2 response

```
(H,F,b) -> (3H+2sigma, 3F, 3b−2sigma),  b=F−H.               (2)
```

Both prepared targets below satisfy the sole-out-neighbour hypothesis.

## 3. The two preparations on the actual old names

For the gap2 preparation reverse, in the displayed order, the six old pairs

```
(1,8), (1,9), (2,3), (4,5), (7,10), (8,9).
```

The resulting literal target is

```
(1790,1412,1664,1782,1702,1926,1718,256,1117,386,640).
```

The auxiliary is old vertex7, and its only out-neighbour is old root8.

For gap1, first permanently add name11 with sole out-neighbour old root8;
every other old vertex beats11. Then reverse

```
(2,3), (3,8), (5,8), (8,9), (8,10), (9,10).
```

The literal target is

```
(3838,3716,3712,4086,3718,3734,3766,3328,119,3456,2304,256).
```

Its auxiliary is new11 and again its only out-neighbour is old8. The
preparatory birth and every displayed adjacency-stage tournament are strong.
The preparation as a **whole** is H-neutral; neither the old flips nor the
permanent birth is asserted to preserve H individually.

|Branch|Target order|H|F_8|S_a|L|K|G|sigma|tau=G−F|
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
|gap2|11|1429|1493|1431|469|1024|1493|2|0|
|gap1|12|1429|1923|1430|899|1024|1923|1|0|

These are literal adjacency counts, independently reproduced by forward and
backward subset recurrences. In particular F is allowed to change: the
preparation sends root defect134 to64 or494. Applying C2 then gives

|Digit|Actual child H|Actual child F_8|Child root defect|Child order|
|---|---:|---:|---:|---:|
|+4|4291|4479|188|13|
|+2|4289|5769|1480|14|

This proves two actual nonzero next branches from the same physical H1429
output, including the changing inner/root defect. It proves no return of
either child to the same prepared class.

## 4. Native compilation and every intermediate

The audited native pair compiler is the one used in
[NATIVE_SIX §§2–4](../../multiport_grammar_20261103/component/NATIVE_SIX.md),
ultimately routed there to the complete pair-block guard construction.
For each old pair it chooses a directed Hamilton word in one of the two pair
orientations, temporarily births one helper, performs two Hamilton
complements, then deletes **only that helper**. The inverse recipe is used
when the selected word has the opposite orientation. Thus an old-pair flip
is compiled, not admitted as a primitive native operation.

The new standalone `surgery_native.py` copies that proved compiler and
independently replays every command with named adjacency sets. For each of
the two branches the receipt records six temporary helpers, twelve guarded
Hamilton complements, six helper-only deletions, all paths and cuts, all
intermediate matrices, and all intermediate `(H,F_8,F_8−H)` readouts. All18
native states in each branch are strong. Both C2 births are also strong.
The peak orders, including permanent C2 births, are13 and14 respectively.
The gap1 branch's preparatory permanent birth has `(H,F_8)=(6209,6939)`;
the subsequent surgery returns H to1429.

The finite checks suffice for these two fixed programs: every command and
every matrix is explicitly present. No claim about arbitrary native
programs preserving strongness is used. Invalid reversed or repeated-name
guards are separate negative controls. A Hamilton complement is not
silently treated as a fixed-port differential update; the complete graph
and guard word are retained throughout.

## 5. Unbounded response, and the exact neutrality boundary

For original clock t>=3 put k=t−2. The prepared input has TT_k, the fixed
six-core D, terminal u, root r, and the last two C2 names x,y. Write these
outside-prefix names as D0,...,D5,u,r,x,y in that order. The same six-flip
recipes are independent of k:

```
gap2: (D0,r),(D0,x),(D1,D2),(D3,D4),(u,y),(r,x);
gap1: birth a with sole out-neighbour r;
      (D1,D2),(D2,r),(D4,r),(r,x),(r,y),(x,y).
```

All adjacency stages are strong for every k>=1. Indeed the prefix is a
module throughout: root beats it and it beats every other outside vertex.
Replace it by one vertex. This is precisely the k1 matrix at the same
stage, which is strong. Lifting an outside return path through this module
connects any two prefix vertices even when their internal direction is
wrong, proving strongness of every larger expansion. This is a graph-level
statement; the finite native receipts above are specifically the k1 lift.

Delete r and the prefix from each final target. In the remaining physical
order the suffix and its incoming-root mask are

```
gap2: rows=(322,448,507,467,451,475,0,65,192), ports=209;
gap1: rows=(962,960,1019,963,971,987,768,832,512,0), ports=964.
```

Their exact native one-port signatures `(h,p,f,g,j,e,s)` are respectively

```
(5,1024,744,93,64,5,0),
(1,1024,961,527,494,1,1).
```

The inherited ordered-prefix formula `H=h+2^k f−j`, `F=h+2^k f` gives

```
gap2: H=744·2^k−59,   F=744·2^k+5;
gap1: H=961·2^k−493,  F=961·2^k+1.                         (3)
```

In each target the auxiliary has sole out-neighbour r. After deleting it,
every Hamilton path ending at r must place the entire prefix first: only r
has an outside edge into that prefix, and r occurs last. Thus its ending
count is the fixed suffix ending count2 or1, independent of k. The exact
one-ready/two-ready inverse of READY §1, or its elementary deletion proof,
gives S−H=2 or1. Equation(1) supplies G=F at every clock.

The original prepared H is `756·2^k−83`. Hence the two surgery changes are

```
gap2: −12(2^k−2),      gap1: 205(2^k−2).                   (4)
```

Both are neutral exactly at k=1, namely original t=3. This is an analytic
uniqueness statement, not a fit to four clock controls. Equivalently the
gap2 funnel family `93·2^(q−1)−59` meets the original family
`189·2^t−83` at the matching-order clock q=t+2 only for t=3. If q and t
are independent, the equation is `31·2^(q−1)−63·2^t=−8`;
Equal exponents would give 2-adic valuation at least8, impossible for−8.
For unequal exponents the smaller is3. If q−1=3<t, division by8 would
give `63·2^(t−3)=32`, impossible. Otherwise t=3<q−1 and division by8
gives `31·2^(q−4)=62`, hence q=5. Thus `(q,t)=(5,3)` is the unique
solution in q>=2,t>=3.

## 6. Exact finite-bank inverse and the remaining closure problem

`surgery_ready.py` exposes reusable pure functions `core_data`,
`tail_coefficients`, `ready_count`, `ready_target`, and `search`. The bank is
34 explicitly listed actual entry alternatives of orders1,3,4,5, with
entry sum1 or2 verified by adjacency DP. It is a declared finite bank, not
an assertion that these are all possible ready cores.

For core response data `(h,p2,e,F,G,b,C)` and a binary tail w ending in1,
the inherited exact ready-target count is

```
H(w)=L_2(w)(G+p2)+B_1(w)h−D_1(w)e−C−F−b.                 (5)
```

The direct supplier is [FUNNELS §§F3–F4](../../polarized_path_calculus_20261103/FUNNELS.md):
the funnel parent has `H_parent=L_2 G−C+B_1 h−D_1 e`,
`F_parent=L_2 p2−F`, and `s_parent=b`; its indegree-one birth has
`H_new=H_parent+F_parent−s_parent`. These identities give(5), with converse
preserving H. No numerical interpolation is used in the all-tail inverse.

All coefficients are calculated from actual core path tables. `B_1>=0`,
`0<=D_1<=1`, and a nonempty length-n tail has `L_2>=1+2^(n−1)`.
Consequently, for a fixed finite target set with largest target X, lengths
satisfying

```
(1+2^(n−1))(G+p2)−(e+C+F+b)>X
```

are impossible, at that length and every larger length. The selector is
therefore complete over **all finite tails of this bank**, without an
arbitrary tail cutoff. Empty tails are included exactly when the entry set
is proper. The271 short-word controls compare (5) directly with adjacency
counts and the actual gap; the full requested target search uses21,247
words and finds28 literal hits.

The principal targets are obtained without search-dependent relabelling:
gap2 uses C4 `(2,4,9,3)`, entries{1,2}, tail01111; gap1 uses the singleton,
entry{0}, tail001111101. Build the funnel, add an indegree-one newborn,
then take converse. This produces actual start-ready targets; the explicit
six-flip matrices above provide the retained-name implementation.

The child scout finds both gaps at order13 for H4291, and a gap1 target of
order14 for H4289. Its gap2 targets at order13 for H4289 are **too small**
to receive the fourteen-name child without deleting an old name. They are
scalar witnesses only. The present bundle does not compile the next round,
does not allow old-name deletion, and does not prove an unbounded selector
closed under both nonzero ternary digits. That is the remaining obligation.

## 7. Reproduction and status

Run from the worktree root:

```
python3 04-computation/nonzero_preparation_20261103/surgery_ready.py --output /tmp/ready.json
python3 04-computation/nonzero_preparation_20261103/surgery_native.py --output /tmp/native.json
```

Repeat with `python3 -O`. Both programs use always-active checks. Frozen
receipts in this directory are `ready.json` and `native.json`; full output
and hash pins are recorded in `MANIFEST.json`. The independent audit passed; the manifest retains its audited prepromotion proof pin. The finite bank, two actual programs,
unbounded response formulas, and OPEN closure obligation have separate scope.
