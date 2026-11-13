# Uniform target return by moving the repeated module

Status: **PROVED**, with two independently implemented finite module-word certificates and literal **VERIFIED** cases k=1,...,5. An infinite family of exact source returns is obtained. Repetition of the same return on its output remains **OPEN**.

## Result and literal operation

Let P=(78,12,24,16,35,79,30). Replace vertex 1 by a transitive tournament TT_k, k>=1, retaining the old vertex as its first element. Append a root z beaten only by vertex 4. The resulting strong source has

    H=45*2^k+7,  b=1,  d=21*2^k+3.

The full-count lane separately proves the same-parent two-arc switch between entry ports 0 and 4. At unchanged H it changes the defect from 6*2^k+1 to 1. Thus this source already comes from a genuine unbounded-defect repair.

Apply the usual two successive sole-predecessor births x,y, first x beaten only by z and then y beaten only by x. The full count becomes 135*2^k+23, and the marked y has defect 24*2^k+5. We now repair its defect back to one while preserving that full count and the same physical marked vertex y.

At k=1 the repaired target is

    Q=(14,12,8,64,175,79,535,367,127,447),

root 9, with sole predecessor 6. Its uniform family replaces **vertex 5**, not vertex 1, by TT_k. This change of the repeated vertex is essential.

Use source names 0,...,6 for the original P vertices, clone names 7,...,k+5, then z=k+6, x=k+7, y=k+8. Starting with the raw two-birth output, toggle the following arcs simultaneously:

- the six fixed-name pairs (0,6), (2,4), (3,4), (3,6), (6,y), (x,y);
- for each clone j=7,...,k+5, the three pairs (0,j), (1,j), (6,j).

The output is exactly Q with vertex 5 replaced by TT_k, under the label map Q's 7,8,9 to z,x,y and Q's additional clones to 7,...,k+5. There are 3k+3 changed arcs. No existing vertex is deleted; the operation retains y literally and changes its unique predecessor from x to 6. Both endpoints are strong. The assertion is macro-level equality of full counts; individual intermediate arc flips are not asserted count-neutral.

The repaired output has

    H=135*2^k+23=3(45*2^k+7)+2,
    b=1,
    d=27*2^k+27.

The frozen vertex-1 inflation of Q instead has H=132*2^k+29, agreeing only at k=1. Its discrepancy is -3(2^k-2). Moving the module, rather than extrapolating the base witness, is what repairs this quantified failure.

## Why the identities hold at every height

Collapse consecutive visits to the replaced module into a single module letter. A Hamilton path of the expanded tournament determines a word using every outside vertex once and the module letter j times, never twice consecutively. Conversely such a word, together with an ordered partition of the module into j directed paths, gives exactly one expanded Hamilton path. This is a bijection.

For TT_k, the number of ordered j-path partitions is j!S(k,j). Independent exact word enumeration gives the following coefficients in this basis; all j>=3 coefficients vanish:

| quantity | j=1 | j=2 |
|---|---:|---:|
| source H | 97 | 45 |
| source root deletion d | 45 | 21 |
| source defect b | 1 | 0 |
| migrated target H | 293 | 135 |
| migrated target root deletion d | 81 | 27 |
| migrated target defect b | 1 | 0 |

Since 1!S(k,1)=1 and 2!S(k,2)=2^k-2, the claimed formulas follow immediately. Defect is counted by paths starting at the entry port in the root-deleted parent, using the indegree-one birth identity. Both independent compilers enumerate finite collapsed words, not a polynomial guessed from sampled k.

There is a stronger arbitrary-module version. For any nonempty tournament module M, let h_M=H(M) and P_M(2) count ordered two-bin directed path covers with empty bins allowed. The same word bijection gives

    H(source)=7*h_M+45*P_M(2),     b(source)=h_M,
    H(target)=23*h_M+135*P_M(2),   b(target)=h_M.

Thus H(target)=3H(source)+2h_M, retaining the module's exact defect contribution. TT_k has h_M=1 and P_M(2)=2^k, giving the ready return above. This is a transfer law, not a statement that the target quotient is the same as the source quotient.

## What has and has not closed

The construction now composes two uniform operations: same-count unbounded-defect repair at H=45*2^k+7, followed by a full-target-preserving source return at H=135*2^k+23. The second quotient differs from the first. We have not produced a uniform successor source at 405*2^k+71, nor a return of Q to the original quotient/interface family. Hence this does not yet give indefinitely repeated affine target closure or prove Chan-Pak.

The exact isomorphism from Q names to the root lane funnel names is 0->4, 1->5, 2->6, 3->7, 4->0, 5->3, 6->8, 7->1, 8->2, 9->9. It carries module vertex 5 to the funnel module vertex 3, root 9 to root 9, and predecessor 6 to predecessor 8, so it extends over every substituted module, including nontransitive ones. This is a graph isomorphism, not merely equality of count observers.

The root lane independently recovered the same second family via a C3 funnel and produced another label map preserving both y and its predecessor x, at a different arc cost. That is a distinct retention guarantee, recorded in `../compiler/MODULE_RETURN.md`; this note claims the smaller 3k+3 arc recipe retaining y and retargeting its predecessor.

## Reproduction and finite search scope

    python -B 05-knowledge/results/cycle_opening_six_return_20261112/returns/module_migration.py
    python -B 05-knowledge/results/cycle_opening_six_return_20261112/returns/sparse_repair.py --radius 4
    python -B 05-knowledge/results/cycle_opening_six_return_20261112/returns/relocate_modules.py

`module_migration.json` includes the literal source/raw/repaired states and exact changed-arc lists for k=1,...,5. `sparse_repair_radius4.json` examines all 1+36+630+7140+58905 subsets of at most four arcs of the 9-vertex raw parent, every latest-root sole predecessor, target H=293, and defects at most two. This universe has 112 hits (including defect zero). Exactly three positive-small-defect hits are strong; they all have b=1 and frozen-module H=132*2^k+29. Relocating the target module yields H=135*2^k+23 with b=1 in all three. `repair_module_relocation.json` records these three target-module choices. These are bounded neighborhoods and module choices, not an exhaustion of all ready tournaments.

The radius-three frontier has no strong positive defect <=2 hit. The first successful source-return route at radius four is therefore minimal within this precise parent-flip radius and latest-root retarget search, before the module relocation. This does not claim a universal minimum number of changed arcs for arbitrary constructions.

## Bounded next-step stopping result

Starting with Q at k=1 (H=293,b=1), the next two-birth target is H=881. Delete the latest newborn to obtain the 11-vertex raw parent. `next_repair.py --radius 3` examines every subset of at most three of its 55 arcs, then every possible sole predecessor of a replacement latest root. Across 1+55+1485+26235=27776 parent variants, there are **zero** H=881 children with defect 1 or 2, even without requiring strongness. Thus this exact local repair grammar fails before any choice of target inflation module is considered. The receipt is `next_repair_radius3.json`; radius two is retained as an independently run smaller control.

This is a finite obstruction at one specified base source and one arc radius. It does not exclude larger sparse repairs, vertex deletion/replacement, a different H=293 source, modules with several boundary types, or a nonlocal uniform next return. The positive 45*2^k+7 to 135*2^k+23 construction and the negative next-neighborhood result have different quantifiers.

    python -B 05-knowledge/results/cycle_opening_six_return_20261112/returns/next_repair.py --radius 3
