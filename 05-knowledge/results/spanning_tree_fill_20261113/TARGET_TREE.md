# A decreasing target tree with five finite bypasses

**PROVED arithmetic organization and finite witnesses; CONDITIONAL
tournament lift. Chan–Pak remains OPEN.** This note transfers the explicit
macro-parent from the enlarged Collatz graph to the desired odd target
set. It does not assume that those targets are already Hamilton counts.

## 1. The exact map and its scope

Let D be the decreasing integer parent in [the Collatz proof](collatz/PROOF.md).
Map a positive integer m to the odd label H=2m−1. The induced parent p(H)
for H>1 is

    (H−1)/2       if H=3 mod4;
    (H−2)/3       if H=5 mod12;
    (H−4)/3       if H=1 mod12;
    (2H−3)/3      if H=9 mod12.                         (1)

The H=1 root is excluded from the third clause. These cases partition all
odd H>1. Every output is positive and odd, and

    (p(H)+1)/2 ≤ (2/3)(H+1)/2.

Consequently (1) is a rooted tree on all positive odd labels. It is an
arithmetic address system, not yet a graph constructor for Hamilton counts.

Remove the labels7 and21. Whenever p lands on a removed label, continue
the parent map until it leaves that pair. This still yields a smaller
parent and a spanning tree on

    O={positive odd H : H is not7 or21}.                (2)

Only five edges change: the old children of7 were15,25, and its parent
was3; the old children of21 were33,43,65, and its parent was13. The new
edges are therefore

    3→15, 3→25,       13→33, 13→43, 13→65.            (3)

All other child edges have one of four shapes, with their exact guards:

|Parent h|Child H|Guard in the target tree|
|---|---|---|
|h|2h+1|h odd; omit child7|
|h|3h+2|h=1 mod4|
|h|3h+4|h=3 mod4|
|h|(3h+3)/2|h=5 mod8; omit child21|

Rows whose parent is removed are replaced by (3). This tree has no claim
of novelty for the four affine expressions themselves. Two-birth laws,
binary shells and signed Collatz work already contain them in different
guarded settings. The contribution is their exact parent partition and
finite bypass bookkeeping derived from the present macro-descent.

## 2. What a perfect-fill proof would have to lift

Suppose there is a class of actual source certificates (tournament,
retained marks, admissibility data) containing a certificate at1, and for
every edge of the tree in(2) there is a procedure that accepts every
admissible parent certificate produced by the induction, or first converts
it by a certified procedure, and returns an admissible certificate at its
child. Then all labels in O are Hamilton counts, by strong induction.

This is a useful sufficient construction criterion. Adjoining the next
allowed target maintains exactly one selected witness at every earlier
target; its parent has already been supplied. A final tournament need not
contain its parent's vertices as an induced subgraph. The actual
constructors may replace a core, change arcs, or use a finite preparation.
The maximum order of the saved witnesses is a separate nondecreasing size
budget. No uniform order bound follows from the arithmetic tree alone.

The five bypass targets themselves are already attained. The script
retains literal witnesses:15 as C3→R4,25 as R4→R4,33 as R7,65 as R8,
and43 from the inherited complete order-seven bank. Here R_n is the
transitive tournament with its first/last edge reversed. Literal subset
counts and permutation counts independently verify these five graphs.
These base witnesses do not automatically restore arbitrary live marks
from a parent certificate. A typed induction must either accept them as
fresh source roots or prove the required port conversion.

## 3. The current interface obstruction lands on an actual tree branch

There is already a stronger, small obstruction to an induced-growth lift:
the tree contains13→27, but no H27 tournament contains any H13 tournament
as an induced subtournament, at any order. This is the inherited
[November2 induced obstruction](../metagraph_primegame_20261102/INDUCED_OBSTRUCTIONS.md),
recovered and independently audited in the concurrent
[fill-tree packet](../fill_tree_matrix_20261113/INDUCED_BOUNDARY.md).
Even arbitrary H13-preserving preparation followed solely by vertex births
cannot implement that edge. A successful implementation must permit changes
to old arcs during the transition, loss of old vertices, or a different
arithmetic parent. Thus the broad certificate formulation in Section2 is
essential; a pure-birth version is false.

For an indegree-one root with state(H,s,q,F), let δ=q−s=F−H. The inherited
two births NR give the balanced count3H+2δ. Thus the second and third
rows suggest δ1 andδ2 source certificates. But the arithmetic guards do
not supply those certificates or rearm them at the child.

The attained count159252697 is1 modulo4, so the tree assigns its ternary
child477758093. Yet the previous packet's complete inverse excludes every
indegree-one unit-ready source at159252697. The NR/δ1 implementation of
that edge therefore fails even though its arithmetic guard passes.
The smaller known unit-ready hostile1997 is also1 modulo4.

At53084231, which is3 modulo4, a defect-two source does exist, and its NR
exit lands exactly on159252697. This is one physically realized edge of
the target tree followed by a genuine rearming obstruction. It explains
why growing ordinary multiplicity alone cannot prove the induction.
Both the positive defect-two source and the all-order ready exclusions
are inherited from [AT_TARGET](../target_directed_cores_20261112/AT_TARGET.md).

The next task is a certificate language that admits changed interfaces or
selective source replacements at these children. The general adjacency
cube supplies all tournaments and the arithmetic tree supplies all target
labels, but a surjective correspondence between them is still the missing
mathematics. No convergence assertion for the absent-count set is used
as a premise.

## 4. Reproduction

`python3 -B 05-knowledge/results/spanning_tree_fill_20261113/target_tree.py`
checks the parent guards, complete bypass list and decreasing depths on
the declared finite head. It independently counts the five literal base
witnesses in two ways. The infinite tree and conditional lifting statements
are proved above rather than inferred from that finite head.
