# Full-count port fibres, a growing-defect return, and its continuation boundary

**PROVED:** the exact full-count port-fibre classification; a universal
module-substitution identity giving a same-target two-arc return; its
specialization from unbounded defect to defect one; and the all-size failure
of the next return if only the latest root's predecessor may change.
**FINITE-EXACT:** the explicitly bounded catalogue of 225 substitution
charts, 78 with strong parents. **OPEN:** unlimited target-preserving
continuation and Chan–Pak. This file does not claim either open problem.

## 1. Inheritance and the six-concept board

The closest proved mechanisms are the endpoint-cover birth calculus in
[polarized path calculus](../../polarized_path_calculus_20261103/CALCULUS.md),
the actual five-core and indegree-one source work in
[positive compiler](../../chanpak_five_core_bridges_20261110/algebra/POSITIVE_COMPILER.md),
and the ordered-run substitution method used in
[target returns](../../berggren_target_return_20261111/reset/TARGET_RETURN.md).
The directed-cycle-opening formula is inherited from
[triples/PROOF.md section 9](../../triples_ports_sieve_20261112/squares/PROOF.md).
Its canonical hostile example has one fixed cycle-opened path but full
counts 9,3,17,5,17,5,23,5 across its eight star completions. Controlling that
one sector does not control the total. The corrected near miss is the
425 rebase: its positive third-run channel prevents the proposed all-size
identity. The least-used sidecar here is the complete endpoint-cover
response, including every interior split.

Live concepts are (1) full Hamiltonian count, (2) physical root identity,
(3) indegree-one defect, (4) five-core module substitution, (5) exact
ordered-run covers, and (6) a reusable return after a target-changing
excursion. The new positive result couples the first five. Section 7
records exactly where the sixth still needs a parent change.

The H9 pair obtained from the four-vertex parent `(2,4,9,3)` is recovered
old five-core work, not a new seed: entries 2 and 3 have equal full counts
9 and defects 1 and 2. Inflating vertex 0 gives their entry-score
difference `2^k-2`, so that seed alone is not uniform. The construction
below is a different, coefficientwise neutral chart.

## 2. Definitions and the complete full-count law

For a finite nonempty tournament T, let H(T) count directed Hamiltonian
**paths**, with H(empty)=1. Let s_a(T) count such paths starting at a.
Let E_a(2;T) count ordered pairs (P1,P2) of disjoint directed paths covering
T, where P1 is nonempty and ends at a, and P2 may be empty. Put

    Q_a(T) = E_a(2;T) - s_a(T).

Append a named vertex z with a→z and z→v for every old v≠a. Denote the
result B_a(T), always rooted at the same physical z. Then

    H(B_a(T)) = H(T) + Q_a(T).                         (1)

This is a bijection-level full-count identity. Paths starting at z count
H(T)-s_a(T). All other paths have the consecutive pair a,z. Removing z
splits such a path into a first path ending at a and a possibly empty
second path. Conversely, insert z at that split. This sector counts
E_a(2;T). The two sectors are disjoint and exhaustive. In particular,
interior cycle-opening contributions are not silently discarded.

For a rooted tournament (U,z), write

    delta(U,z) = E_z(2;U) - H(U).

For this birth,

    delta(B_a(T),z) = s_a(T),
    s_z(B_a(T)) = H(T)-s_a(T),
    H(B_a(T)-z) = H(T).                                (2)

Indeed an ordered two-cover ending its first bin at z either has that bin
[z], contributing H(T), or has penultimate a, contributing E_a(2;T) after
deleting z. Thus E_z(2;B_a(T))=H(T)+E_a(2;T), proving (2).
For any indegree-one root, these identities also give
`delta=H(U-z)-s_z(U)`.

### Complete classification of a fixed-parent slice

Fix T and the named z. Partition V(T) into level sets of Q_a(T). By (1),
these are **exactly** the full-count fibres among all extensions in which
z has indegree one. Each fibre is a clique under the operation that changes
its predecessor from a to b: reverse a–z and b–z and nothing else. Its
vertex labels are the defects s_a(T). These are actual two-arc involutions,
not arithmetic substitutions of equal numbers. Every old named vertex,
every old arc, and the physical marked z are retained.

This classification is all-order but only for a fixed parent and an
indegree-one apex. It does not classify arbitrary parent rewrites or all
vertices of the tournament metagraph. It deliberately allows marked
response to change; full-response-frozen fibre equivalence would be too
restrictive for a defect-changing return.

## 3. Exact run coefficients and arbitrary modules

Choose a non-port vertex e of a fixed parent. Replace e by an arbitrary
nonempty tournament M, with homogeneous comparisons to all other vertices.
A Hamiltonian path or ordered two-cover breaks its visits to M into
ordered, nonempty, maximal runs. Let K_j(M) count ordered j-tuples of
nonempty directed paths partitioning M. Each skeleton with j runs lifts
in exactly K_j(M) ways. Thus every observation has the exact form

    observation = sum_j A_j K_j(M).                   (3)

For H or a specified start there are at most n runs when the base parent
has n vertices; for an ordered two-cover there are at most n+1. The finite
skeleton dynamic program in `inflation_probe.py` enumerates all of them.
It retains an endpoint, the used fixed vertices, and the run count. For
E_a(2), split the fixed vertices between two ordered bins, constrain the
first endpoint to a, and convolve their run counts. Empty second bins are
included. This proves the compiler's finite recurrence and its degree
bound; the formula is not inferred from a finite fit.

For M=TT_k, a path on any chosen subset has its unique increasing order,
so `K_j(TT_k)=j! S(k,j)`, where S is a Stirling number of the second kind.
These functions are independent on positive integers k: evaluating at
k=1,2,... is triangular with positive diagonal. Consequently a vanishing
complete run coefficient vector is necessary and sufficient for an
all-k transitive-block identity, and is sufficient for an identity for
**every** nonempty tournament module M.

Write h=H(M), p=P_M(2), where P_M(2) counts ordered two-covers with either
bin allowed empty. Then K1=h and K2=p-2h. Higher K_j retain information
that h and p alone generally lose; their absence below is a proved
property of this particular interface, not a universal two-state claim.

## 4. A five-core return from growing defect to defect one

Take the outgoing-bit parent

    P = (78,12,24,16,35,79,30),

replace vertex 1 by M, and call the result P[M]. Compare its entries
0 and 4. The complete nonzero run coefficients are:

| Observation | K1 coefficient | K2 coefficient |
|---|---:|---:|
| H(P[M]) | 45 | 21 |
| s_0(P[M]) | 13 | 6 |
| s_4(P[M]) | 1 | 0 |
| E_0(2;P[M]) | 65 | 30 |
| E_4(2;P[M]) | 53 | 24 |
| Q_0 = Q_4 | 52 | 24 |
| Either child H | 97 | 45 |

Every higher coefficient through the proved run bound is zero. Therefore,
for **every nonempty tournament M**, the two actual children satisfy

    H(B_0(P[M])) = H(B_4(P[M])) = 7h + 45p,
    delta(B_0(P[M]),z) = h + 6p,
    delta(B_4(P[M]),z) = h.                            (4)

The only changed arcs are 0–z and 4–z. This is a full-count-preserving
involution with an unchanged named root; it can be toggled in either
direction any number of times at that target. It returns the growing
extra defect 6p to the module's baseline defect h.

For M=TT_k, k≥1, h=1 and p=2^k. Thus (4) becomes

    H = 45·2^k+7,
    delta_high = 6·2^k+1,
    delta_low = 1.                                    (5)

The first values are H=97,187,367,727 and high defects 13,25,49,97.
This is an unbounded-defect-to-ready return at an absolute target, with
no count normalization. On the low side the full useful state is

    H=45D+7, E_z(2)=45D+8,
    s_z=21D+2, H(U-z)=21D+3,  D=2^k.

Both parent and children are strong. The parent core P is strong, and
nonempty module substitution preserves reachability between its modules;
vertices inside any module can reach one another by a path through the
strong outer core. In either child, old vertices reach z through its
predecessor, and z reaches the old strong parent through any other vertex.

### The connection to five is structural

P compresses to the five-vertex tournament

    C = (22,4,8,3,14),       H(C)=9.

Its modules, in the original P labels, are the ordered pair (5,0), the
singleton 1, the ordered pair (2,3), singleton 4, and singleton 6.
Thus the whole family is the explicit substitution
`C[TT2, M, TT2, 1, 1]`, with entries the last vertex of its first TT2
and its fourth module. This is a genuine map preserving homogeneous
comparisons. No connection to other five-element collections is inferred
merely from their cardinality.

## 5. Complete finite chart search and hostile controls

Input is the complete ordinary-parent bank through order 8 in
[`returns/port_collisions.json`](../returns/port_collisions.json), produced
by the independent return-search lane. For every unequal-start neutral
pair and every non-port vertex chosen for substitution, the complete
run vectors were compared. There are exactly 3,615 such charts, of which
225 are neutral for every nonempty substituted module; 78 have strong
parents. These are charts, not isomorphism classes of infinite families.
This is complete only in the declared chart universe.

Among the strong charts, four have a constant defect 1 or 2 in the
transitive specialization. They yield three distinct (H,defect) formulas:

| H | Larger defect | Smaller defect |
|---|---|---|
| 45·2^k+7 | 6·2^k+1 | 1 |
| 75·2^k+7 | 10·2^k+1 | 1 |
| 90·2^k+7 | 12·2^k+1 | 1 |

The third row is a shift of the first value family; graph charts are not
silently identified by equal counts. The receipt lists each parent,
entries, substituted vertex, complete coefficients, and source-bank hash.
No completeness claim beyond one substituted non-port module in that
finite parent bank is made.

Independent controls in `audit_counts.py` check all 225 charts at k2 by
literal subset path dynamic programming, including both full children and
their defects. The selected family is checked at k1 through k8. Its first
two parents are also enumerated independently by permutations and every
ordered two-bin cut. Nontransitive hostile modules C3, C4, and the H9
five-core verify (4) directly; this tests the stronger statement that
transitive run coefficients were not mistaken for scalar-count sufficiency.

## 6. What the full-count return repairs in cycle opening

The cycle-opening identity isolates a selected internal sector. Equation
(1) sums all start, end, and interior sectors before comparison. Equality
of Q is therefore an exact criterion for a physical compensation of the
entire birth count. The five-core family exhibits such compensation with
positive growing defect, and (3) explains why it survives arbitrary module
size and internal structure. It does not assert that every cycle opening
can be compensated, nor that arbitrary formal differences of counts have
literal tournaments realizing them.

The useful new object is a **port fibre with defect labels**: its level is
H(T)+Q, its edges are actual two-arc toggles, and its labels record the
source readiness exposed by each entry. This supplies known local pieces
of the metagraph across size; it does not classify all of that metagraph.

## 7. Exact all-size boundary after two ordinary births

Start on the low side of (5), rooted at z. Append x beaten only by z,
then y beaten only by x. The literal target is

    H = 135·2^k+23,       delta_y = 24·2^k+5.           (6)

Let X_k=P[TT_k]+z+x be the parent of y. Its H is 69·2^k+13, and the
existing entry x has Q_x=66·2^k+10. All fixed alternatives have the
following complete nonzero run coefficients for Q_a-Q_x:

| Old base label a | K1 | K2 | K3 |
|---|---:|---:|---:|
| 0 | 348 | 156 | 0 |
| 2 | 1224 | 1156 | 346 |
| 3 | 1638 | 1372 | 346 |
| 4 | 68 | 30 | 0 |
| 5 | 176 | 78 | 0 |
| 6 | 524 | 261 | 0 |
| z (base 7) | -36 | -18 | 0 |
| x (base 8) | 0 | 0 | 0 |

For every k≥1, K1(TT_k)=1 and all K_j≥0. Hence every alternative fixed
entry has strictly different Q. This fixed-entry conclusion even holds
for every nonempty module M.

To include every internal vertex of TT_k, write it as
`TT_l → v → TT_m`, where l,m≥0 and k=l+m+1. The exact two-block compiler
in `two_block.py` expresses Q_v-Q_x in
`A_i(l) A_j(m)`, with `A_i(t)=i!S(t,i)`, including
`A_0(t)=1` for t=0 and zero otherwise. Its coefficient matrix is

| i \ j | 0 | 1 | 2 |
|---|---:|---:|---:|
| 0 | 1026 | 1948 | 903 |
| 1 | 3122 | 5988 | 2847 |
| 2 | 3134 | 6116 | 2982 |
| 3 | 1038 | 2076 | 1038 |

All other coefficients vanish. The four positive entries at
(0,0),(0,1),(1,0),(1,1) establish strict positivity in the four cases
l=m=0, l=0<m, m=0<l, and l,m>0. This includes both ends of every inflated
block; no growing internal-port set is omitted.

Therefore the Q-level containing x in X_k is a singleton for **every**
k≥1. The next target (6) has no different indegree-one apex extension
with this exact fixed parent. This is a mechanism-specific obstruction,
not an obstruction to a ready representative at (6). Parent substitution,
additional arc changes, and retained-vertex replacement remain available;
the parallel compiler lane proves a larger physical return below.

The two-block recurrence permits zero-size blocks and is independently
checked on all 16 pairs l,m∈{0,1,2,3}. A finite radius-two parent-arc search
also found no defect-one repair at k1 among the 667 parent graphs examined
for non-inflated entries; that auxiliary scope is recorded in
`repair_search.json`, not promoted to an all-order obstruction.

### A parent-changing return survives this obstruction

The [uniform module return](../compiler/MODULE_RETURN.md) now constructs a
ready representative at the successor target for every nonempty M. Its
full count is `23h+135p` and defect is h. In the transitive case it restores
defect one at exactly `135�2^k+23`. It keeps the same physical newest root
y, its sole predecessor x, every internal M arc, and every root-star arc;
it changes six fixed parent arcs and three homogeneous module boundaries.
Thus section 7 excludes only predecessor toggles with a frozen parent,
while this positive result supplies the necessary physical compensation
outside that slice. The [separate sparse migration](../returns/MODULE_MIGRATION.md)
is another exact realization that preserves y and changes its predecessor.
Both have been independently checked here on five modules; see
[the root-return audit](ROOT_RETURN_AUDIT.md) and `cross_audit.json`.

## 8. Reproduction, status, and stopping point

From the repository root:

```
python 05-knowledge/results/cycle_opening_six_return_20261112/counts/inflation_probe.py
python 05-knowledge/results/cycle_opening_six_return_20261112/counts/two_block.py
python 05-knowledge/results/cycle_opening_six_return_20261112/counts/audit_counts.py
python -O 05-knowledge/results/cycle_opening_six_return_20261112/counts/audit_counts.py
python 05-knowledge/results/cycle_opening_six_return_20261112/counts/cross_audit.py
python -O 05-knowledge/results/cycle_opening_six_return_20261112/counts/cross_audit.py
```

The receipt includes script/source hashes, the exact tested universes,
positive examples, nontransitive controls, and the hostile continuation
certificate. The proof is the explicit birth bijection and finite
run-skeleton recurrence; literal enumerations check the implementation.
The full-count port fibre is exactly classified at every fixed parent,
and the selected family works at every size and for every nonempty module.
One uniform target-changing excursion now has a proved physical ready
return. The unresolved step is indefinite closure after further excursions.
Repeatedly toggling a single fibre, or performing that one excursion, is not
a proof of indefinite closure.
