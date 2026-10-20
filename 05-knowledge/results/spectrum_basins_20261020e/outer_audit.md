# Independent audit of the outer-SCC black detour

**PASS.** The corrected [outer proof](outer.md), including its restricted
prefix lemma, recursive termination, local ledgers, and endpoint-hub
branch, proves the stated all-order result. Together with the inherited
singleton-padding theorem, every blue line incident to a nonstrong mixed
tournament has a black-only detour of length at most five. This does not
settle blue lines with two strong endpoints or global black connectivity.

## 1. Dependency and type audit

The needed rigid-strong alphabet is `{1,C3}`, not the stronger assertion
that every pure-blue strong tournament is C3 or R5. I independently read
Section 3 of the incoming [finite-atoms note](../merged_metagraph_finite_pure_atoms_20261020/FINITE_ATOMS_AND_RIGID_CLOSURE.md)
using `git show origin/main:...` and checked the primary
[Benhocine–Wojda publisher abstract](https://onlinelibrary.wiley.com/doi/abs/10.1002/jgt.3190070413).
Its tournament theorem supplies a spanning directed-cycle orientation
with one reversed edge at every order at least seven. That is exactly a
Hamiltonian path whose first vertex beats its last.

If raw `H=|Aut|`, the free automorphism action has one Hamiltonian-path
orbit. Camion supplies a closing path, so every Hamiltonian path closes,
contradicting that published bypass theorem for order at least seven.
Successive cuts of a Hamiltonian cycle also force its full cyclic
successor to be an automorphism. Tournament automorphisms have odd order,
excluding orders four and six. At order five this makes the tournament a
circulant; its four skew connection sets are equivalent, and the path
`(0,2,3,4,1)` in the `{1,2}` circulant is a bypass. Thus only C3 remains.
The published theorem is CITED, not reproved; the deductions and small
cases are checked. No uncertain pure-blue alphabet is used.

The outer theorem retains literal labelled path targets before passing
to converse-merged classes. Blackness is checked at the current source
of each move. The half-score test is sufficient only: a zero endpoint
score difference is never treated as a proof of blueness.

## 2. Prefix repair and recursive invariants

Let the endpoint-rigid path be `Q=(a=v0,b=v1,...,z)`, with `z->a` and
`b->z`. If a proper prefix rotation is another directed path ending at z,
endpoint rigidity makes that exact rotation an automorphism fixing z.
Its nontrivial cycle contains a and b, impossible because z distinguishes
them. This is a coordinatewise automorphism conclusion, not merely an
equality of path counts.

A transition `v_i->a->v_(i+1)` would make the left rotation of the prefix
through `v_i` legal. Hence a's outneighbors along Q form an initial
interval `v1,...,v_t`. If `t<=r-3`, the absence of a consecutive cyclic
triangle after a gives `v_t->v_(t+2)`. The right rotation through
`v_(t+1)` is then legal, another contradiction. Therefore a dominates
every vertex except z. The no-consecutive-triangle guard is essential to
this proof; the separate triangle ledger handles precisely that branch.
The earlier attempt to prescribe both endpoints by cutting a Hamiltonian
cycle is not used.

Every path ending at z now starts at a. Deleting a gives a bijection of
terminal paths, and fixing a identifies the two endpoint stabilizers.
Thus `S-a` is endpoint-rigid. If it were strong, closure of its terminal
path would give `z->b`, contrary to hypothesis. Its ordered SCCs therefore
have a proper final component D containing z. Varying only one preceding
component proves that it has one path orbit, so it is 1 or C3. Varying D
alone proves inherited endpoint rigidity there.

Each recursive D is a contiguous terminal suffix. All earlier right-hand
vertices beat `D-{z}`; they may lose to z. A preceding deleted vertex
beats the later suffix except z, and every intervening SCC beats the
whole later suffix, which preserves this invariant at every depth.
Consequently the terminal-ear and near-transitive ledgers have the exact
cross-arcs they require. A C3 block is consecutive and was already
handled by the first triangle branch. In the branch with no such block,
a terminal singleton yields the required near-transitive suffix with a
nonempty interior. The selected suffix strictly shrinks; the recursion
selects one final ledger, so its depth does not multiply the move bound.

## 3. Ledgers, score bounds, and the flexible branch

I checked each displayed path-edge cancellation and the successive
orientation changes in the triangle, terminal-ear, and near-transitive
ledgers. Empty L or R in the triangle ledger is allowed. The terminal-ear
middle C is nonempty because its strong component has order at least
four. Every unchanged interior block is traversed forward on every step,
so the boundary identities lift to arbitrary block length and internal
arcs.

For the three-step ledgers the subsequent endpoint-score differences are
`n-2-s_c-s_a` and `s_c+s_z+2-n`. Distinct r-vertex scores sum to at most
`2r-3`; with `n>=2r` these differences are strictly positive and negative.
The first move is black because a nonautomorphic terminal rotation cannot
be the restriction of an automorphism of the whole source SCC.

For the near-transitive suffix, `k>=1`, `r>=k+3`, `q<=r-k-2`, and
`0<=t<=k-1`. The four differences

    n-k-q-3,  k+3-n,  n-k-t-3,  t+4-n

have signs `+,-,+,-`: the first is at least `n-r-1`; the third is at
least `n-2k-2>=4`; and the other two are at most `k+3-n<0`. The first
rotation is black because the suffix has unequal internal scores, so a
full cyclic shift is not its automorphism. These bounds do not assume
that earlier right-hand vertices beat z.

For a flexible terminal endpoint, change only the final SCC path to a
different orbit while retaining its last vertex. The new whole path is
black, since two positional anti-reversals would otherwise identify the
terminal path orbits. The global first a, second b, and last z stay fixed.
Rotate both complement partners to begin at b: they reach the same
labelled hub `U xor {az,ab}`. The two remaining steps have opposite
nonzero score differences `n-d_U(a)-d_U(b)` and its negative. Indeed the
initial SCC has r vertices, so both endpoint degrees include `n-r`
external wins and their internal scores sum to at least one; with
`n>=2r`, their sum exceeds n. This independently checks the inherited
common-hub use in the outer theorem.

The final SCC-size split is exhaustive: size one uses the previously
audited singleton-padding theorem; size three is C3, whose single path
orbit forces a mixed middle and admits the mixed-middle common hub;
size at least four uses this new result. There is no strong size-two SCC.

## 4. Independent exact controls

[outer_audit.py](outer_audit.py) imports no producer module. It enumerates
all standard-Hamiltonian-path presentations of orders four through six:
8, 64, and 1,024 matrices. Independent mutual reachability finds 5, 50,
and 903 strong matrices. Direct complete terminal-path comparison finds
3, 10, and 29 endpoint-rigid presentations.

Each is padded with the independently checked blue centers of orders
0, 1, 3, and 5. The checker verifies 3,664 flexible common-hub detours,
168 endpoint-rigid recursive detours, and 3,200 general triangle ledgers.
Recursive controls include depth two. Every check requires directedness
of every path, literal blackness at its current source, the exact final
labelled matrix, and the displayed numerical score formulas. This is a
finite control of the universal proof, not its all-order justification.

Reproduce from the repository root:

    python 05-knowledge/results/spectrum_basins_20261020e/outer_audit.py
    python -O 05-knowledge/results/spectrum_basins_20261020e/outer_audit.py

Both runs pass and write the same LF-normalized
[receipt](outer_audit.json), SHA-256:

    aeb95eb713e28bfcbd3de278644d61ccfe9c1446aa0a5dc643467439cfa853d2
