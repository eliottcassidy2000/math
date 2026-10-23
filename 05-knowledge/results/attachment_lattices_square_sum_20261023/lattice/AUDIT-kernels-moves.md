# Independent proof audit of the larger-kernel and coupled-move lanes

23 October 2026. Reviewer: lattice-attachments agent, independently of the
author of [KERNELS.md](../KERNELS.md) and [MOVES.md](../MOVES.md).

**PASS, scoped proofs.** No assertion about arbitrary source recovery,
unbounded-cell efficient canonicalization, or native electrical weights is
needed by either result.

## Kernel factorization

For every incident block stabilizer, its one-sided kernel is normal in its
projection. Intersecting these kernels remains normal in the intersection
of the projections. The product kernel is a subgroup of the full compatible
group because each of its factors preserves every incident block while all
other coordinates are identities. Hence every quotient-compatible tuple
has the full product of kernel sizes as its number of lifts. This statement
does not require separately chosen edge witnesses to agree: the quotient
tuple is explicitly one shared element at each cell.

In the source/target comparison, a nonempty intersection of projected cosets
is a coset of the intersection of their subgroups. Translating by a chosen
actual representative therefore gives the stated domains, with all
constraints still checked on the common translated variables. The proof
retains possible global emptiness after local restriction.

The source/target coset direction is correct: if `f0` maps source to target,
all other maps are `f0*h` with `h` in the SOURCE stabilizer. After choosing
the projected representative `a_i`, write `f_i=a_i*p_i`; the source kernel
still acts by right multiplication `p_i*k_i`. Because it is normal in
`P_i`, a quotient tuple has exactly the claimed lifts. An observer descends
to the quotient exactly when `P(f*k)=P(f)` on every such fibre. Literal
owner/path observations need not satisfy that condition. Replacing their
fibre enumerator by its cardinality would be invalid; the theorem expressly
retains the enumerator or the actual lift instead.

The marking lower bound is necessary because an unbroken local kernel acts
while fixing every other cell. After local bases have killed their product,
the remaining stabilizer injects into the quotient. It has odd order; each
nontrivial actual orbit has at least three vertices. Successively fixing
one owner in a nontrivial orbit gives the logarithmic upper bound. This is
a base bound, not an assertion that ordinary equitable refinement always
recovers a general rigid tournament.

## Nested ternary tournament family

The displayed common-out-neighbor formula is correct. At the first
differing digit, if `x->y`, the common outgoing vertices are the outgoing
neighbors of `y` within its child; higher equal-prefix ancestors contribute
their common outgoing sibling blocks. This is the precise interpretation
of the split-level count. Its distinct values recover every prefix depth,
and each ternary quotient is `C3`, giving the full iterated wreath action.

Equality-prefix cross blocks require equality of the actions on the
corresponding prefix level. Counting source-prefix rotation choices by the
components of each threshold graph therefore gives
`3^(sum_r 3^r*c_(r+1))`. The coordinates have a semidirect composition law;
the proof correctly does not replace that law by coordinatewise addition.

Each bottom-triple prefix over each deepest-support component has an
independent rotation with disjoint actual support. It needs an actual mark,
giving the lower bound `3^(d-1)c_d`. Selecting one marked leaf per bottom
triple in a representative of each deepest component kills the entire
representative-cell action. It also gives distinct directed profiles:
within a bottom triple use its mark; across bottom triples use a mark in
the third branch at their lowest common ancestor. Identity matching blocks
then propagate singleton refinement. Thus the same upper bound proves
both base minimality and the claimed directed-refinement sufficiency.

The special `d=2` reduction from `3m` to `3c_2` depends on actual identity
matching edges at the deepest level. Coarse connectedness alone does not
provide it. Twisted transports are explicitly outside the count-only rule.

## Wide moves and finite-state compression

The native triangle-chain path formula is valid: the last entry in cell
`i` is `t-i+2`, equal modulo three to the first entry `t-(i+1)` in the next
cell. The block matching thus contains every cross step. Shared rotations
are exactly three, independently of the number of cells.

In the equality/one-hot matrix kernel, all cell triples agree and sum to
zero. A conformally primitive nonzero zero-sum triple is exactly one unit
transfer, so every primitive full move has support and one-norm `2m`.
The root-state description still has three possibilities. This is a valid
encoding-relative support lower bound, not an intrinsic hardness claim.

For the weighted block summary, composition transports the root state
through the first permutation before evaluating the second cost. The
displayed formulas are associative by direct substitution. Retaining the
decomposition and lengths allows prefix-address queries; keeping only a
final permutation and total cost does not. The theorem appropriately
requires source-dependent alternatives to be part of a cache key.

## Cross-lane boundary

Our [cyclic CRT elimination](README.md) solves a special unweighted
comparison/counting problem with scalar kernel multiplicities. The full
ternary wreath groups in KERNELS.md are nonabelian. Their exact group
factorization, actual ports, and conditional observation enumerators cannot
be replaced by two congruences per attachment without an authenticated
cyclic reduction. Conversely, an unbounded primitive move in one-hot
coordinates does not obstruct a small-state arithmetic or root-state
description. These two results are compatible and explain different
sources of information cost.

## Checker audit and executed replay

I read `check_kernels_moves.py` and independently replayed both

    python -B 05-knowledge/results/attachment_lattices_square_sum_20261023/check_kernels_moves.py
    python -O -B 05-knowledge/results/attachment_lattices_square_sum_20261023/check_kernels_moves.py

Both returned **PASS, 4,949 checks**, with identical frozen receipt content.
The computational universe is:

- Every common-out-neighbor pair at depths one through four (3,9,27,81
  vertices), distinguishing mark profiles at those depths, and 30 seeded
  actual tree-action samples per depth.
- All 81 source-prefix actions at depth two; all 6,561 action pairs on each
  of the uniform/coarse/fine blocks; all `81^3` chart triples grouped across
  all 27 edge-depth assignments on three named cells. These are full exact
  enumerations within this declared native family, not all nine-tournaments.
- A mixed coarse/fine path with 2,187 frames, local kernel sizes 27,1,1,
  81 quotient frames each with 27 lifts, a verified six-owner base, and six
  disjoint-support actual symmetries proving its lower bound.
- The three actual Hamiltonian paths in every chain of lengths one through
  eight,16,31,64; exhaustive chart assignments through eight cells; and the
  bounded-box root primitive-move control. The unbounded move-size statement
  uses the accompanying zero-sum decomposition proof, not extrapolation.
- Every literal prefix through substitution depth 12, plus exact selected
  addresses at depth 100, whose word has 927,372,692,193,078,999,176 slots.

An unpublished first implementation lacked the full-block return in the
prefix evaluator and tried to expand huge subwords. The author added the
necessary cached-summary return before final replay; the corrected checks
complete in under a second here. This repair matters for the claimed
compressed-address operation, though the concatenation identity itself
was unaffected.

The illustrative Fibonacci letter `A` acts as a transposition on three
states, so that particular finite-state example is not a native internal
automorphism of `C3`. It tests the separate noncommutative adapter; the
native triangle-chain chart and path controls are different declared
tests. Keeping this distinction prevents a representation test from being
misreported as an additional tournament automorphism.

## Additional actual hostile for the quotient and observation caveats

The later [non-normal kernel test](NON_NORMAL_KERNEL.md) independently
checks all `9!` permutations and all chart pairs for `diag(J3,C6)` with
circulant offsets `{0,1,3}`. It yields full group 81, block stabilizer nine,
projected domains and local kernels three, and one quotient frame with
nine lifts. The kernel is not normal in the full group. A literal target
transport shows that the source kernel must act on the RIGHT of a
comparison map; using it on the left fails. The owner image has three
values with multiplicity three each, proving that unweighted kernel
multiplicity is not a substitute for an actual-port observation.

Both normal and optimized replays of `audit_non_normal_kernel.py` returned
**PASS, 18 checks**. This is a concrete adversarial control for exactly
the projection-domain and observation-invariance conditions audited above.
