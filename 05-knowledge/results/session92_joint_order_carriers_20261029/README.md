# Joint order, compositional memory, and hidden return phases

2026-10-29. Source snapshot `b07ec51b204e5c0e09e52eed206cb2bd3cc29ebc`.
This session develops exact connections among tournament paths, weighted
observers, divisor fibres, LCA regrouping, signed operators, simplicial
complexes, Fano incidence, and LRC phase controls. Classical mechanisms are
credited. No external-priority or major-conjecture closure is claimed.

**Original witnesses recovered and verified.** The source arrived during this
session in commit `eea15f5d968c8e493a809862bbfec3989d09c436`, followed by the
selector-scope clarification `c72bf26e15e362bc3f53872698ff396369ea6d80`.
The supplied64→24→16 and named-vertex35 counts are independently reproduced.
All four Hamilton counts, 2,483,993 / 2,491,877 and2,899,675 / 2,906,091,
are now verified with two subset recurrences and a separate meet-in-the-middle
implementation. See [the exact recovered inputs and guards](runs/RECOVERED_RUN_WITNESS.md).
Earlier independent substitute controls remain explicitly separate.

## 1. A native quadratic query detects what linear continuation misses

The [folded observer census](observers/FOLDED_OBSERVER_COSTS.md) reconstructs
W=C3[C3,C3,C3], its anti-involution and all64 used-pair/current-vertex states.
Their96 transitions give33 folded Hamilton paths. The centralizer quotient
has24 states, and the exact minimal weighted partition has16. Adding
“current vertex is3” gives35.

Minimality is for partitions preserving level, final flag and weighted outgoing
transitions. It is not a lower bound against arbitrary linear algorithms.
The base linear continuation space has rank12; its quadratic products span
all16 partition functions. Its12 coordinates already separate the16 cells,
so this first rank gap concerns linear dependence, not pairwise blindness.

With the named3 observer, there is a genuine native hostile. The linear space
has rank21 and34 distinct signatures, but the exact partition has35 cells.
States(mask8,current7) and(mask4,current8) agree on every power of the transition
operator applied to every output-colour indicator. One chosen future count
has values(1,1,1) on the first state's successors and(0,0,3) on the second's.
Their sums agree; their squared sums are **3 and9**.

This is a shared-successor branching query, not two unrelated samples. The
claim concerns these linear continuation observables, not an unrestricted
Hankel-minimal realization of all interleaved observation words.

All256 menus of named current-vertex tests were checked. Exactly six minimal
six-name menus recover all64 states; every seven-name menu does. The observer
cost is neither submodular nor supermodular, so generic diminishing-return
sensor heuristics have no guarantee here.

## 2. The algebraic reason joint observations need another closure

[Joint-closure algebra](algebra/JOINT_CLOSURE_ALGEBRA.md) gives the precise
finite-state rule: the coarsest weighted equitable partition corresponds to
the smallest rational function space containing the observations and1, closed
under both transition pullbacks and pointwise multiplication.

Transition closure alone is a linear-response model. Products test joint
properties of the same state. Classical partition-algebra/bisimulation theory
is cited, and the proof is given directly.

A six-state example makes the composition failure explicit: either of two
binary target coordinates separately has an exact three-cell quotient. Their
common refinement has five cells but is not stable; joint closure needs six.
The missing query is the product of the two coordinate indicators. This is
why separately complete responses cannot generally be combined by simply
intersecting their partitions once.

## 3. Exact tournament seams and arbitrary-order dependence

The [recovered strong-child witness](runs/RECOVERED_RUN_WITNESS.md) uses
A=C3[C3,1,1] with fixed outer vertices3,4 and middle orders(0,1,2) and(0,2,1).
All individual rooted cursor responses agree, including on the converse sheet,
because the middle vertices lie in one C3 orbit. But no single automorphism
implements their transposition. Along the same13-vertex C3 quotient walk,
native complementation detects exactly the two relative-parity classes.
All six middle orders were replayed, for both the transitive-four control and
the all-three-children-strong construction, recovering every quoted count.
The one-bit claim concerns this fixed fibre, not every run-cover interface.

[The seam carrier](runs/SEAM_CARRIER.md) first protects the inherited ordinary
substitution theorem: for fixed children and quotient, reordering descriptive
runs cannot change the underlying substituted tournament. A changed Hamilton
count must belong to a path-dependent partner or another specified operation.

For native complement, which reverses every arc except the chosen Hamilton
path edges, the skew matrix identity is exact:

    B_P = -A + 2D_P.

For a fixed run family, all internal edges cancel between two choices, leaving
the joint directed seam relation. For fixed labelled T, the labelled partner
actually reconstructs the whole Hamilton path from its agreement edges.
This is not a minimum-information theorem for H alone.

An independent strong five-vertex seed has partner counts11 and15. A shared-start
padding construction preserves those counts at every larger order while
keeping the original tournaments strong. The padded partners have a universal
source, so this is not an all-partners-strong or equal-cursor extension theorem.
After recovery, both original quoted pairs also satisfy the shared-start
hypothesis: the same padding gives strong-source families preserving all four
counts at every larger order, with exactly these limitations.

A separate all-order result enlarges the operation menu to independent edge
restorations. The Hamilton-count multilinear polynomial has exact degree n−1
at odd n>=3 and n−2 at even n>=2. For odd n>=3 and even n>=4, uniform even/odd parity ensembles on a
suitable path have every proper marginal identical but different expected H.
Thus the needed correlation order is unbounded in this richer language.
The restoration ensembles need not satisfy native Hamilton-seam guards;
that distinction is retained throughout.

## 4. Arithmetic regrouping has an exact source-preserving repair

[Divisor clone associativity](arithmetic/CLONE_ASSOCIATIVITY.md) gives a concrete
positive answer to the LCA/tree-shape problem. Multiplication maps named divisor
tuples to a native divisor. Uniform native-fibre correction preserves sums,
but resetting that correction independently at every internal merge is not
associative.

For three copies of a prime, conditioned on total divisor p, two binary trees
produce source probabilities(1/4,1/4,1/2) and(1/2,1/4,1/4), instead of the
correct(1/3,1/3,1/3).

Retain the number C_v(d) of descendant tuples producing d. At an internal node,
choose a split a,b with weight C_left(a)C_right(b)/C_v(d). These factors telescope,
so every tree gives the same full named leaf law, preserving all joint guards.
Counts can be replaced by extension masses for a declared product source law.

The two-factor corrected weights are independent **iff the factors are coprime**.
For6×6, restoring native mass yields a flag angle9/49, not36/91. Correct total
mass does not restore the old product-spectrum proof. Nor can a conditional
mean replace a Boolean guard inside a later product: conditional variance is
the exact discrepancy.

## 5. Hidden boundary phases realize arbitrary face-count obstructions

[The signed-operator construction](geometry_operator/HOLONOMY_FACE_CARRIER.md)
starts with the same matching constraint on m modules, each with three feasible
states. It changes only signs of feasible/forbidden flip couplings. Every
absolute transition magnitude is unchanged, and the entire legal operator
PBP is literally unchanged.

Yet the two-step leakage Gram can change. Given any hypergraph of forbidden
occupied-module subsets, its kernel has one line for each subset containing
no forbidden hyperedge. Consequently every simplicial complex containing its
named vertices occurs as these surviving sector labels:

    nullity = number of faces,       rank = 3^m − number of faces.

The proof is gain-incidence propagation around exchange cubes. A nontrivial
square holonomy kills that sector's kernel; Boolean monomials ensure that a
higher-order twist cannot disappear through an unjustified cancellation.
These are spectator-dependent couplings, not the onsite weights of Session91.

A single full-support hyperedge reverses one Hermitian boundary coupling.
Every proper occupied-module sector agrees with the untwisted model, but the
full leakage rank increases by one. The matrix difference has rank two; a
rank gap of one must not be confused with a rank-one perturbation.

The [Fano control](geometry_operator/FANO_FACE_CONTROL.md) makes the connection
to earlier coding work explicit. Seven forbidden Fano triples give face counts
1,7,21,28,7, hence64 surviving kernel lines and rank2123. The seven largest
faces are exactly the nonzero simplex-code supports. The whole64-face family
is not linear and is not identified with the folded64-state carrier. The rank
follows from the proved theorem and checked sectors, not a dense2187-by-2187
elimination.

Face count is not homotopy type. At four modules, the full simplex and the
independent-set complex of P4 are both contractible, yet have16 and8 faces,
respectively, and therefore different leakage nullities. This is the precise
boundary when comparing the incoming `lucas_mirror_gap_spaces_20261028` gap
topology with this operator construction.

This realizes a broad family of high-order composition failures without
claiming a fast counting algorithm: the matrices are exponentially large.
No actual tournament-to-quantum-model or LRC-allocation equivalence is asserted.

## 6. Separate source controls: selectors and LRC phases

[Independent source replays](phase_controls/SOURCE_PHASE_CONTROLS.md) verify
the newer joint-port packet's smaller explicit examples. On the old seven-pair,
native complement followed by a selected exceptional sink gives307/293;
selected C3 substitution gives777/981. Direct permutation enumeration agrees
with subset Hamilton dynamic programming on all four.

For LRC, velocity rows{1,3,6} and{2,3,6} at threshold1/3 have identical complete
local length germs, but their isolated times are{4/9,5/9} and{2/9,7/9}.
Adding speed7 kills only the first. Exact rational interval calculations check
the whole stated local threshold intervals, not just sample points.

These are inherited independently replayed controls, not new LRC theorems and
not the supplied multi-million Hamilton counts. They demonstrate the same
methodological boundary with an actual arithmetic phase: a scalar response
can be exact and still fail under a later joint constraint.

## Verification, dependencies and remaining limits

[Independent audit](audit/AUDIT.md) uses separate folded signatures, Hamilton
permutations, integer kernel witnesses, modular rank lower bounds, divisor-tree
laws and exact algebra checks. The audit explicitly identifies its dependencies.
Producer programs use the Python standard library; two audit programs use
SymPy for independent exact rational elimination. Other audits use the standard
library, with one sharing an audit rank routine. All released scripts are run
normally and with `python3 -O`; active exceptions remain in both modes.

The session supplies precise closure criteria, infinite hostile constructions,
a native nonlinear observer witness, and an associative arithmetic repair.
The original large-count run witnesses are now independently certified in
their precise read-only-cursor/native-partner scope. Unrestricted tournament classification, LRC(14), positive Collatz convergence,
and other major conjectures remain open in the scopes considered here.
