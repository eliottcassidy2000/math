# Session 41 independent audit: shared owners and zero-slot adjacency

Status: PROVED elementary obstructions; FINITE-EXACT independent controls.
This report does not certify a complete tournament-wiring rule. No Git operation
was performed by this audit lane.

## 1. The zero-length slot is geometrically invisible but combinatorially visible

Consider labelled weights A=1, B=1, C=1, then 2,3,5,... . A legal support
contains no two adjacent slots. In the compressed order A,B,C,2,3,...,
A+B is forbidden. In the retained order A,0,B,C,2,3,... it is allowed.
The zero slot is never selected in these conventions: it is an adjacency spacer,
not an extra binary choice. Allowing it to be selected creates another model.
The first changed fibre is N=2: compressed supports are {2} and {A,C};
retained supports add {A,B}. Thus “zero contributes no length” cannot justify
silently deleting its incidence position.

The sequence 1,0,1,1,2,3,5,... genuinely satisfies the Fibonacci recurrence
from its first three entries. Its parities are odd, even, odd, odd, even, odd,
with the usual period three. The compressed sequence 1,1,1,2,3,... does not
satisfy the recurrence at its initial triple, although it does on its ordinary
Fibonacci tail. The three colored units must not be equated with standard
Zeckendorf uniqueness.

Let r(n) indicate use of the smallest unit in the standard greedy expansion
in weights 1,2,3,5,... . For positive n the compressed fibre has size
3-r(n-1). For n>=2 the retained fibre has size 4-r(n-1)-r(n-2), and its n=1
fibre has size three. These counts split by whether A and/or B are selected;
the remaining tail has its unique standard expansion. Independent literal
subset enumeration confirms both formulas for every 1<=n<=100.

The Fibonacci lane's self-expansion coloring needs a specified selection rule:
proper nonadjacent support alone does not select between A+2 and B+2 at total
three. A proved selected rule may be valuable even where the unselected fibre
is nonunique. Equality of its color prefix is not a uniqueness theorem for all
allowed decompositions.

## 2. Individual attachment feasibility is not joint feasibility

On F3 the two affine reflections h0(x)=-x and h1(x)=-x+1 each admit a fixed
owner, respectively 0 and 2. Their common fixed-owner set is empty. Hence a
solver that tests each closed-walk attachment separately produces a false
positive for a shared owner. Across all 36 ordered pairs of affine maps
x->epsilon*x+delta (epsilon=1 or -1), the common-fixed-set sizes have census:
26 empty, 9 singleton, and 1 of size three. The checker enumerates this directly.

Likewise, put translations with shifts 0,0,1 on an oriented triangle. Every
edge relation individually admits all three choices for its source owner.
No three vertex gauges simultaneously flatten all edges, because their
oriented shift sum is 1, not 0. Any claim about local attachment classes must
retain one actual gauge at each shared vertex and impose the joint cycle law.
These are connection-model hostiles, not claimed actual tournament pairs.

## 3. The correct general compatibility carrier

For finite domains D_v and actual relations R_e on the participating named
ports, compatible assignments are the relational join of all R_e, with repeated
port names identified, not the product of their separate nonemptiness tests.
For bijections on a connected graph, choose a spanning tree and transport one
root owner to each vertex. Chord constraints then demand that this SAME root
owner be fixed by EVERY resulting based holonomy. This proves necessity and
sufficiency for a parallel section. It does not classify the connection under
gauge, demand identity of all holonomies, or authenticate native tournament
path reversal. Those are distinct predicates.

For an abelian translation connection, joint cycle shifts vanish exactly when
a vertex potential exists. For affine nonabelian transport, a scalar residual
potential cannot replace the based holonomy action. For capacity-constrained
linear transport, residual-cycle criteria apply to that declared feasible
polytope; added integer or shared-owner constraints must be retained in the
feasible set before concluding realizability.

## 4. Provenance and inheritance guards

The corrected summand record is
`05-knowledge/results/arithmetic_braids_20260917_summand.md`, rather than an
uncorrected use of `07-reflections/metagraph-summand-recursive.md`. It retains
paired-parent generation and corrects Fibonacci-depth claims. A graph shadow
may retain lengths yet lose the companion; analogously, deleting a zero slot
retains sums yet changes admissible supports. This is an exact loss-of-data
comparison, not an isomorphism between the arithmetic and tournament problems.

Inherited Session 39 affine-connection results distinguish common sections,
centralizers, and full frame trivialization. Session 40 complete marked-target
results do not determine a literal partner from a source and endpoints alone.
All subsequent audits in this directory retain those boundaries.

## Reproduction

Run `python check_independent.py` and `python -O check_independent.py`.
All checks use explicit exceptions, not assertions. The receipts record the
first adjacency discrepancy, full affine-pair census, triangle obstruction,
and the finite fibre verification range. No external package is required.

## 5. Reviewed packets and secondary independent controls

The Fibonacci packet `RESULTS.md` has been proof-audited for its complete
fibre comparison, proper-expansion count, lexicographic convention, coherent
row induction, and substitution frequencies. An explicit n>=4 qualifier in
the concatenation prefix argument was requested and incorporated; N=1,2 are
handled separately. Its canonical versus older row-35 discrepancy is retained.
The word formulas |D(N)|=3-1_A(N), |R(N)|=3+1_B(N-1) (N>=2 for the latter)
follow from every C being followed by A and every A being preceded by C.
Their means are phi^2 and 2phi, respectively.

The residual packet's bounded conformal-atom theorem has been proof-audited.
The basepoint-only criterion is sufficient: every conformal atom in a
decomposition of y-x individually lies within the coordinate interval between
x and y, so each can be applied already at x. It is not necessary to test
only some algorithmically discovered list without coverage proof.
Independent enumeration verifies that the two binary 2x2 tables constrained
by X00=Y00 have exactly two feasible states and exactly two nonzero bounded
kernel vectors, negatives of each other and each with l1 norm eight.
The same checker confirms rotational centralizers 3 and 1 for holonomy
subgroups C3 and S3, while neither admits a parallel owner.

One proof gate was sent to the residual author: reference internal cell
isomorphisms in a marked comparison must themselves respect retained marks,
or the allowed chart set must be the corresponding transporter coset rather
than a mark-fixing subgroup after an arbitrary reference. This is a genuine
coordinate boundary, not a problem with the conformal theorem.

The residual reference-map correction has now been incorporated by its author.
Both lane checkers were replayed independently in normal and optimized Python:
outputs are identical and PASS. The residual replay covers 33,750 basepoint
criteria over 250 matrices and all 4,096 binary triangle relation systems.
It finds 54 simultaneous-conjugacy false positives if separate root conjugators
are allowed. Replay receipts are included here. These finite checks complement,
rather than substitute for, the all-size proofs above.

