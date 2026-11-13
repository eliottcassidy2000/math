# Positive birth depth and the cost of compiling an outside model

**PROVED elementary corollaries of inherited mechanisms; VERIFIED finite
controls.** This note changes the requested construction question. It does
not add an all-order Hamiltonian-count realization or prove Chan–Pak.

## Inheritance and the comparison

The positive odd-cycle formula is [THM-002, ocf](../../../01-canon/theorems/THM-002-ocf.md).
The exact deletion kernel and core budget are inherited from
[Session133, core_reduction/THEORY](../session133_fractran_spectrum_refinement_20261101/core_reduction/THEORY.md),
which credits the October20 spectrum work and Busch. The corrected neutral
birth criterion is [THM-1900, the-insertion-response-calculus](../../../01-canon/theorems/THM-1900-the-insertion-response-calculus.md).
The [previous matrix packet](../fill_tree_matrix_20261113/matrix/PROOF.md)
already explains why a formal odd-cycle inventory need not be a tournament
inventory. The hard-core interpretation itself is not new.

Here a birth adds one vertex while retaining every old arc. A repair that
changes an old arc is a different operation. Counts are full directed
Hamiltonian vertex-word counts, H(empty)=1. Let c(T) be the number of vertices
in nonsingleton strongly connected components (SCCs). The inherited exact
budget is

    K(M) = max {3a+4b : a,b>=0, 3^a 5^b <= M},
    c(T) <= K(H(T)) <= 3 log_3 H(T).

The last bound is attained by ordinal sums of directed triangles. It bounds
cyclic core size, not total size: arbitrarily many singleton SCCs cost no H.

## G1. A pure-birth branch has only logarithmically many positive steps

For a nested birth chain T0,T1,...,Tm, let k be the number of indices i with
H(Ti)>H(Ti−1). Then

    c(T0)+k <= c(Tm) <= K(H(Tm)) <= 3 log_3 H(Tm).       (G1)

**Proof.** Every old cyclic vertex remains cyclic because its old cycle is
retained. By the inherited deletion-kernel theorem, a birth increases H
exactly when its newborn belongs to a nonsingleton SCC of the child. Each
positive step therefore supplies a distinct new cyclic vertex which survives
to Tm. Apply the core budget to Tm. QED.

This works with arbitrarily many neutral births between the positive ones,
arbitrary choices of old attachments, and unbounded order. If the successive
positive-step counts were H0+2j, or any polynomial in j, (G1) would fail for
large j. More generally their logarithms cannot be o(j). Thus an infinite
ordinary +2 chain of induced extensions is impossible even without a
unit-ready condition. This is a corollary of an existing core theorem, not
a new obstruction to Chan–Pak or to exponentially growing digit recursions.

**Failure boundary.** A branching family can still realize linearly many
different counts on different branches. A birth followed by rewiring can
destroy old cycles and is outside (G1). Neutral padding is not the missing
resource for a persistent additive chain; changing the old cyclic structure
or abandoning the chain is necessary. The one-shot absorption question in
[the absorption lane](absorption/PROOF.md) deliberately does not require the
output to be usable for another correction.

## G2. The outside hard-core model has the arithmetic but loses realizability

For an undirected graph G set Z(G)=sum over independent sets I of 2^|I|.
Adjoining an isolated vertex sends Z to 3Z. Adjoining a vertex adjacent to
every old vertex sends Z to Z+2. Starting from the empty graph Z=1, these
two operations realize every positive odd integer by threshold graphs.

**Proof.** Work backward from odd h>1. If 3 divides h, divide by 3; otherwise
subtract 2. Every value stays positive and odd and strictly decreases. Before
each division there are at most two subtractions. Reversing the resulting
word gives the graph. The number of vertices is at most
3 floor(log_3 h)+2, with h=1 using zero. QED.

This is an elementary model calculation, not an external novelty claim.
It is useful precisely because it exhibits an efficient arithmetic compiler
whose direct tournament analogue fails. Under OCF, H(T)=Z(Omega(T)), where
Omega(T) has a vertex for each directed odd cycle and an edge when two cycles
share a tournament vertex. The assignment T->Omega(T) preserves H but loses
the shared arc incidence needed to realize the cycle inventory.

There cannot be an all-odd scalar-preserving lift of this entire threshold
graph compiler: its outputs include 7 and 21, which tournaments exclude.
Even restricting the desired outputs does not justify lifting an isolated
or universal *conflict-graph vertex* to a single tournament birth. G1 rules
out an indefinitely nested sequence of universal-vertex +2 moves. The
previous packet's three-page example already forces an extra five-cycle
where an abstract conflict graph predicted 7. No new omission proof is
claimed here.

**Revised question.** Which changing, physically realizable cycle inventories
implement a whole arithmetic word, with an auditable full count? Do not ask
for each graph operation to have the same vertex-preserving tournament lift.
The block repairs in this session supply one limited example of changing
the inventory while retaining H; they do not create new prime counts.

## G3. An easy general-digraph witness can be too expensive to preserve

For n>=3, the directed cycle C_n has H(C_n)=n. Any tournament containing its n arcs
on the same n distinct vertices has those vertices in one SCC. The cited
strong floor, together with SCC multiplicativity, gives

    H(T) >= 5^((n−1)/3).                                (G3)

Extra vertices do not help: that SCC has at least n vertices. Thus no such
completion compiler has H(T)<=A n+B for all n with fixed finite A,B. For
n=5, the exact strong minimum is 9 although the target 5 is realized by
a four-vertex tournament. This is a small witness to the cost of insisting
on retaining the input's geometry.

The general-digraph reduction is nevertheless useful with the right output
predicate. The inherited [CT2 matching-response compiler](../digraph_matrix_boundary_20261110/completion/PROOFS.md)
already realizes every nonnegative integer as a selected-matching response
of a tournament, using logarithmic order for positive targets. It does not
realize it as the full H. The unresolved step is controlling the unwanted
paths, not obtaining an arbitrary number in an intermediate coefficient.

**Revised question.** Can a compressed, marked response construction be
compiled into a full-count construction with a prescribed, nonnegative
correction? A signed extraction identity by itself is not such a compiler.
This question must be attacked at the interface with the full path count,
not by enlarging the unary cycle completion.

## G4. A generous rewiring budget can make the wrong question automatic

Suppose two counts h,h'<=M are already known to be attainable. Delete all
singleton SCCs from chosen witnesses and pad both resulting cores to order
N=K(M) with source or sink singletons. The resulting tournaments retain h,h'.
Identify their vertex sets arbitrarily and reverse all disagreeing arcs in
one batch. This changes at most binomial(N,2)=O((log M)^2) arcs.

Thus if h'=h or h±2, an endpoint-neutral or endpoint-unit batch repair between
*some chosen witnesses* exists automatically. This deduction does not promise
small intermediate counts when the batch is performed one arc at a time.
Nor does it give an efficient construction without the target witness.

The important quantifier is existential choice of both witnesses. This says
nothing about every physical witness in an order-capped graph, preservation
of an existing boundary, or a uniformly constant arc radius. It explains why
a future claim of logarithmic-size batch rewiring, on its own, would not
advance exact range: it can merely restate that the endpoint already exists.
A useful repair theorem must supply the target witness from an explicitly
available source, certify a much narrower local operation, or preserve a
specified interface needed by a later construction.

## Verification and sources

`python -B 05-knowledge/results/strategy_fibres_absorption_20261113/growth_check.py`
checks every labelled tournament through order6, all prefix chains and their
starting points, independently checks H by permutations through order5, and
checks the threshold compiler on all odd h<=10001 and three large targets.
The finite universe and controls are in `growth_receipt.json`; the universal
arguments above do not rely on extrapolating those computations.

Primary sources re-opened November13,2026:

- Arthur H. Busch, [A Note on the Number of Hamiltonian Paths in Strong Tournaments](https://www.combinatorics.org/ojs/index.php/eljc/article/download/v13i1n3/pdf),
  Theorem1 and the opening SCC factorization: the inherited strong floor.
- John Irving and Mohamed Omar, [Revisiting the Rédei–Berge Symmetric Functions via Matrix Algebra](https://arxiv.org/html/2412.10572v3#S3.SS4),
  Corollary20, crediting Grinberg–Stanley: the inherited positive cycle formula.

No assertion about a hard-core phase transition is needed. The useful
outside comparison is an exact construction and the exact information its
transfer loses.
