# Decoding the unequal cyclic reversal cube

## Statement

Let A_1,A_2,A_3 be fixed nonempty named tournaments satisfying the six
cross-profile inequalities from ../unequal_covered_response/PROOF.md:

q_(k+2)(A_i)q_k(A_j) <= q_(k+1)(A_i)q_(k+1)(A_j)

for k>=1 and i!=j, with q_0=0 and entries past the block order zero.
Let r be the number of blocks having at least two vertices. Consider the
states obtained by choosing the orientation theta of the cyclic quotient
and, independently, whether to converse each nontrivial child interior.

**Corollary.** The 2^(1+r) resulting named full two-port responses are all
distinct. At most 3+r adaptive actual queries, or 3+2r fixed actual queries,
suffice to identify the state. In particular, for three nontrivial children
all sixteen states are distinct, with six adaptive or nine fixed queries.

These are upper bounds, not optimality claims. Each actual query uses the
same physical old names and the fixed arc x->y. A sum of three query values
is not counted as one query.

## Proof

For each theta choose a baseline in which the three children have their
given original orientations. Query the three block-uniform assignments
(a,c)=(1_A1,1_A2),(1_A2,1_A3),(1_A3,1_A1).

Every child is a module in each actual two-port extension: its incidences
to the other old blocks and both new vertices are uniform. Reversing that
child therefore preserves the Hamilton count exactly by the maximal-run
bijection. Consequently the vector of these three query values, and in
particular its sum, depends on theta alone, not on any child-converse bit.

The unequal cyclic theorem says that the sum for theta=0 is strictly less
than the sum for theta=1. Thus these three queries determine theta.

For that fixed theta and a nontrivial child A_i, use the existing single-
module isolation theorem: choose one physical two-port assignment that
separates A_i from its converse in the chosen baseline ambient tournament,
with both ports having A_i's exterior incidences. Its contrast is 2 modulo
4. This is the same lifted isolation used by the pairwise-disjoint family
theorem in ../simultaneous_response/PROOF.md; it applies to every nontrivial
module, including one of order two.

For j!=i, both ports see A_j uniformly because arcs between disjoint
modules are uniform. Hence A_j remains a module in this actual extension.
Conversing any of those other children preserves the query count exactly.
The chosen child query therefore has two distinct possible values depending
only on A_i's converse bit. It decodes that bit independently of every
other child bit. Repeat for each of the r nontrivial children.

Adaptive strategy: take the three root queries, determine theta, then use
its r corresponding child queries. This uses 3+r actual queries.

Fixed strategy: take the three root queries and all r child queries for
each of the two possible quotient orientations. After determining theta,
read the r appropriate child coordinates. This uses at most 3+2r queries;
coincident assignments can be merged but need not be. Thus full response
equality cannot occur between any two states.

## Operations, state count, and boundaries

The whole-parent reversal toggles the quotient orientation and all child
interior orientations. A child reversal toggles just its own interior.
For the nontrivial child bits this is an invertible change of coordinates
between operation bits and the state parameters. A singleton converse
does nothing and contributes no independent bit. Distinct states really
are distinct named tournaments: quotient orientation changes interblock
arcs, and each nontrivial converse changes at least one internal arc.

Every three blocks of orders two or three satisfy the cross inequalities,
so all sixteen responses are distinct in all these unequal as well as
equal choices. No quotient reconstruction is inferred for arbitrary block
profiles. This corollary concerns the indicated cyclic two-level family,
not an arbitrary strongly connected outer ambient substitution or an
arbitrary legal module-reversal sequence.

The incoming isolation theorem is used as a proved dependency. This packet
does not replace its mod-four proof with the finite controls below.

## Reproduction

verify.py independently constructs the sixteen states for mixed blocks
(TT2,TT2,C3) and (TT2,TT3,C3). It finds literal lifted child queries, verifies
their independence from every other child bit, and checks the nine fixed
query signatures are pairwise distinct. It also checks the singleton
boundary and counts only nontrivial child bits. The full mixed7 polynomial
contrast, including cubic and quartic data, remains in the separately
audited ../unequal_covered_response packet.

    python verify.py > normal.json
    python -O verify.py > optimized.json
    cmp normal.json optimized.json
