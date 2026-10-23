# The static rooted-motif boundary extends to every large order

**PROVED, scoped observer lower bound.** This separate addendum preserves all earlier frozen packets. The two-vertex operation is the r=1 case of the inherited cyclic lift in `05-knowledge/results/session38_incidence_ports_electrical_20261023/theory/INCIDENCE_COLLISION_PROPAGATION.md`. The new transfer statement here concerns ALL singly-rooted induced motifs and preserves self-converse, intrinsically recoverable examples.

## 1. The two-vertex lift is strong, self-converse, and injective on strong types

For a strong tournament T of order N>=3, define

    L(T) = C3[T,{a},{b}], with T -> a -> b -> T.

The lifted tournament is strong. If f is an anti-isomorphism of T, extend it by f(a)=b and f(b)=a. This reverses every arc of L(T), so self-converse is preserved.

Every old vertex has its degree increased by one. Because T is strong, its old degrees lie in [1,N-2], hence the lifted old degrees lie in [2,N-1]. The new degrees are

    d+(a)=1,    d+(b)=N.

Thus a and b are respectively the unique minimum- and maximum-degree vertices. Deleting these two intrinsically identified vertices recovers T. Consequently L(T) and L(U) can be isomorphic only if T and U are isomorphic. Repeated lifting preserves nonisomorphism by repeated deletion of the unique extremes. No Hamiltonian-count difference is needed for this injectivity argument.

## 2. Every static singly-rooted profile transfers

Suppose tournaments T and U have a fixed vertex correspondence under which the complete static singly-rooted induced-motif profiles through size k agree pointwise. Then the correspondence, extended by a_T<->a_U and b_T<->b_U, gives equal profiles through size k in L(T) and L(U).

For an old root v, divide subsets containing v by the subset I of {a,b} that they contain. Their old part is a rooted induced tournament on t-|I| vertices. Its rooted type determines the extended rooted type because every old vertex meets a and b uniformly and the arc a->b is fixed. Summing these four deterministic extension maps gives the new profile. The old rooted sizes required are at most k.

For a new root a or b, divide subsets by whether the other new vertex is present. Their old part is an UNROOTED induced tournament on t-1 or t-2 vertices, and again its type determines the extended rooted type. The empty old part is one fixed object. For every nonempty old unrooted type F on j vertices, its count is recovered from the old rooted vectors by

    count(F,T) = (1/j) sum_v sum_{tau: forget_root(tau)=F} count(tau,T,v).

Each induced copy contributes once for each of its j possible actual roots. Thus equality of the pointwise rooted profiles implies every old unrooted count needed for the new-root extension. This proves the transfer, including the new roots, without assuming vertex transitivity.

If a phase-blind equitable partition of T makes the static profile constant on each cell, it extends by the two identified singleton cells {a},{b}. All new cell-to-cell degrees are fixed by the displayed block directions. Its equitable quotient remains phase-independent, so ordinary unindividualized vertex refinement remains blind at every depth just as before.

## 3. Linear necessary cutoff at every order n>=11

ROOTED_CUTOFF_BOUNDARY.md constructs, for every g>=3, a pair of distinct strong self-converse tournaments T_(g,+1),T_(g,-1) of order 3g whose full static singly-rooted profiles through g agree at every corresponding vertex and whose ordinary unindividualized refinements agree forever.

For any n>=11, choose r in {0,1,2} so that n-2r is divisible by three, and put

    g=(n-2r)/3.

Explicitly r=0,2,1 for n congruent to 0,1,2 mod3 respectively. Then g>=3 and n=3g+2r. Applying L exactly r times to the pair gives distinct strong self-converse tournaments of order n with the same static rooted profiles through cutoff g and the same complete ordinary-refinement sequence.

Therefore a class separator on ALL strong self-converse tournaments of order n using only full static singly-rooted induced profiles through k(n), followed by ordinary directed vertex refinement, necessarily satisfies

    k(n) >= (n-2r)/3 + 1,     n>=11.                       (1)

In particular k(n)>=(n-1)/3 at every such order, and

    liminf_{n -> infinity} k(n)/n >= 1/3.                   (2)

These equalities are HISTOGRAM equalities, not indexed hypomorphy. Already the same three actual vertices {(0,0),(0,1),(g-1,0)} induce a transitive triple in T_(g,+1) and a cyclic triple in T_(g,-1). The root-fixed local bijections may move the other selected vertices and depend on the occupied cells; they do not preserve every indexed subset under one common vertex correspondence. No reconstruction theorem for indexed small-subset data is being challenged, and no external priority claim is made.

These are necessary conditions for the stated observer, not sufficient cutoffs or algorithmic complexity lower bounds. They do not restrict actual individualization, multi-root overlap data, path-decorated flags, or complete adjacency. The incoming prime-cyclic sqrt(log p) result has a different subclass and is reconciled in ROOTED_CUTOFF_BOUNDARY.md; (2) makes no claim about prime-cyclic examples.

## 4. A related exact-E consequence, using the inherited lift theorem

For g>=5, GIRTH_PHASE_HOSTILES.md gives three distinct strong self-converse T_(g,delta), delta=0,1,2, with equal exact E. The inherited Session38 theorem proves that L preserves exact-E equality. Sections1 and3 above preserve self-converse and recover the original phase class intrinsically.

Every n>=17 has a representation n=3g+2r with g>=5 and r in {0,1,2}. Hence there are THREE distinct strong self-converse tournaments with the same exact E at every order n>=17. In particular they remain separate converse-merged classes. This is a class-signature and global-forest obstruction, not a proof that any actual endpoint-constrained port ledger has alternative wiring.

## 5. Finite controls

`check_lifted_cutoff.py` imports only the separately frozen standard-library rooted canonicalizer `check_rooted_cutoff.py`, whose exact source hash is recorded in its receipt. For g=3,4,5 and r=0,1,2, it constructs both phase tournaments, enumerates the full rooted profile through cutoff g at every actual vertex, and checks equality and stable ordinary refinement. It checks literal anti-isomorphisms, strong connectivity, unique degree extremes, and recovery of the previous tournament at every lift. These controls cover orders9,11,13;12,14,16;15,17,19. A separate arithmetic check covers the representation and cutoff inequality for every n from11 through100; the all-order assertion is the proved residue calculation, not an extrapolation.

    python marked-cells/check_lifted_cutoff.py
    python -O marked-cells/check_lifted_cutoff.py

Both modes reproduce `lifted_cutoff_receipt.json`. These finite checks verify the new rooted transfer directly; the exact-E consequence cites the already-proved inherited lift instead of silently treating a rooted-motif check as an E check.
