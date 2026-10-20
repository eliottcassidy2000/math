# The nonzero Walsh support is a robust matching relation

**Status: PROVED; FINITE-EXACT support controls through n=6.** This gives a canonical probability coupling and strong matching flexibility independently of any selected odd prime.

Put a_T=|Aut T|, b_X=|Aut X| and m=binomial(n,2). Define

    U[T,X]=sqrt(b_X/(2^m a_T)) Phi_X(T),
    B[T,X]=U[T,X]^2.

The two Walsh orthogonality relations make U orthogonal. Hence B is doubly stochastic: every row and every column sums to one. Its entries are nonnegative rational numbers and its support is exactly Phi_X(T)!=0. Unlike the signed matrix, B is independent of every graph-orientation gauge.

## 1. Every supported pair extends

The cofactor identity is especially transparent:

    (M^-1)[X,T]=b_X M[T,X]/(2^m a_T).

Thus an entry M[T,X] is nonzero if and only if its complementary determinant is nonzero. Expanding that complementary determinant proves that every nonzero entry lies in a complete support matching. Any prescribed supported pair can therefore be the first step of rational cofactor recursion.

In standard nonnegative-matrix terminology B has total support. The terminology and its classical context can be checked in [Sinkhorn and Knopp, *Concerning nonnegative matrices and doubly stochastic matrices*, Pacific J. Math. 21 (1967), 343–348](https://msp.org/pjm/1967/21-2/pjm-v21-n2-p14-s.pdf). No new general theorem about doubly stochastic matrices is claimed here.

## 2. Strict Hall expansion from the empty graph

For any nonempty proper set A of tournament rows, let N(A) be the even-graph columns with a nonzero coefficient in at least one of those rows. Then

    |N(A)| >= |A|+1.                                    (1)

Indeed the total B-mass of A is |A| and is contained in N(A), so |N(A)|>=|A|. If equality held, these columns' entire mass would come from A. But the empty-graph column has positive entry in every row, since Phi_empty=1. It belongs to N(A) and receives positive mass from rows outside A, a contradiction.

Thus deleting any one tournament row and any one graph column still leaves a support graph with a perfect matching: every subset of the remaining rows loses at most one neighbor from the strict bound (1). This also holds when the deleted pair's own coefficient is zero; in that case the remaining numerical determinant vanishes, although its support still admits a matching. This is a useful separation between structural support and cancellation.

The strict bound is optimal as a general statement: at n=4 the single odd-edge column (the star) is adjacent to only the two mutually converse non-self-converse tournament types. Deleting both those rows and two other columns can isolate it. Therefore arbitrary two-row/two-column deletion is already false at n=4.

The p-nonzero support can be much less flexible. For n=3,p=3 the matrix [[1,1],[1,-3]] has a nonzero upper-left entry modulo 3 that lies in no complete matching modulo 3. Its valuation filtration explains the failure; real total support does not descend to arbitrary residue supports.

## 3. An explicit distribution on support bijections

Starting with B, repeatedly choose the first support-perfect matching in a fixed canonical order. Subtract from its entries the smallest current entry on that matching. The residual has identical row and column sums, so while it is nonzero its support still satisfies Hall by the same mass argument. Each subtraction removes at least one positive entry; after at most E-N+1 steps, where E is the initial number of positive entries, the residual is zero.

The subtracted positive rational weights sum to one. They express B as a convex combination of support-permutation matrices and give an exact sampler of class bijections with marginals

    Pr(beta(T)=X)=b_X Phi_X(T)^2/(2^m a_T).               (2)

This is the elementary Birkhoff decomposition algorithm, stated here to make the probability construction explicit. The marginal B is canonical; the finite decomposition depends on the canonical tie order. A deterministic support bijection is obtained by taking any one matching in the decomposition or by the rational cofactor recursion.

For the empty graph, b_empty=n!, so (2) says

    Pr(beta^-1(empty)=T)=n!/(2^m |Aut T|)
                        = |orb(T)|/2^m.

Thus the tournament sent to the empty graph has exactly the type distribution of a uniformly random labelled tournament. This is a concrete probabilistic interpretation of the canonical coupling.

## 4. Exact finite scope

The complete atlas verifies stochasticity and, independently by augmenting paths, every one-row/one-column deletion through n=6. Nonzero entry counts are 4/4, 12/16, 122/144, 2329/3136 for n=3,4,5,6. The script does not sample the Birkhoff distribution; (2) is established by the proved decomposition, not a simulation.

The matching relation preserves nonzero signed-copy imbalance. It does not provide a local graph transformation, preserve all automorphism information, or eliminate the need for a complete class atlas. Rational matchings, modular matchings, and parity-block modular matchings are different valid choices with different guarantees.
