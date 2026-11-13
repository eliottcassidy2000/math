# Tight suffixes give marked generating functions and blockwise inverse mechanics

## Exact sources

- [Acyclic upper-tail two-color exchange](../graceful_upper_tail_two_color_connectivity_20261019/UPPER_TAIL_TWO_COLOR_CONNECTIVITY.md): graphic contraction, rank staircase, closure-preserving double moves, and a connected shortest-color buffer.
- [Polynomial diameter](../graceful_polynomial_diameter_20261019/POLYNOMIAL_STAIRCASE_DIAMETER.md): all tight sets in a rank staircase are suffixes, and consecutive tight suffixes form strict-rank blocks after deleting one full-rank buffer. This is stronger than the older exponential diameter.
- [Frozen-fiber obstructions](../graceful_fiber_obstructions_20261019.md): a single-color fiber can be disconnected, and a fixed band of low colors is not a universal single-exchange buffer.
- [Metagraph Hall/Bott](../metagraph_hall_bott_20261023/README.md): inverse-free extraction of a marked harmonic phase and linking unit, not a Bott recurrence for metagraphs.
- [Collatz marked Bott](../collatz_bott_marked_periodicity_20261007d.md): a marked tensor corner recovers the original operation; forgetting to a stable or even similarity class loses native divisibility guards. Its simple example is that valuation words 12 and 21 have similar matrices, but their actual odd-source cells are 11 and 9 modulo 16.

The supplied order-five polynomial is from [Pairing Wiring](../joint_metagraph_dyadic_20261113/PAIRING_WIRING.md). That source establishes the actual independent flips for its merged-wiring model; the inverse test below applies to other candidate binary models. See [Voltage lifts](VOLTAGE_LIFTS.md) for what happens when the missing wiring is restored.

## Structural theorem: fixed-flat factorization, with explicit inverse

Let a finite matroid have disjoint classes \(A_1,\ldots,A_s\), satisfying the rank staircase \(\rho(A_i)\ge s-i+1\). A state is a choice of one actual element of every class, with all chosen elements independent. For a color set \(K\), put \(A_K=\bigcup_{i\in K}A_i\). The existing rank-staircase proof gives \(\rho(A_K)\ge |K|\); equality implies that \(K\) is a suffix, by applying the bound to \(\min K\). List the tight suffixes as
\[
\varnothing=K_0\subsetneq K_1\subsetneq\cdots\subsetneq K_t,
\]
and write \(F_j=\operatorname{cl}(A_{K_j})\), \(B_j=K_j\setminus K_{j-1}\). The residual prefix is \(B_*=[s]\setminus K_t\).

For \(j\ge1\), define \(X_j\) as the independent transversals of the classes \(B_j\) in \(M/F_{j-1}\), with original elements marked. Define \(X_*\) the same way on \(B_*\) in \(M/F_t\). Then restriction of a state to the blocks gives an exact bijection
\[
X\ \longrightarrow\ X_1\times\cdots\times X_t\times X_* .
\]

**Proof.** In any state the choices on \(K_j\) have cardinality \(|K_j|=\rho(F_j)\) and lie in \(F_j\), so they span exactly \(F_j\). Thus each restricted block is independent in the stated contraction. Conversely take one allowed choice per block. The first choice spans \(F_1\). If the choices accumulated through \(K_{j-1}\) span \(F_{j-1}\), the next allowed block extends them independently and makes the union a basis of \(F_j\), since the block size is the difference of the two flat ranks. Finally the residual block is independent after contracting \(F_t\). Taking the union of the **same original marked elements** is the inverse. In particular no parallel contacts are identified.

Legal changes in separate factors are independent: a change in one factor stays legal regardless of all other choices. An exchange changing one coordinate in each of two distinct factors is equivalent to two legal single-coordinate moves in either order; it provides a commuting square, not a new obstruction or rescue. A coupled two-color move is genuinely necessary only within a single block. This sharpens the distinction between independence of exact families and graph connectedness.

### Weighted and equivariant form

Let the weight of a state be an additive sum of integer weights on its chosen actual elements, and let \(g\) be a matroid automorphism carrying the family of classes to itself and preserving weights. It preserves rank and thus permutes the tight color sets. They are nested and have distinct cardinalities, so each is setwise fixed. Therefore \(g\) acts on each factor separately, preserving the full ordered factorization. For either the identity or any such \(g\), let
\[
P_{X,g}(x)=\sum_{z\in X:gz=z}x^{w(z)}.
\]
Then
\[
P_{X,g}(x)=\prod_j P_{X_j,g}(x)\,P_{X_*,g}(x).
\]
This statement is about fixed individual states, not fixed orbits. A constant baseline is simply an overall power of \(x\). The proof is the product bijection: a tuple is fixed iff all coordinates are fixed.

For an involution, \(x^8(1+x^2)^2\) is therefore a plausible trace of two independently variable fixed binary blocks of relative grades 0 and 2, with baseline 8. A checkable certificate is **four marked states** indexed by the complete binary rectangle, with the two block changes commuting, the two grade increments equal to 2, and no other fixed choices. Factorization of the polynomial alone does not supply that certificate. For example, take the four vertices of a path with weights 8,10,10,12 and attach two same-weight leaves at one endpoint; an involution exchanges the leaves and fixes the path pointwise. Its fixed-weight polynomial is identical, yet its fixed graph is a path with no commuting square. Even for a genuine product, the fixed-weight image is \(\{8,10,12\}\), not a consecutive interval. In particular connectedness and even grade parity are different conclusions.

### How this can transfer, and the precise obstruction

For a tournament or metagraph family the analogue must first provide:

1. a native legal-state space, together with its original interface markings;
2. a certified nested chain of boundaries at which *every* legal choice spans the same closure/interface;
3. a recombination law: independent choices on opposite sides of those boundaries glue to a legal native state;
4. if using the fixed-weight polynomial, both equivariance and additivity of the actual grading under the glue.

Items 2–3 are the matroid rank proof's content, not consequences of observing equal counts. Items 3–4 cannot be imported from Bott periodicity: marked tensor stabilization is a faithful bookkeeping map, not an existence theorem for the recombined native states. This narrows the strongest test on any proposed tournament family to an exact **rectangle completion** with all readiness sidecars, rather than matching numerical factorization or Smith summaries. If rectangle closure holds, the above product and commuting exchanges follow without a global enumeration. If a mixed corner fails, the first failed boundary coordinate is the appropriate obstruction to refine.

## A second inverse: signed cube involutions

For a Hamming-weight-graded binary cube, possibly with fixed coordinates or a fixed weight baseline, suppose the involution is a coordinate permutation followed by bit flips. Decompose the coordinate permutation into fixed points and transpositions. On a fixed free coordinate, the local fixed-point polynomial is \(1+x\) if unflipped and zero if flipped. On a free transposition the involution condition forces equal flip flags at its two ends. The local polynomial is \(1+x^2\) for two unflipped ends (solutions \(00,11\)), or \(2x\) for two flipped ends (solutions \(01,10\)). Fixed prescribed coordinates only supply their weight monomial if consistent; inconsistent prescribed coordinates give no fixed states. A one-free/one-prescribed transposition is impossible for an action preserving the whole subcube. These orbit checks give the exact product
\[
x^b(1+x)^f(1+x^2)^p(2x)^q
\]
when the fixed set is nonempty, where \(b\) includes the consistent prescribed baseline. In this particular **certified cube** model, the polynomial \(x^8(1+x^2)^2\) forces \(q=0\) by its least-degree coefficient, \(f=0\) by evaluation at \(-1\), \(b=8\) by the least degree, and then \(p=2\) by the degree span. Thus the factorization is reversible *once the cube and action are proved*, and then gives two explicit independent degree-two flips. The weighted-path hostile above shows why the same inference cannot be made for a general invariant state graph.
