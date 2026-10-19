# Gold partitions: terminal-chain contexts, flat equality defects, and endpoint-locality

Date: 2026-10-19. Status: **PROVED scoped** (arguments below); **FINITE-EXACT** controls shipped separately. The Gold Partition and general poset 1/3–2/3 conjectures remain **OPEN**. This is a new proof attempt and exact transport result, not a rerun of a conjecture census. No claim is made that the example families escape previously known Gold classes.

## 1. Precise target, inheritance, and what changed

For a finite nonchain poset P, write t_0=e(P). The Gold Partition Conjecture asks for two comparisons such that every realized outcome path has

    t_0 >= t_1 + t_2.

Here t_i counts the extensions after i comparisons. Comparisons may already be determined; if the first leaves a chain, t_2=1. The original weak inequality is essential. The two-element antichain has 2=1+1. See Peczarski, *Comments on the Golden Partition Conjecture* (2017), pp. 106–109:
https://cdm.ucalgary.ca/article/download/62349/46785/178291

Two comparisons fixed in advance suffice. After the first split H>=L, the light branch automatically satisfies L+t_2<=2L<=H+L, so only the heavy branch dictates the second comparison. Thus adaptive and fixed second-comparison formulations are equivalent here. The strict sign in the displayed conjecture in arXiv:2410.12494v1 is not the original conjecture.

Inherited results used, rather than reverified:

- `05-knowledge/results/gold_heavy_child_20261014j.md`: exact heavy-child defect H-L-B(heavy), signed condition 3a+|b+r|<=1, and the 14-vertex failure of maximally-balanced and conditional-best-response selectors.
- `05-knowledge/results/gold_partition_context_transport_20261014j.md`: weighted four-cell transport, strict-slack upper-context buffers, and a first zero-slack hostile.
- `05-knowledge/results/golden_order_bridge_20261014j.md`: Fibonacci thin-order extensions, all-length boundary certificate, and the cyclic seam that produces the eleven-state non-poset law.
- Peczarski's low-slave and balanced-triple sufficient certificates, cited in the 2017 source. The July 2026 order-14 computation is accepted as working inheritance, not recounted: https://arxiv.org/html/2607.23926v1 .

The new question was whether equality certificates can be transported using sufficiently accurate asymptotic uniformity, and whether a successful pair of disjoint comparisons can be replaced locally by overlapping ones. The answers are sharper than merely saying “keep the weights”: every positive-degree coefficient can agree while the same certificate fails forever; and even a genuine five-vertex poset can require a fifth vertex for endpoint-local uncrossing.

The active concepts were four-cell budgets, exact extension fibres, equality versus strict slack, finite-state Fibonacci boundaries, and hidden contextual comparisons. The positive result is an exact context classifier, with explicit strictification and certificate relocation on an infinite family.

## 2. Terminal-chain fibre theorem

Let Q be a finite poset with n vertices and D an order ideal of Q. For integers r>=1 and m>=0, form P_{m,r}(Q,D) as follows:

1. Add z_1<...<z_r, with precisely the additional core relations d<z_r for d in D.
2. Add an m-element chain above every vertex of Q and incomparable with the z-chain.
3. Add no other cross-relations, apart from transitive consequences.

For sigma in L(Q), let d_sigma be the position of the last D-element, with d_sigma=0 if D is empty. Put

    T = binom(m+n+r,r),
    h_r(d) = binom(r+d-1,r), with h_r(0)=0.

**Theorem 1.** The exact number of extensions of P_{m,r}(Q,D) restricting to sigma is

    w_{m,r}(sigma) = T - h_r(d_sigma).                         (1)

**Proof.** Fixing sigma orders Q followed by the m-chain into one backbone of length n+m. Without the attachment to z_r, its shuffles with the r-chain number T. An illegal shuffle has z_r before the last member of D. Since z_r is the final z-element, all r z-elements then occur before backbone position d_sigma. They can be interleaved with its first d_sigma-1 vertices in exactly binom(r+d_sigma-1,r) ways; the remaining backbone is forced. These are all illegal shuffles, and no others. Subtract them. For empty D none is illegal. ∎

This formula is not a module argument: the terminal z_r can distinguish arbitrarily many different ideals D inside Q. The fibre defect is independent of m.

### Exact classifier for a fixed weak certificate

Fix two comparisons on Q, and write C_{ij} for their four extension cells, N_{ij}=|C_{ij}| and N_i=N_{i0}+N_{i1}. Put

    Delta_{ij} = N_{1-i}-N_{ij},
    S_{ij}(r) = sum_{sigma in C_{ij}} h_r(d_sigma),
    K_{ij}(r) = S_{ij}(r)-S_{1-i,0}(r)-S_{1-i,1}(r).

Restriction preserves both answers, so the exact path deficit in the enlarged source is

    T Delta_{ij} + K_{ij}(r).                                (2)

Suppose the source certificate is Gold, so all Delta_{ij}>=0. Then for fixed r:

- If Delta_{ij}=0 and K_{ij}(r)<0 for any cell, the same certificate fails for every m.
- Otherwise it succeeds precisely when T>=H_r, where

    H_r = max(0, max_{Delta_{ij}>0} ceil(-K_{ij}(r)/Delta_{ij})). (3)

Since T increases strictly with m, this is an exact terminal interval m>=m_*, not only a sufficient eventual bound. A monotone integer search finds m_*. This is a complete all-m classifier for the supplied marked certificate. It does not optimize over alternative comparisons.

### Why arbitrarily accurate uniformity still misses equality

For fixed n,r, all fibre polynomials in m have the same coefficients of m^r,m^{r-1},...,m. Only their constants vary. Their normalized restriction law approaches the uniform law on L(Q) at order O(m^{-r}). Therefore matching an arbitrarily long finite asymptotic expansion is insufficient to transport an equality certificate. Section 3 gives an explicit failed budget for every m, with a constant positive integer defect.

This strengthens the inherited r=1 hostile. It does not refute Gold on the enlarged poset, and it does not say that all near-uniform contexts fail. It says that a weak inequality cannot be inferred from a limiting equality without the signed constant retained in (2).

## 3. A two-boundary Fibonacci family

Let Q_N, N>=3, have vertices 0,...,N-1, with i<j exactly when j>=i+2. Write F_0=0,F_1=1. Its extensions are identity orders with disjoint adjacent swaps. Equivalently they are independent sets of the path of N-1 swap positions. This inherited bijection gives the left boundary certificate

    compare (0,1), then (1,2),
    core table [[a,b],[a,0]],
    a=F_{N-1}, b=F_{N-2}.                                   (4)

Rows and columns below put the natural increasing answer first.

### 3.1 Attachment above vertex 0: wrong boundary, exact defect, repaired boundary

Take D={0} in Theorem 1 and put T=binom(m+N+r,r). Vertex 0 has position 1 unless the first swap is occupied, when its position is 2. Hence

    w(sigma)=T-1-r*1{first swap occupied}.

The left table becomes

    [[a(T-1), b(T-1)], [a(T-r-1), 0]].                        (5)

The first row is heavy and its largest leaf has size a(T-1). Its Gold-budget violation is exactly

    t_1+t_2-t_0 = r F_{N-1} > 0.                             (6)

It fails for every m, however nearly uniform the law becomes. In fact

    max w / min w = 1 + r/(T-r-1) = 1+O(m^{-r}).             (7)

But the source itself has an explicit Gold certificate at the other boundary:

    compare (N-2,N-1), then (N-3,N-2).

Its table is [[A,B],[A,0]], where, using F_{-1}=1 only for the N=3 endpoint convention,

    A=(T-1)F_{N-1}-r F_{N-3},
    B=(T-1)F_{N-2}-r F_{N-4}.                               (8)

**Proof of the right certificate.** In the cells where the penultimate swap is absent, toggling the last swap is a bijection. It preserves the first-swap flag and hence the fibre weight. These are the two A cells. The B cell injects into an A cell by deleting the occupied penultimate swap, without changing the first-swap flag for N>=4; its weighted contribution is no larger. For N=3 the explicit table is [[T-1,T-r-1],[T-1,0]], so B<=A directly. Thus 0<B<=A and the heavy row A+B has largest leaf A. The total is 2A+B, giving exact equality on that leaf. Formula (8) is the usual independent-set count with and without the first position occupied. ∎

At N=4 the B-to-A deletion can involve the first swap only if one uses the wrong boundary indexing; the penultimate swap is position 1, and the first flag is position 0, so it remains unchanged. The direct N=3 exception is necessary.

The selector lesson is constructive: locate the perturbation and move the equality injection to a boundary that preserves its weight. A global count or a near-uniform approximation does not identify that boundary.

### 3.2 Attachment above vertex 1: strictification

Now take D={1} and put q=binom(r+2,2). The position of vertex 1 in the three nonempty left cells is respectively 2,3,1. Thus the exact left table is

    [[a(T-r-1), b(T-q)], [a(T-1), 0]].                        (9)

**Theorem 2.** Table (9) is Gold for every N>=3,r>=1,m>=0. Every realized path inequality is strict except when (N,r,m)=(3,1,0) or (4,1,0).

**Proof.** Denote its three entries A,B,C in row-major order. We need A<=C, B<=C and C<=A+B. The first follows from C-A=ar>0. The second follows from b<=a and q>=1. The third is B>=ar.

For N=3, a=b=1 and

    T-q >= binom(r+3,r)-binom(r+2,r)=binom(r+2,3)>=r.

For N>=4, a<=2b, and

    T-q >= binom(r+4,r)-binom(r+2,r)
          = r(r+1)(r+2)(r+7)/24 >= 2r.

These prove B>=ar. Equality in the required bound occurs only at N=3 or 4, r=1,m=0. In all other cases B>ar, and the other actual path deficits are strictly positive as well. ∎

This gives a constructive passage from an equality certificate to strict slack, without a uniform-fibre assumption. The inherited strict-slack upper-context theorem can then attach arbitrary upper context after a sufficiently long further buffer.

**Prior-art boundary.** Peczarski's low-slave theorem already covers these Gold examples: the selected boundary has at most one slave and the displayed table proves its required majority orientation. The contribution here is the exact nonuniform kernel, all-parameter strictification, quantitative failure at the wrong boundary, and arbitrarily high-order flatness. It is not a claim to have settled a previously open general poset class.

## 4. A genuine endpoint-locality obstruction

Take a four-vertex core with a,b below c,d. Add z with a<z<c. Restricting the eight extensions to the four core vertices gives the following fibre table, by the orders a/b and c/d:

    [[2,3],[1,2]].                                          (10)

Comparing a with b first and c with d second is Gold: its heavy branch has size 5 and its largest leaf size 3, giving 8=5+3.

However no Gold certificate using overlapping pairs from {a,b,c,d} exists. The only uncertain core comparisons are {a,b} and {c,d}; these are disjoint and neither is half-balanced. A determined first comparison cannot work, and a determined or repeated second comparison leaves the heavy branch of size 5 unchanged, violating 8>=10. Thus even for a genuine poset a local uncrossing rule may need a fifth vertex.

There is an exact unbounded continuation. Replace z by an ell-chain between a and c. The table becomes

    [[ell+1, binom(ell+2,2)], [1,ell+1]].                    (11)

This follows by distributing the ordered interior chain among the slots before c in each fixed core extension. The two disjoint core comparisons, in either order, satisfy Gold exactly when ell<=1. For ell>=2, no two-comparison strategy using only the four core endpoints works at all. The full poset has width two, so its Gold property is already known; the failure is endpoint locality, not the conjecture.

This exhibits three distinct levels of information:

- the induced core order, unchanged with ell;
- the exact four-cell weights, changing quadratically with ell;
- the extra comparisons involving an interior z_j, needed once ell>=2.

A core-only comparison search can fail even when it retains the exact core law. Merely restoring weights does not automatically restore a sufficient comparison alphabet.

## 5. A useful new proof target, and its limits

Theorems 1–2 explain one way to transport equality: preserve the injection pointwise, or change the comparisons so that the contextual weight becomes invariant. They also explain why matching finitely many asymptotic coefficients is not an adequate substitute.

A sharply stated remaining target is this:

> Given a uniform linear-extension law and a proposed four-cell certificate on an induced core, find a bounded boundary enlargement that either transports its nonpositive budget defect or supplies a different two-comparison certificate.

No universal bound is claimed. Section 4 shows that “no enlargement” is false. Section 3 shows that “use the same comparisons once the law is close enough to uniform” is false. The terminal-chain formula gives a tractable testing ground in which both the exact loss and an alternative certificate can be calculated.

For the requested golden-lane connection: the same Fibonacci transfer matrix controls the swap path, but the proof uses open boundary conditions and full extension fibres. Closing the path into the eleven-state cyclic no-11 law still deletes extensions without a poset-realizable closure, as the inherited golden-order bridge proves. The present results retain that seam warning and add a new boundary selector; they do not turn a periodic Collatz parity law into a uniform poset law or prove a Collatz statement. No square-sum/graceful implication is asserted here.

## 6. Reproduction and declared universe

Files:

- `gold_terminal_chain.py`: exact ideal-DP counts of the actual enlarged posets, rather than merely checking the symbolic formulas against themselves.
- `gold_terminal_chain.json`: the finite receipt in `05-knowledge/results`; the script prints a fresh receipt.

Run:

    python3 05-knowledge/experiments/gold_terminal_chain.py
    python3 -O 05-knowledge/experiments/gold_terminal_chain.py

The script checks every N=3..10, r=1..7, m=0..5 for both marked attachments, all left tables, the negative family's repaired right tables, exact failure defects, and strictness boundaries. It separately verifies the fibre theorem over every naturally labelled poset on 2..5 vertices, every nonempty ideal, every core linear extension, and (r,m) in {(1,0),(2,1),(4,2)}. The bridge family is checked for ell=0..14, including every ordered pair of core comparisons. This is a finite identity and consequence audit, not an unlabeled-poset census or proof of the universal conjecture.
