# The exact padding complex for tournament embedding parity

2026-10-19. **PROVED elementary identities**, with finite exact controls and [independent review](review/INDEPENDENT_REVIEW.md), plus a separate [algebra review](REVIEW_ALGEBRA.md). No external-priority claim. This gives an operator explanation for padding cancellations; it is not a classification of all good trees.

## 1. Operators and their domains

Let C_n be the vector space over F_2 of isomorphism-invariant functions on n-vertex tournaments, including the unique empty tournament at n=0. Put C_(-1)=0 and let d_(-1) be the zero map. Let t_n=dim C_n, the number of tournament isomorphism classes.

Define d_n:C_n→C_(n+1) by

    (d_n f)(T)=sum_{v in V(T)} f(T-v),

and h_n:C_n→C_(n-1), n>=1, by

    (h_n f)(T)=f(T with one new universal source).

A universal source points to every old vertex. Put h_0=0. All sums are in F_2. These are well-defined on isomorphism classes.

## 2. Contracting identity

**Theorem.** For every n>=0,

    d_(n+1)d_n=0,
    d_(n-1)h_n+h_(n+1)d_n=identity on C_n.

The first equality holds because every unordered pair of deleted vertices occurs in two orders. For the second, evaluate h_(n+1)d_n f on T: deleting the new source contributes f(T); deleting an old vertex contributes f((T-v)+source). The sum of these latter terms is exactly d_(n-1)h_n f, so they cancel. At degree0 the source deletion is the only term, and the same identity holds with C_(-1)=0.

Consequently

    kernel(d_n)=image(d_(n-1)).

Moreover, if d_n f=0, the explicit preimage is h_n f. Thus the kernel is not merely given by a dimension count: it has a canonical certificate using source adjunction.

This is a contracting homotopy on the full graded tournament-function spaces. It must not be confused with a claim that every graph-level primitive retains the shape of a tree.

## 3. Exact action on embedding parities

For an oriented graph S on n vertices, write phi_S(T) for its number of spanning injective embeddings into T modulo2. Then

    d_n phi_S = phi_(S disjoint union one isolated vertex),
    h_n phi_S = sum_{u with indegree_S(u)=0} phi_(S-u).

For d, choose the image of the isolated vertex. For h, the preimage of the added universal source must have indegree zero; every such choice leaves precisely an embedding of S-u. Isolated vertices are counted among these sources.

These identities explain why retaining an isolated vertex in a deletion–reversal calculation matters. Two isolates give zero by d²=0. A single isolate can erase nonconstant parity information, exactly the functions in the image of d at the preceding degree. The leaf-pair and padded-double-path examples in the accompanying tree work are concrete graph-level manifestations.

Tournament embeddings themselves span the entire function space: if S is a tournament, phi_S is the indicator of its isomorphism class. Indeed, its nonzero value is |Aut(S)| modulo2, which is1. A tournament has no automorphism of order2, because a transposed pair would reverse its directed edge. Cauchy's theorem then makes its automorphism group odd. Thus no extra subspace restriction is hidden in the function-space theorem.

## 4. Exact kernel and image dimensions

Let r_n=rank d_n and r_(-1)=0. Exactness gives

    t_n=r_(n-1)+r_n,
    r_n=sum_{j=0}^n (-1)^(n-j)t_j.

In particular dim kernel d_n=r_(n-1). The formula is in ordinary integers, although the matrices are over F_2. The finite verifier independently constructs all tournament classes through five vertices and checks the operator identity and matrix ranks; it does not assume the formula in its elimination routine.

The first class counts are t_0,...,t_5=1,1,1,2,4,12. Thus ranks d_0,...,d_4 are1,0,1,1,3, and the formula predicts rank d_5=9. The rank d_5 prediction is a theorem from exactness and t_5; the supplied direct matrix audit only reaches d_4 because it stops at five-vertex hosts.

## 5. A sharp parity-position restriction

Let 1_n denote the constant-one function. Then

    d_n 1_n=(n+1 mod2)1_(n+1).

If n is odd, d_n f cannot be the constant-one function: applying d_(n+1) would give zero on the left but one on the right, because n+2 is odd. Therefore an even-order oriented graph with an isolated vertex can never be a Rédei graph (odd embedding count in every same-order tournament). This particular obstruction also has the simpler direct witness of the transitive host: all one-vertex deletions are isomorphic and there are an even number of them. The contracting identity supplies the broader kernel description, rather than being necessary for this elementary corollary.

For even n, suppose d_n f=c·1_(n+1), with c in F_2. The contracting identity gives the exact decomposition

    f=c·1_n+d_(n-1)h_n f.

Thus a padded-odd even-order source differs from the constant-one function by a precisely specified padding-invisible term. This statement allows f itself to be nonconstant, as the eight-vertex double-path example demonstrates.

The obstruction is sharp in order: a directed path on any positive even number n of vertices has constant odd embedding parity by Rédei's theorem; see the [primary modern treatment](https://arxiv.org/html/2510.10659v2), Sections 1 and 3, and [THM-4533](../../../01-canon/theorems/THM-4533-redei-graphs-parity-of-shaved-tournaments-and-u9.md). Adjoining one isolated vertex multiplies the subset parity by n+1, which is odd. Hence odd-order padded Rédei graphs exist at every odd order at least3. The one-vertex isolated graph supplies the remaining order1 example from the empty graph.

## 6. What the connection preserves

Source: parity functions of tournament embedding counts. Target: a graded linear complex with source-adjunction contraction. Preserved: every parity value, exact deletion/padding operators, kernels and dimensions. Lost if one insists on a single graph representative: h produces a sum of vertex-deleted graphs, and a primitive need not be a tree. No integer counting identity beyond modulo2 is inferred from the cancellation.

This is an exact operator connection, rather than a numerical analogy between a graph and a differential equation. It gives a reusable way to organize padding certificates and to distinguish automorphism-forced zeros from tournament-specific cancellations. General good-tree classification and extremal Rédei size remain open.
