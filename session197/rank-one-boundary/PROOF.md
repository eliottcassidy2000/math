# Signed root-deletion boundary and its rank-one endpoint factorization

Session 197. All-order theorem for arbitrary rooted tournaments; no symmetry assumption. This supplements, and does not change, the previously frozen mirror-source boundary results.

## 1. Root-gap polynomial and a signed Hamilton sum

Let Q have root r, exterior D=Q-r, O=N^+(r), and I=N^-(r). For k>=1 let a_k count vertex words using every exterior vertex once and r exactly k times, with no adjacent r's and every consecutive pair an arc. Put a_0=H(D), including H(empty)=1, and define R_Q(z)=sum a_k z^k. For positive transitive clone length m,

    H(Q_m)=sum_(k>=1) a_k k! S(m,k).

Thus R_Q(-1) is a boundary correction, not the number of Hamilton paths of a tournament with a negative number of clones.

On the exterior vertices define a signed matrix M with zero diagonal, the usual tournament adjacency matrices within O and I, and cross blocks

    M_(O,I)=B,    M_(I,O)=-B^T,

where B_(u,v)=1 exactly when u->v for u in O,v in I.

**Lemma.** If D is nonempty, R_Q(-1) is the sum of weights of all Hamilton words of M starting in I and ending in O, where a word's weight is the product of its consecutive matrix entries. If D is empty, R_Q(-1)=0.

To prove this, fix an ordering v_1,...,v_N of D and sum over possible insertions of r in the N+1 gaps. Write o_v=[r->v] and i_v=[v->r]. The start and end factors at z=-1 are 1-o_(v1) and 1-i_(vN). The factor for an internal gap u,v is

    [u->v] - i_u o_v.

The first two factors force the word to start in I and end in O. The internal factor is ordinary adjacency except from I to O, where it equals [u->v]-1=-[v->u]. These are exactly the entries of M. This expansion counts both the no-root term a_0 and every admissible root-gap choice, with its correct sign. For a singleton Q, R_Q(z)=1+z directly.

This identity applies to every crossing matching rank, without an assumed cancellation of exponential bases.

## 2. Exact rank-one factorization

Let G be the bipartite graph of O-to-I arcs, and suppose its maximum matching size is one. It has a one-vertex cover c: all its edges meet c. This also follows directly by choosing an edge; if there were edges avoiding each of its endpoints, they would either be disjoint from it or yield two disjoint edges with each other.

Write s_v(T) and e_v(T) for the counts of Hamilton paths of T starting and ending at v, respectively.

**Theorem.** For any such cover center c:

- If c is in O, then

      R_Q(-1)=-s_c(Q[O]) sum_(v in I, c->v) e_v(Q[I]).   (1)

- If c is in I, then

      R_Q(-1)=-e_c(Q[I]) sum_(u in O, u->c) s_u(Q[O]).   (2)

In particular R_Q(-1)<=0 for every rank-one source. If more than one cover center is available, all the displayed products agree.

For (1), every nonzero cross entry of M meets c. A Hamilton word of M starting in I can first enter O only at c. After this entry it cannot return to I and later end in O, because that would require visiting c again. Hence it consists of a Hamilton path through all I, then c, then all remaining O. Its unique cross transition is v->c, with weight -1 exactly when c->v in Q. The I-side endpoint and O-side starting path are independent, giving (1). Reversing the roles proves (2). Either factor may be zero; no unjustified positive-endpoint assertion is used.

For matching rank zero, G has no edges and there is no nonzero I-to-O entry of M, so R_Q(-1)=0. This includes empty O or I and the singleton boundary separately handled above.

## 3. Binary-clock interpretation and limitations

At rank one the only positive root-run coefficients are a_1,a_2. Put A=a_2>0 and beta=a_1-2a_2. Then

    H(Q_m)=A 2^m+beta,    m>=1,
    R_Q(-1)=H(D)-A-beta.

Writing K for either nonnegative endpoint product in (1) or (2), the exact deletion identity is

    H(D)=A+beta-K,    beta=H(D)-A+K.                    (3)

This extends the negative-square boundary of the root-fixed self-converse rank-one construction to arbitrary rank-one sources. In a symmetric single-cross-edge construction, the two endpoint factors are dual and K is their square. Symmetry is not needed for (1)-(3), and arbitrary endpoint products need not be squares.

Nonpositivity does not extend to rank two. Two literal hostile examples, with adjacency rows and distinguished root, are:

- Rows [16,5,9,3,14], root 1: (a_0,a_1,a_2,a_3)=(5,15,11,2), so R_Q(-1)=-1.
- Rows [2,20,9,3,13], root 1: (a_0,a_1,a_2,a_3)=(1,11,17,4), so R_Q(-1)=3.

Both have crossing matching rank two. Thus the general signed identity does not imply a sign-definite rank-two boundary, a zero boundary for unit ternary clocks, or a unit-offset congruence. Those questions remain separate.

## 4. Standalone checks

`verify.cpp` has no external code or data dependency. One DP counts all exterior permutations and allowed root gaps to obtain R_Q. A separate signed adjacency DP counts Hamilton words starting in I and ending in O directly. A third endpoint DP on each side evaluates (1) and (2) for every available cover center.

The verifier regenerates all unlabelled tournaments through order 7 by all vertex extensions and exact degree-class permutation canonicalization, independently lists every labelled tournament through order 5, and adds twelve seeded graphs at each order 8 through 13. The bank sizes 1,1,2,4,12,56,456 are checked. These universes give:

- 9,774 roots checked against the general signed identity.
- 4,430 rank-one rooted tests and 6,241 center-factorization checks.
- 1,428 rank-zero boundary checks.
- Both signs explicitly present at rank two, including the two displayed examples in the complete small bank.

The overlap between labelled and unlabelled controls is intentional. The larger seeded controls are not an exhaustive census. Optimized and debug C++17 runs produce identical receipts. All checks use exceptions, not assertions. Counts at orders at most 13 fit signed 64-bit integers, using the bound (n-1)! 2^n for the sum of absolute root-gap contributions. Compile and run `verify.cpp` directly; no symbolic algebra package is required. The finite checks support the all-order gap expansion and path decomposition above rather than replacing them.
