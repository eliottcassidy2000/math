# A composable endpoint-typed path-cover algebra

Session 195, 2026-11-05. **PROVED elementary composition and representation theorems; FINITE-EXACT independent controls.** The full tournament Hamilton-path spectrum with exceptions `{7,21}` remains **OPEN**. This packet supplies a physical recursive interface, not a proof that its attainable values cover the required spectrum.

## 1. Inheritance, scope, and the information retained

The starting guide is [Session 194](https://github.com/eliottcassidy2000/math/blob/188bb29fa08b6dfd51ab94cc8e494d6d56f5e9ec/00-navigation/S194.md). Its [module barrier](https://github.com/eliottcassidy2000/math/blob/188bb29fa08b6dfd51ab94cc8e494d6d56f5e9ec/session194/module-barrier/PROOF.md) proves that an intact nontransitive input module cannot produce literal `3H+2` or `3H+4`. Its [endpoint floor](https://github.com/eliottcassidy2000/math/blob/188bb29fa08b6dfd51ab94cc8e494d6d56f5e9ec/session194/endpoint-floor/PROOF.md) excludes rearming a unit endpoint while retaining the specified induced strong cores. Its [two-apex analysis](https://github.com/eliottcassidy2000/math/blob/188bb29fa08b6dfd51ab94cc8e494d6d56f5e9ec/session194/two-apex/PROOF.md) explains why periodic growth can retain nonconstant lower-order counters.

The closest inherited positive mechanism is the ordinary path-cover/walk-content substitution kernel, THM-3121, `01-canon/theorems/THM-3121-path-cover-walk-content-substitution-kernel.md`, named in the module-barrier proof. Here the inter-piece arcs need not be uniform on an entire operand. Instead, the actual endpoints of each path fragment retain their **cross-neighborhood types**. The new construction interface is deliberately finite in number of types but unbounded in polynomial degree. No literature-priority claim is made.

Live concepts: maximal runs, endpoint types, fragment connectivity, nonmodule growth, persistent counters, and physical color recoloring. A numerical root statistic is useful only to the extent that the next operation can read it. The hostile examples are an induced prime H13 bank, arbitrary retained strong cores, and a directed cycle of fragments mistakenly admitted by degree-only bookkeeping.

A directed path is a nonempty ordered vertex word with every consecutive arc directed forward. A path cover is an **unordered set** of vertex-disjoint directed paths covering the whole vertex set. Singleton paths are allowed; empty paths are not. There are no cycle components. The empty digraph has the unique empty cover.

Give every vertex one color in `[k]={0,...,k-1}`. Define

    F_T(z) = sum over path covers P of T, product over paths p in P z_(color(first p),color(last p)).   (1)

Thus a coefficient `f_T(m)` counts covers with exactly `m_ij` paths of endpoint type `(i,j)`. For nonempty T, the sum of all coefficients of total degree one is `H(T)`. Forgetting all colors by replacing every `z_ij` by z gives the ordinary unordered path-cover polynomial. The full profile, not only its degree-one part, is retained.

## 2. The physical gluing operation

Take disjoint k-colored tournaments A and B. Choose a binary matrix `C` of size k by k. Add every cross arc by the rule

    x in A, y in B: x -> y iff C_(color x,color y)=1.       (2)

The opposite arc is used when the entry is zero. No old arc is changed. The rule is uniform on **color classes across this particular cut**, not necessarily on either whole operand.

After creating the cross arcs, recolor A by a map `alpha:[k]->[l]` and B by a possibly different map `beta:[k]->[l]`. The two maps may merge colors; they cannot split a color class inside an operand. Operand side is remembered while this recoloring is performed. A parent may thus retain only A's root, only B's root, or both roots in distinct colors, as its chosen output palette permits.

Write this tournament operation as `G_(C,alpha,beta)(A,B)`. Binary gluing plus leaf creation and such recolorings are actual graph constructions. They are not merely formal arithmetic operations on requested Hamilton counts.

### Theorem 1: exact all-cover binary closure

There is an explicit bilinear operation on profile polynomials, depending only on `k,l,C,alpha,beta`, such that

    F_(G(A,B)) = F_A star_(C,alpha,beta) F_B.             (3)

In particular, two inputs with the same full colored profile have the same profile in every subsequent expression made from these operations. This is **counting-context completeness**, not an assertion that the profile determines the tournament up to isomorphism or is a minimal possible signature.

**Combinatorial kernel.** Fix local cover multiplicities `m` for A and `n` for B. Treat every actual local path as a distinct fragment. A fragment has type

    t=(side, start color, end color),   side in {A,B}.

There are `2k^2` possible types. On the distinct fragments put an arc from fragment t to fragment u exactly when their sides differ and the endpoint-to-endpoint cross arc is available under (2). Thus, for `t=(A,i,j)` and `u=(B,p,q)`, the arc t->u exists iff `C_jp=1`; the arc u->t exists iff `C_iq=0`. These are independent tests on different original endpoint pairs. In particular, the fragment digraph need not be a tournament and can contain a pair of opposite arcs. This causes no problem: only acyclic path covers of the fragment digraph are counted.

A path cover of this fragment digraph gives a path cover of the glued tournament by concatenating each listed sequence of local paths. Its output path type is determined by the recolored start of the first fragment and recolored end of the last. Let `K_(m,n)(z)` be its profile polynomial. Then

    F_(G(A,B))(z) = sum_(m,n) f_A(m) f_B(n) K_(m,n)(z).   (4)

**Proof.** Given any path cover of G(A,B), cut each path precisely at the cross edges. The maximal resulting A-runs form an unordered path cover of A, and the B-runs form an unordered path cover of B. Each original global path becomes a directed sequence of distinct fragments with alternating sides. These sequences form a path cover of the fragment digraph. Conversely, local covers and such a fragment path cover concatenate legally into a global cover. Cutting and concatenating are inverses because every interfragment arc crosses the cut, while all internal fragment arcs remain inside one operand. The covers and fragments themselves are actual named combinatorial objects; repeated endpoint types do not identify them. The type multiplicities determine the fragment digraph up to an endpoint-tag-preserving isomorphism, so its count is independent of which particular local cover was chosen. This proves (4) and hence (3). QED.

The associativity statement is physical: two expression trees that build the same colored graph have the same output profile. Arbitrary cross matrices with incompatible edges do not define a single unqualified associative product.

### A finite-variable rational/exponential specification

Introduce one variable `x_t` for each fragment type. Let D be the `2k^2` by `2k^2` allowed-fragment-transition matrix above, and `X=diag(x_t)`. Let `s(t)` and `e(t)` be the start/end colors after alpha or beta. For output colors i,j put

    W_ij(x) = sum_(t:s(t)=i, u:e(u)=j) [X(I-DX)^(-1)]_(t,u),
    W(x,z) = sum_(i,j) z_ij W_ij(x).                    (5)

The inverse is a formal power series with constant matrix I. W enumerates nonempty type words of directed fragment paths, recording their output endpoint type. Distinct fragment names are restored by factorials:

    K_(m,n)(z) = (product_ij m_ij! n_ij!)
                [product_ij x_(A,i,j)^m_ij x_(B,i,j)^n_ij] exp(W(x,z)).   (6)

The exponential makes an unordered set of paths. For a multiplicity pattern its `1/r!` factors remove the ordering of path components, while the displayed factorials distribute the named fragments of each type among their occurrences. This is the multitype labelled-set counting identity; it can also be proved by expanding W into words and counting each labelled cover. Although rational coefficients occur during expansion, (6) is an integer polynomial by the bijection.

There is no cycle series in (5)-(6). Imposing only indegree/outdegree at most one would wrongly admit directed cycle components. For example, a directed triangle has profile `3z+3z^2+z^3`, not an additional constant contribution from its directed cycle. A closed fragment cycle cannot later become part of a spanning path without deleting an already selected arc, which the operation does not do.

For fixed k and l the number of variables is fixed, and all coefficient ranges are bounded by the input order. Consequently brute formal coefficient arithmetic in these fixed numbers of variables gives a polynomial-in-order arithmetic algorithm, with a degree depending on k,l. The full profile has at most `(n+1)^(k^2)` terms; this is a loose bound, not a claim that all terms occur. The exact integers have `O(n log n)` bits since the number of path covers is at most `n! 2^(n-1)`. The companion verifier intentionally uses small explicit labelled-fragment subset DP rather than implementing the asymptotically polynomial fixed-type coefficient engine.

## 3. An exact structural interpretation of the interface width

This section **defines** the width used here; it does not appeal to a potentially differently normalized standard graph invariant.

For a nonempty vertex set U in a tournament T, two vertices of U have the same external type if they have the same orientation toward every vertex outside U. Let `b_T(U)` be the number of those equivalence classes; set `b_T(V(T))=1`. A binary leaf decomposition of T is a rooted full binary tree whose leaves are the vertices, with each node representing its descendant vertex set U. Its **boundary-type width** is `max_U b_T(U)`. The boundary-type width of T is the minimum over such trees.

### Theorem 2: tree description equals physical colored construction

A tournament has boundary-type width at most k if and only if it can be built by binary operations (2), with at most k colors at the output of every operation and side-specific color maps as above.

**Forward proof.** At every node U color its vertices by their external types relative to T minus U. For children A,B, the arcs between A and B are uniform on their external-type classes: if x,x' have the same child-A type they agree toward every y in B; if y,y' have the same child-B type they agree toward every x in A. Thus (2) holds. A child's external type determines its parent external type because the set of external test vertices only shrinks. This gives the required side-specific color maps. Leaves and the root satisfy the stated conventions. QED.

**Reverse proof.** Fix the binary expression tree. At a node U, vertices having the same output color will always receive the same arcs toward every vertex introduced outside U: later color maps cannot split them, and every later cross arc depends only on colors. Thus each current color class lies inside one final external-type class. Hence `b_T(U)` is at most the number of current colors, at most k. QED.

This criterion says exactly when a chosen decomposition can pass only k types at each cut. It is not a small vertex separator: the underlying undirected graph of a tournament is complete. The boundary types encode dense external incidence patterns rather than a set of k exceptional edges or vertices.

**Representational limitation.** For each fixed k, this construction class does not contain every tournament. For n leaves, there are fewer than `4^n` ordered binary tree shapes, at most `k^n` leaf colors, and at most `2^(k^2) k^(2k)` choices of cross rule and side-specific recolorings per node, using a fixed k-label palette. Labelling the leaves gives at most `n! C_k^n` labelled tournaments for a constant C_k. This is less than the `2^(n(n-1)/2)` labelled tournaments for all sufficiently large n. Extra unary recolorings can be absorbed into a neighboring binary operation, so they do not evade the count.

This does **not** prove that bounded-width tournaments miss any Hamilton count: many nonisomorphic tournaments have the same count, and a much smaller family could still meet the entire numerical spectrum. Nor does the width bound allow arbitrary edge data to be hidden in a node: a node stores only the binary matrix C and the two color maps. Allowing arbitrary cross edges would invalidate both the theorem's definition and the counting argument.

## 4. Unbounded linear state rank already at width two

The polynomial's unbounded degree is not merely an inconvenient implementation choice if one requires a universal **linear** continuation interface.

### Theorem 3: no fixed-dimensional bilinear counting algebra

There is no fixed d and assignment of rational vectors `s(A) in Q^d` to all colored tournament inputs for which every allowed gluing context has its exact Hamilton count as a linear function of s(A). Consequently there is no fixed-dimensional rational vector algebra representing every width-two gluing operation bilinearly and reading H linearly. The result already holds on the very simple inputs `A_n=Tr_n`, with all n vertices in one color.

**Proof.** Use the identity context `C_0[A]=A`. For m>=1 use

    C_m[A] = C3[A,Tr_m,Tr_m],

with cyclic interblock arcs `A -> B -> D -> A`. These are width-two contexts on monochrome transitive inputs: first join A to B while retaining their two distinct colors, then join the monochrome D using the rule `B->D->A`. No operand's color needs to be split.

Let `c_r(A)` count ordered r-path covers, so it is r! times the coefficient of z^r in the one-color F_A. Cutting a Hamilton path into maximal block runs gives

    H(C_m[A]) = sum_(r=1)^(m+1) kappa_(m,r) c_r(A),
    kappa_(m,m+1) = (m!)^2 > 0.                      (10)

Between consecutive A-runs the cyclic quotient requires a nonempty B-run and a nonempty D-run. Since each exterior block has m vertices, there are at most m+1 A-runs. At equality the quotient word is uniquely `A(BDA)^m`; both exterior blocks are split into m singleton paths, each of which can occur in m! orders. This proves the leading coefficient. All remaining coefficients are nonnegative integers obtained by the same maximal-run bijection. For m=0 set `kappa_(0,1)=1`.

Every set partition of a transitive n-tournament into r nonempty parts gives exactly one unordered path cover, so

    c_r(A_n) = r! S(n,r),

where S(n,r) is the Stirling number of the second kind. It vanishes for r>n and has diagonal value `c_n(A_n)=n!`.

For any N form the N by N matrix with entries `M_(n,m+1)=H(C_m[A_n])`, for `1<=n<=N` and `0<=m<N`. Equation (10) factors M into the matrix `(c_r(A_n))` and the transpose of `(kappa_(m,r))`. Both factors are triangular with the stated nonzero diagonals. Hence

    det M = product_(n=1)^N n! ((n-1)!)^2 != 0.       (11)

If every column were a linear readout of the same d-dimensional signatures s(A_n), every such matrix would have rank at most d. Since N is arbitrary this is impossible. A bilinear gluing algebra with linear final readout would make every fixed context a linear readout of the input vector, so it is excluded as well. QED.

The first context rows are `c_1`, `3c_1+c_2`, and `11c_1+17c_2+4c_3`. Direct vertex DP and independent cyclic quotient-word counts check the first five by five matrix, whose determinant is 2,866,544,640.

This is a rank obstruction with specified quantifiers. It does not forbid a fixed scalar encoding with arbitrary nonlinear decoding, a nonlinear transfer law, or a finite-dimensional linear state for a restricted future language such as one particular fixed-gadget attachment process. The fixed-apex matrices and fixed-bank rooted matrices are therefore entirely compatible with it. Nor does it imply numerical spectral incompleteness at width two.

## 5. Infinite prime growth with only two types

Here **prime** means module-indecomposable: there is no module M with `2<=|M|<|V(T)|`. This is not numerical primality of H.

For `r>=1`, define Q_r on `a,x_1,...,x_(2r)` by

    x_i -> x_j whenever i<j;
    a -> x_i for odd i;
    x_i -> a for even i.                              (7)

Retain color 0 for a and color 1 for the whole chain. Start from singleton a. Append each new chain vertex with every old chain vertex beating it; the apex beats odd additions and loses to even additions. This uses (2) with two colors. Appending a pair gives Q_(r+1) while retaining Q_r as an induced subtournament.

### Theorem 4: indefinite prime recursion and its exact counter

For every r>=1, Q_r is prime and strong, has boundary-type width at most two, and

    H(Q_r) = (4^(r+1)+6r+5)/9.                         (8)

With `h_0=H({a})=1`, the actual recursive growth satisfies

    (h_r,r) -> (4h_r-2r-1, r+1).                       (9)

Thus the values begin `1,3,9,31,117,459,1825,...`. The construction grows prime tournaments indefinitely rather than repeatedly substituting intact nontrivial modules. Nevertheless its scalar update carries a nonconstant counter.

**Primality.** A module not containing a must be an interval of the increasing chain: an exterior chain vertex between two module vertices would distinguish them. An interval of at least two chain vertices contains consecutive opposite-parity vertices, which a distinguishes. Hence no such nontrivial module exists.

Suppose a proper module M contains a and at least one chain vertex. Let p and q be the smallest and largest chain indices in M. Any even chain vertex j outside M beats a, so it must beat every chain vertex in M, forcing `j<p`. Any odd chain vertex j outside M loses to a, so it must lose to every chain vertex in M, forcing `j>q`. If p>1, exterior odd vertex 1 contradicts the second requirement; hence p=1. If q<2r, exterior even vertex 2r contradicts the first; hence q=2r. Now neither exterior condition is possible at any chain index, so M is the whole tournament, a contradiction. QED.

**Strongness.** The directed Hamilton cycle `a,x_1,...,x_(2r),a` is present. In particular, every vertex reaches every other along this cycle. QED.

**Direct count.** Every Hamilton path consists of an increasing chain subword P, then a, then its increasing complementary subword R. If P is nonempty its largest index must be even; if R is nonempty its smallest index must be odd. Empty P and empty R supply one path each. If both are nonempty, write `i=min R` and `j=max P`.

When i>j they must be consecutive, and the valid splits are `j=2,4,...,2r-2`, contributing r-1. When i<j write `i=2s-1`, `j=2t`, `1<=s<=t<=r`. All indices below i are forced into P, all above j into R, i and j are fixed in their respective subwords, and each index strictly between them is free. This gives `2^(j-i-1)=4^(t-s)` paths. Therefore

    h_r = r+1 + sum_(d=0)^(r-1) (r-d)4^d
        = (4^(r+1)+6r+5)/9.

This also gives (9) by exact subtraction; no fitted samples or unproved recurrence are used. QED.

The formula is an illustrative elementary control, not a claim of a newly discovered tournament family. Its counts form a sparse sequence and do not supply arbitrary prime or 7p atoms.

## 6. Three-type births and heterogeneous joins

A physical two-apex tournament has named a,b with a->b and increasing chain `x_1,...,x_m`. A letter `c_i in {0,1,2,3}` records `a->x_i` in bit 0 and `b->x_i` in bit 1. Give a,b,chain colors 0,1,2. Appending a chain vertex is a three-color gluing step, whatever its letter. Hence the entire two-apex class has boundary-type width at most three in the present definition. This does not remove the Session 194 fixed-macro counter barrier.

In outgoing-bitmask notation, where bit j of row i means i->j, the physical prefix `120` is

    [6,8,26,17,3], H=13.

Appending letter 1 gives

    [38,8,58,49,35,2], H=23.

Both tournaments are prime by exact module enumeration. This particular example is **FINITE-EXACT**, not an all-prefix primality claim. The retained H13 input is not a module in its child: the new vertex sees its two apex colors differently. The entire profile transfer, not just H, is checked.

As an explicit heterogeneous composition, take Q_2 with H=9, root a=0 and rows

    A=[10,28,25,16,1], colors (0,1,1,1,1),

and the two-apex H13 bank

    B=[6,8,26,17,3], colors (0,1,2,2,2).

Orient every cross arc A->B except the root-root arc, which points from B's vertex 0 to A's vertex 0. This is C with every entry one except `C_00=0`. The actual ten-vertex tournament has

    H=279.

Its full typed profile has 946 monomials in the displayed three-color palette. The gluing count equals an independent vertex-subset Hamilton DP on the constructed ten-vertex graph. These inputs came from different growth descriptions; the interface composes them without requiring a common scalar recurrence. This is one illustration, not a density or coverage result. Recoloring after this join can retain one root for another join, or retain more colors if the next operation needs them.

The same operation handles a root-root arc reversal between arbitrary rooted operands by marking the root separately from the bulk. The rooted-transfer analysis in the parallel Session 195 work may compress some fixed-gadget continuations further. The full cover profile here is a universal sufficient interface for all color-uniform gluing continuations, not a claim that every restricted continuation needs all of it.

## 7. What remains for the spectrum proof

The algebra removes a false dichotomy: recursive pieces need not be whole modules, and a prime tournament need not be an irreducible object for a bounded-type construction. It also separates three tasks that cannot be conflated:

1. **Physical representation:** build an actual graph expression of bounded or explicitly growing type width.
2. **Exact evaluation:** transport the full endpoint-typed path-cover profile, or prove a smaller closed quotient for the chosen future operations.
3. **Inverse targeting:** for every required integer atom, find an expression whose degree-one count is exactly that integer.

Theorems 1, 2, and 4 accomplish the first two for their stated classes, not the third. No all-order inverse-target theorem is proved. The module threshold and induced-core endpoint floor remain valid wherever their hypotheses hold; nonuniform cross arcs escape the module hypothesis but do not escape an endpoint floor while retaining its obstructing core. In particular, a new expression need not have a unit endpoint merely because its color boundary is small.

The precise next question is whether one can identify a finite collection of physically realizable profile classes with return operations that control the additional coordinates and permit exact prime/7p targeting. A numerical affine identity without a closed profile class remains insufficient. Conversely, our counting bound against representing every tournament does not rule out such a spectral classification.

## 8. Reproduction and verification scope

Run from this packet directory:

    python3 verify.py > normal.stdout
    python3 -O verify.py > optimized.stdout
    cmp normal.stdout optimized.stdout
    python3 audit_literal.py > literal_normal.stdout
    python3 -O audit_literal.py > literal_optimized.stdout
    cmp literal_normal.stdout literal_optimized.stdout
    python3 check_rank.py > rank_normal.stdout
    python3 -O check_rank.py > rank_optimized.stdout
    cmp rank_normal.stdout rank_optimized.stdout

The main checker verifies all k=2 colorings of every labelled tournament of orders one through four, all nonempty proper unordered cuts whose cross edges are color-uniform, the full all-cover composition, Hamilton readout, and color-forgetting. It separately checks 256 side-specific recoloring cases, recursive full profiles through seven chain births, Q_r scalar/primality controls through r=6, the H13-to-H23 full-profile birth, and the heterogeneous ten-vertex join.

The independent literal audit enumerates every vertex permutation and every cut pattern, selecting the lexicographically ordered path list to count each unordered cover once. It matches all 1,098 binary-colored tournaments through order four, selected two-apex banks through six vertices, and 16 tagged fragment kernels. It includes the directed-cycle hostile. This independent expected-profile path does not use the subset-cover DP. All validity gates are explicit exceptions and remain active under `python3 -O`.

The rank checker independently compares 25 actual cyclic three-block tournaments through order thirteen against Stirling cover vectors and cyclic quotient-word enumeration. It verifies each leading coefficient and all five leading principal determinants.

The all-order conclusions are established by the proofs above. The finite checks are controls with an explicit universe, not substitutes for proof. The scratch directory is separate from the unrelated dirty checkout; no existing repository file was edited.
