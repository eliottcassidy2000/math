# Consecutive unique-start roots and the exact start-two grammar

Session 197. **All-order refinement of the C=1 structure and a complete recursive classification of rooted start count two.** This does not classify every C=2 tournament. The completed C=1 theorem is in `C1_COMPLETE.md`; the one-seam conventions and Fibonacci enumeration are in `ROOTED.md` and `FIBONACCI_CLASS.md`.

Let s_v(T) count Hamiltonian paths starting at v. A nontrivial strong tournament satisfies `1<=C(T)<=s_v(T)` for every vertex v, by cutting Hamiltonian cycles. Therefore a strong graph with a start-two root has C=1 or C=2.

## 1. The start-one vertices form one short interval

**Theorem 1.** If T has a unique Hamiltonian cycle and order n>=4, the vertices with s_v(T)=1 form a consecutive block of either two or three vertices on that cycle. All but the last vertex of this block have outdegree one, and the last has outdegree at least two. These are all the outdegree-one vertices of T. The triangle is the exception: all three vertices have start count one and outdegree one.

The proof is by the canonical rooted grammar, not an inference from the Fibonacci census. We give the additional path-counting details needed for the induction.

Write a canonical construction as K(A,a;G,b), with A a singleton or strong with s_a(A)=1. The new gadget is b plus a transitive tail W=(w_1,...,w_m); if m>0 its last vertex points to b. Put epsilon_i=1 when b->w_i. Thus epsilon_m=0. Every cross arc is A->G except b->a. The new unique cycle is b, the old unique a-start path, W, b.

For x in A other than a, let p_(x,a)(A) count Hamilton paths starting at x and ending at a. Let D_(x,a)(A) count ordered spanning covers by two nonempty paths, the first starting at x and the second at a. Cutting the old Hamilton cycle at x and a proves D_(x,a)>=1 when A is nontrivial. Let g=H(G), t=s_b(G), and L be the gadget's root-ending two-cover count from `ROOTED.md`. Directly splitting paths according to use of the seam gives

    s_x(K)=g*s_x(A)-t*p_(x,a)(A)+L*D_(x,a)(A).        (1)

Without the seam, the A-path is followed by a G-path, except for the forbidden a->b join. With the seam, the two A-runs start at x and a and the G-run before the seam ends at b. This proves (1), including the empty final G-run convention. At x=a, the seam cannot be used; its start count is g when A is nontrivial and g-t when A is a singleton.

### Nonempty tail

Suppose m>=1. Then L>=2. If t=0, (1) is at least s_x(A)+L>=3. If t=1, G contains a triangle, so g>=3; using p_(x,a)<=s_x(A), (1) is at least (g-1)s_x(A)+L>=3. Thus no x in A\{a} remains a start-one vertex.

The old root a has new start count one exactly when all epsilon_i are zero. Indeed that is exactly when G is transitive and g=1,t=0; if any epsilon_i is one, the final zero gives a directed triangle in G and g-t>=2.

The new root b always has start count one. For a tail vertex, every path starts with a G-fragment from w_i to b, then uses the unique a-start path in A, and finishes with the remaining increasing tail. Therefore

    s_(w_i)(K)=[epsilon_i=0]
              +sum_(j>i,epsilon_j=0) 2^(j-i-1).      (2)

The last tail vertex has start count one. The penultimate tail vertex does too exactly when its epsilon is one. Every earlier tail vertex has start count at least two from the final-zero summand alone.

If all bits are zero, the start-one block is w_m,b,a. Otherwise it is w_m,b, with w_(m-1) prepended exactly when epsilon_(m-1)=1. These are consecutive on the unique cycle. Their outdegrees have the asserted form: w_m points only to b; a penultimate vertex with epsilon one points only to w_m; b has outdegree one exactly in the all-zero case. In that all-zero case a has outdegree at least two unless A is a singleton and m=1, which is precisely the triangle.

### Empty tail

Now m=0, so G is the singleton b and A is nontrivial. The new root b points only to a, and a retains start count one. Formula (1) becomes

    s_x(K)=s_x(A)-p_(x,a)(A)+D_(x,a)(A),   x!=a.      (3)

Since D_(x,a)>=1, an extra start-one vertex requires s_x(A)=p_(x,a)(A). The old Hamilton path obtained by cutting the unique cycle at x then ends at a. Hence x is the cycle successor of a, and every x-start path ending at a closes to that unique cycle; therefore s_x(A)=1. By the inductive block-and-outdegree assertion, such an adjacent pair a,x of start-one vertices occurs exactly when a has outdegree one. This also holds in the triangle.

Conversely, if a has outdegree one, its sole out-neighbor x has start count one: any x-start Hamilton path must visit a last and therefore closes to the unique cycle. In a cover counted by D_(x,a), the path beginning at a must be the singleton a, since its only possible next vertex x belongs to the other path. The other path is uniquely determined, so D_(x,a)=1. Thus (3) is one exactly at this x.

The new block is consequently b,a when a's old outdegree exceeds one, and b,a,x otherwise. The old a has unchanged outdegree; every other old vertex gains an outgoing arc to b. This proves the asserted outdegree pattern. The induction begins at the triangle, and the complete C=1 grammar covers every graph under consideration. Finally any outdegree-one vertex has s_v=C=1, so no such vertices lie outside the displayed block. QED.

## 2. Exact refined Fibonacci counts

Let a_(n,2) and a_(n,3) count unlabelled strong C=1 graphs with two or three start-one vertices, for n>=4. The already proved C=1 enumeration and rigidity give

    a_(n,2)+a_(n,3)=F_(2n-6),
    2a_(n,2)+3a_(n,3)=F_(2n-4).

Thus

    a_(n,2)=F_(2n-8),       a_(n,3)=F_(2n-7).         (4)

Let d_n count a C=1 graph with an **unordered pair** of distinct start-one vertices, up to isomorphism. At n>=4 the graph is rigid, so Theorem 1 gives

    d_n=a_(n,2)+3a_(n,3)=2F_(2n-4)-3F_(2n-6).

At order three, the triangle's three pairs are one orbit, so d_3=1. Put Delta(z)=1-3z+z^2. The generating function is

    D_2(z)=sum_(n>=3) d_n z^n
          =z^3(1-z^2)/Delta(z).                     (5)

Its first coefficients at n=3,4,5,6,7,8 are 1,3,7,18,47,123. The unordered-pair convention and the triangle orbit correction are essential.

## 3. A complete canonical grammar for start count two

Let b_(n,2) count rooted strong tournaments (T,r) with s_r(T)=2. Such a graph has no order below four.

**Theorem 2.** Every such rooted graph has the following uniquely recoverable description. Delete r. The initial strong component U has order m>=3; all later components are singleton vertices in a transitive tail W. Exactly one of these two cases occurs:

1. r beats exactly one vertex a of U, and (U,a) has start count two.
2. r beats exactly two vertices a,b of U, both of which have start count one. Thus U is a C=1 graph with an unordered marked pair counted by d_m.

Every other U vertex beats r, all U vertices beat W, and if W is nonempty its last vertex beats r. All other incidences between r and W are arbitrary.

Conversely every description of either kind produces a strong graph with start count two.

**Proof.** For the SCC decomposition of T-r, the general root-start identity is

    s_r(T)=q*sum_(a in U,r->a) s_a(U),

where q is the product of H over all later components. Every nontrivial strong component has H>=3. Thus s_r=2 forces q=1 and makes all later components singletons. Every summand in the initial component is positive, so the sum two has exactly the two stated possibilities. A singleton U cannot supply sum two, hence m>=3. In the second case U is strong and has a start-one vertex, so C(U)=1.

Strongness forces the final tail vertex to beat r. For the converse, r has at least one in-neighbor and one out-neighbor in U, because neither allowed entry set can be all of U. If the tail is empty this makes U+r strong. If the tail is nonempty, U->W->r->U gives strongness. The root-start identity then gives exactly two. The ordered SCCs, marked entry set and ordered tail are intrinsic to (T,r), so there is no ambiguity or multiplicity from construction histories. QED.

This is a classification of rooted start count two, not a claim that every C=2 tournament admits such a root. In particular the regular five-tournament has three Hamilton paths from every root and is outside this rooted class.

## 4. Exact all-order enumeration

If d=n-m is the number of added vertices including the new root, the tail choices have multiplicity

    k_1=1,       k_d=2^(d-2) for d>=2.

Theorem 2 therefore gives the exact recurrence

    b_(n,2)=sum_(m=3)^(n-1) k_(n-m) * (b_(m,2)+d_m).  (6)

The two marked-entry cases are disjoint. In particular the pair case uses unordered pairs, not twice their number. The generating function K(z)=z(1-z)/(1-2z), together with (5), yields

    B_2(z)=K(z)(B_2(z)+D_2(z)),

    B_2(z)=z^4(1-z)(1-z^2)/(1-3z+z^2)^2.             (7)

At orders 4 through 10 the counts are

    1,5,18,60,191,589,1775.

Equivalently the recurrence for n>=8 is

    b_(n,2)=6b_(n-1,2)-11b_(n-2,2)
            +6b_(n-3,2)-b_(n-4,2),

with initial counts 1,5,18,60 at orders 4 through 7.

## 5. Verification and remaining boundary

The verifier builds the complete C=1 grammar in unique-cycle coordinates, checks the consecutive-block and outdegree statement on actual graphs, constructs the unordered-pair bank, and builds the canonical start-two grammar through order ten. A rooted start-two graph is canonically labelled by taking the smaller adjacency encoding in its two actual rooted Hamilton-path orders. An isomorphism preserving the root permutes those two paths, so this is an exact rooted-isomorphism key rather than a heuristic.

An independently generated complete labelled universe through order six checks the start-two inventory. The recurrence is also checked far beyond those physical orders. The certificate records all ranges. Finite controls do not replace the induction and canonical-decomposition proofs.

The full C=1 structure is now settled within these packets. The remaining bridge-free C=2 problem is restricted to 2-strong graphs of order at most73 by `C2_CUT_REDUCTION.md`; this new rooted grammar does not silently resolve that residual problem.
