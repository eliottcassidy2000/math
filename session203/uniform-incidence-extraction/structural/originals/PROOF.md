# No fixed letter-local graft extracts an affine marked-edge count

Status: **PROVED CANDIDATE / INDEPENDENT AUDIT PENDING.** The argument is all-order within the exact construction class below. Eventual odd-spectrum coverage remains **OPEN**. This result strengthens a finite failed graft search to an arbitrary-size obstruction; it supplies no new spectrum threshold.

## 1. Scope and inherited source

For every integer u>=2, write u-1 as a binary word b_0...b_(m-1), with b_0=1. The inherited marked-edge tournament E(u) has vertices r,c,p_0,...,p_(m-1), with

- r->c and every p_i->r;
- p_i->p_j for i<j;
- c->p_i exactly when b_i=1.

Its forward edge-order sector has size u, while its total count is u+g(u)-(u mod2), where

    g(1)=1, g(2v)=g(v), g(2v-1)=g(v)+v-(v mod2).

The marked-edge construction and its interpretation are inherited from `session202/positive-overlap/EDGE_SECTOR.md`. The root-amalgamation identity used below is inherited from `session197/path-polynomials/ROOTED_ORDERED_AMALGAMATION.md`. No priority claim is made for either ingredient.

A **fixed letter-local graft** chooses an arbitrary fixed finite set Z of new vertices, arbitrary fixed edges within Z and between Z and {r,c}, and, for each z in Z, two fixed bits eta_z(0),eta_z(1). For every u it retains ALL of E(u) and puts z->p_i exactly when eta_z(b_i)=1. Thus the incidences of a new vertex may be uniformly incoming, uniformly outgoing, copied from c, or opposite to c. The number of new vertices is not bounded by the theorem. Call the resulting actual tournament T(u).

The graft is fixed across targets. It may not edit old edges, delete old vertices, use the position of a private vertex as another input, change its own size with u, or select different grafts for different u. Those are genuine exits from this result. In particular the theorem does not forbid a target-dependent rebuilding scheme or a finite menu whose pieces work on different subsets.

**Theorem.** For no fixed letter-local graft do there exist real constants a,b and an integer U0 such that

    H(T(u))=a u+b for every integer u>=U0.

Consequently this class cannot realize either desired extractor 6u-1 or 6u+1 on a full tail. This holds for every fixed finite graft size, not just the previously screened sizes through three.

Inheritance board: free edge-order sector; transitive-module pumping; physical path insertion; root amalgamation; a fixed-bit mixed difference. The closest mechanism is the exact sector construction. The hostile near miss is replacing freely controlled sector size by freely controlled total count. The least-used sidecar is the number of times a Hamilton path enters a transitive module.

## 2. Three module visits force base-three growth

Let a tournament contain a transitive module M of size t: all its vertices have the same incidences to every vertex outside M. Partition the outside vertices into A and B, where M->A and B->M.

**Lemma 1 (two cross edges).** If there are two vertex-disjoint arcs a_1->b_1 and a_2->b_2 with a_i in A and b_i in B, then

    H(T)>=3^t-3*2^t+3, for t>=1.                         (1)

The useful assertion is for t>=3.

**Proof.** First replace M by three ordered vertices x_1,x_2,x_3. There is a directed path

    x_1,a_1,b_1,x_2,a_2,b_2,x_3.

Insert every remaining outside vertex into this path, retaining the relative order of vertices already present. Such insertion is always possible: insert at an endpoint when legal, and otherwise use a consecutive pair y_i,y_(i+1) with y_i->v->y_(i+1). The resulting Hamilton path still has the three x vertices separated by outside vertices.

Replace x_1,x_2,x_3 by three nonempty increasing subsequences partitioning the t original module vertices. Every assignment to the three labelled subsequences is legal because M is a transitive module. The number of onto assignments is 3^t-3*2^t+3. Different assignments give different Hamilton paths, proving (1). QED.

If a fixed binary context contains a variable block b^t, all vertices of that block form a transitive module in T(u). Its index has form u=A*2^t+B for constants A,B determined by the fixed context, with A>0 in the uses below. A hypothetical affine formula therefore grows as O(2^t), contradicting (1). No spectral interpolation or cancellation of signed coefficients is involved.

## 3. The first pump removes every bit-dependent attachment

Suppose a new vertex z has eta_z(1)=0 and eta_z(0)=1. Use binary word

    1,0,1^t.

The final block is the module M. The fixed arc r->c is an A-to-B arc because M->r and c->M. The earlier zero-bit vertex y satisfies y->M, while M->z and z->y. Thus z->y is a second disjoint A-to-B arc. Lemma 1 contradicts an affine formula on the tail.

If instead eta_z(1)=1 and eta_z(0)=0, use word

    1^t,0.

Again r->c is one cross edge. For the final zero-bit vertex y, M->y, y->z and z->M, so y->z is the second disjoint cross edge. This is the same contradiction.

Both words are valid binary expansions with leading one. Their indices are respectively 3*2^t and 2^(t+1)-1. Consequently every new vertex in any hypothetical affine graft must be uniform toward all private vertices.

Partition the fixed new vertices into B, which beat every private vertex, and A, which lose to every private vertex. These labels agree with the module partition when the private module is any monochromatic block.

## 4. The second pump forces an ordered deletion graph

Use word 1,0^t, and let y be its first private vertex. The zero block M satisfies

    M->c->y->M.

Thus c->y is an A-to-B cross edge for that module. If there were an arc from a fixed new vertex in A to a fixed new vertex in B, that would be a disjoint second cross edge. The same is true of any arc r->z with z in B, since M->r. Lemma 1 excludes both possibilities.

Therefore all B vertices beat all A vertices and also beat r. Put S=A union {r}, with its fixed induced tournament. After deleting c, the entire graft has the ordered decomposition

    T(u)-c = B -> P -> S,                            (2)

where P is the private chain. The sets B or A may be empty, but S is always nonempty. All c incidences with B and S remain arbitrary. This deduction uses only existing arcs and literal module witnesses.

## 5. A positive surviving g coefficient

Let L denote the rooted tournament induced by B union {c}, let W(u) be c union P, and let R denote the rooted tournament induced by S union {c}, all rooted at c. Formula (2) identifies T(u) with the actual ordered root amalgamation

    L * W(u) * R.

For a rooted tournament X use h_X=H(X), d_X=H(X-c), e_X and s_X for Hamilton paths ending and starting at c, and R_X,L_X for labelled two-path covers whose designated nonempty path starts and ends at c. The second path may be empty. Put H(empty)=1. To avoid confusing a tournament name with a cover coordinate, write

    d_L=H(B), d_R=H(S), x=R_L-d_L, y=L_R-d_R.

The unit source W(u) has

    h_W=g(u), d_W=1, s_W=1, e_W=u mod2,
    R_W=u, L_W=ell(b), P_P(2)=2^m,
    ell(b)=1+sum_(i=0)^(m-1) (1-b_i) 2^i.             (3)

The formula for ell follows by splitting the private vertices into two labelled increasing paths, with the one ending at c nonempty beyond c precisely when its last private vertex loses to c.

Applying the exact root-amalgamation law twice gives

    H(T(u)) = d_L d_R g(u)
              + d_L y u + d_R x ell(b) + xy 2^m
              + d_L(s_R-d_R)(u mod2) + K,            (4)

where the fixed constant is

    K=d_R h_L-d_R x-d_L d_R-d_L y+d_L(h_R-s_R).

This formula includes B empty: then L is the rooted singleton, d_L=h_L=1 and x=0. R is nontrivial because it contains r. In particular the coefficient

    kappa=d_L d_R=H(B)H(S)>=1                       (5)

is strictly positive. It has not been chosen independently of a graph.

For completeness, the first amalgamation X=L*W has

    h_X=x(ell-1)+h_L+d_L(g-1),
    e_X=d_L(u mod2), R_X=x2^m+d_Lu, d_X=d_L.

Substitute these into

    H(X*R)=(R_X-d_X)(L_R-d_R)
             +(h_X-e_X)d_R+(h_R-s_R)d_X+e_Xs_R

to obtain (4) directly.

## 6. A fixed-length mixed difference finishes the contradiction

At any fixed private length m>=3, all terms on the right of (4) except g(u) are affine functions of the individual bits b_i. Here u=1+sum b_i2^(m-1-i), ell is given in (3), u mod2=1-b_(m-1), and 2^m is constant.

The g term has a nonzero adjacent-bit mixed difference. One direct path decomposition of W gives

    g(u)=b_0+(1-b_(m-1))
          +sum_(i=0)^(m-2) (1-b_i)b_(i+1)
          +sum_(0<=j<i<m) b_j(1-b_i)2^(i-j-1).       (6)

Paths in the first sum keep the chain order across c. Paths in the second sum have their last pre-c private vertex later than their first post-c vertex; all vertices strictly between those endpoints independently choose one side. These cases are disjoint and exhaustive.

The coefficient of b_(m-2)b_(m-1) in (6) is -2: one contribution from each sum. Fix b_0=1 and every other bit except these final two. Taking the mixed difference of (4) over the four possibilities 00,01,10,11 therefore gives exactly

    -2 H(B)H(S),                                   (7)

which is nonzero. An affine function au+b has mixed difference zero, because u is affine in those bits.

Choosing m sufficiently large puts all four indices beyond any prescribed U0. This contradicts the alleged affine tail and proves the theorem.

A useful explicit quadrilateral uses all intervening bits zero. For every k>=2, every graft surviving the pumping restrictions obeys

    H(T(2^k+4))-H(T(2^k+3))
      -H(T(2^k+2))+H(T(2^k+1)) = -2H(B)H(S).

This is an exact physical obstruction, not a finite asymptotic fit.

## 7. Scope, controls, and the surviving research direction

The proof classifies necessary structural features of an affine graft and then disproves affine behavior for those survivors. It does not classify the full scalar image of any graft. Failure of one fixed extractor is not absence of the target count.

The theorem leaves open grafts chosen from a target-dependent menu, a number of new vertices growing with the target, positional incidences, edits of old edges, replacing old vertices, and two-source C1/C2/C3 rebuilding. Any successful scalar extraction must exploit at least one such exit. It does not follow that every exit is sufficient.

`verify.py` enumerates every labelled graft on zero, one, or two new vertices in the declared class, checks each forbidden configuration's literal three-module-visit witness, and independently tests the exact quadrilateral on every survivor by subset Hamilton-path DP. It also checks (3) and (6) against literal graphs. The finite computation checks implementation and hostile boundaries; the theorem is proved above for arbitrary fixed finite graft size.
