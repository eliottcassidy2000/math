# Two fixed physical gates complete the start-two profile grammar

**PROOF CANDIDATE, pending independent audit.** This supplies an exact physical normal form and target-bounded decision procedure for the start-two fiber Gamma(B). It does not prove the staggered hypothesis F(8,137), all-base occupancy, or an eventual Hamilton-path threshold.

## 1. Inheritance and scope

For a nonempty rooted tournament (U,v), retain the nine physical coordinates

    (h,s,e,L,R,d,C,b,T),

where h=H(U), s and e count paths starting and ending at v, d=H(U-v), L and R count ordered two-path covers with a designated nonempty path ending and starting at v, C=B_U(2), b=B_(U-v)(2), and T counts one designated nonempty root-starting path plus two ordinary labeled paths. Ordinary paths may be empty; H(empty)=1. Put B=h+L-s.

The inherited O/Z/J operations append a global sink, append a vertex beating only v and keep v, or do the latter and mark the newborn. They have the audited nine-coordinate closure in dependencies/ROOTED_COVER_CLOSURE.md. Starting at the rooted edge v->w, their orbit contains exactly the nontrivial start-one sources, by dependencies/PHYSICAL_SOURCE_GRAMMAR.md, Theorem 1.

The other inherited ingredient is dependencies/START_TWO.md: its Theorem 1 says the start-one vertices of a strong unique-cycle tournament form a consecutive block of two or three vertices along its unique cycle; all but the last have outdegree one. The triangle has three start-one vertices, all outdegree one. Its Theorem 2 gives the complete canonical start-two decomposition by deleting the root. These are structural theorems, not inferences from a finite graph bank.

Closest proved mechanism: the complete unit-source grammar and the canonical two-entry start-two decomposition. Hostile case: unit-start B=137 is absent, while that says nothing about start two or shifted bases. Corrected near miss: a finite coupled state and a terminating inverse do not imply all-target existence. The retained sidecar is the complete k=2 cover block, not merely h and s.

## 2. Two literal gates

Input U has at least two vertices and root v. Let A=U-v. Preserve U and direct every vertex of A toward every new vertex below.

Gate K1 adds a and r, with

    r->a, a->v, r->v,

and marks r.

Gate K2 adds c,a,r, with

    c->v, a->c, v->a, r->a, r->v, c->r,

and marks r. Thus before r is added, its old induced graph is J(J(U)), rooted at a; r beats its first and third marked-chain vertices a,v and loses to c.

Both gates satisfy

    s(Ki(U))=2s(U).                                      (1)

For K1 an r-starting Hamilton word either starts r,a,v and then uses an old v-start word, or starts r,v, uses an old v-start word and ends a. In the second case A is nonempty; its final vertex points to a. These two types are disjoint and exhaustive because a has no outgoing neighbor other than v. For K2 the two types are r,a,c,v followed by an old v-start word, and r,v followed by an old v-start word and then a,c. Again they are exhaustive: if v is already used, a,c can only form the terminal pair.

The singleton boundary must be explicit. On U={v}, K1 gives a transitive three-tournament with s=1, not 2. K2 gives s=2 and agrees in rooted isomorphism type with K1 applied to a rooted edge. Our grammar starts at the rooted edge and never calls K1 on a singleton. All cover-coordinate formulas below remain valid even at a singleton; only this start rule needs the exception.

## 3. Why four coordinates at each cover level suffice

This section establishes an all-input finite-dimensional identity, so the finite basis calculation used afterward is not extrapolation.

For an integer k>=1 define

    X_k(U,v)=(C_k,L_(k-1),R_k,D_k),

where C_k=B_U(k), D_k=B_(U-v)(k), L_(k-1) is one designated nonempty root-ending path plus k-1 ordinary paths, and R_k is one designated nonempty root-starting path plus k ordinary paths. At k=1 this is (h,e,R,d); at k=2 it is (C,L,T,b).

**One-port linearity lemma.** Fix a tournament context G containing the shared vertex v and at least one other vertex r. Identify its v with the root of an arbitrary U, direct every U-v vertex toward every G-v vertex, and mark r. For each fixed k, its output X_k is a linear function of the input X_k, with coefficients depending only on G,r,k. The deletion coordinate is included in this statement.

**Proof.** Work in the commutative squarefree algebra on variables x_a for a in A=U-v; x_a^2=0. A path contributes the product of its vertices, with multiplicity when distinct paths use the same set. Put

- f=1 plus the enumerator of all nonempty directed paths in A;
- E=1 plus the enumerator of nonempty A-paths whose final vertex beats v;
- S=1 plus the enumerator of nonempty A-paths whose initial vertex is beaten by v.

Writing [A] for the coefficient of the product of all A variables, the input coordinates are exactly

    C_k = k[A] f^(k-1) E S,
    L_(k-1) = [A] f^(k-1) E,
    R_k = [A] f^k S,
    D_k = [A] f^k.                                  (2)

The first factor k chooses the ordinary bin containing v. Squarefreeness enforces disjointness; the first and last fragments adjacent to v concatenate through v. The other three formulas are the direct designated-path definitions, including empty ordinary bins.

Now enumerate a cover of the output and erase its A vertices, retaining the labeled bins and the G vertices. Since every A vertex points to G-v, an A-run can occur only before the first G vertex in its bin or immediately after v. There is at most one run of the latter kind. Every ordinary bin has an optional initial A-path, contributing f, except that a bin whose first G vertex is v has the optional v-ending prefix E. If a nonempty A-run follows v, it contributes S-1; its final vertex can be followed by any G-v vertex. Empty and nonempty choices of that exceptional run are counted separately, so no arc from v to its next G vertex is silently discarded.

For an output ordinary k-cover, this yields only the four products

    f^k, f^(k-1)E, f^k(S-1), f^(k-1)E(S-1).

The same list applies to an output root-ending designated cover: there are k total bins and its final r allows an optional A-prefix just like an ordinary bin. For an output root-starting designated cover, that bin has no initial A-prefix, leaving exactly k ordinary prefixes. If v is on it, the exceptional A-run may follow v; otherwise v is on an ordinary bin. These again give the same four products. Deleting the new mark r leaves an ordinary fixed context and is covered by the first case. The fixed context determines the finite multiplicity of every erased pattern; the remaining A-count is one of the displayed products. Expanding S-1 and applying (2) proves the lemma. This argument applies also when A is empty. QED.

A convenient basis of four actual rooted inputs, each with root zero, is

    (0), (2,0), (4,1,2), (6,0,2).

They are respectively a singleton, a source-rooted edge, a directed triangle, and a source-rooted transitive three-tournament. Their X_k rows are

    (k,1,1,1),
    (k^2,k-1,k+1,k),
    (k^3+2k,k^2-k+1,k^2+k+1,k^2),
    (k^3,(k-1)^2,(k+1)^2,k^2).

The determinant of this four-row matrix is -4k. It is therefore nonzero for every k>=1. Once one-port linearity is proved, exact counting on these four basis graphs determines the transfer on EVERY U. This is a finite basis identity, not an assertion that a large graph resembles one of the basis graphs.

## 4. The closed gate maps

Writing X_k as a column, the matrices at k=1 and k=2, respectively, are

    K1_1 = [[1,-1,3,-2], [0,0,1,-1],
            [2,-1,4,-1], [1,-1,1,0]],

    K1_2 = [[4,-6,10,-4], [1,-2,2,-1],
            [3,-2,6,-1], [2,-2,2,0]],

    K2_1 = [[3,-2,6,-2], [1,0,1,-1],
            [3,-1,8,-1], [1,0,2,0]],

    K2_2 = [[12,-12,28,-8], [3,-2,4,-2],
            [7,-4,18,-2], [4,-2,6,0]].              (3)

These are obtained by counting the four literal basis inputs and their literal gate outputs, then solving the nonsingular four-row system. `verify.py` reconstructs these matrices with exact rational elimination, independently of the advertised coefficients. Section 3 makes this determination universal.

Equivalently, in nine-coordinate order (h,s,e,L,R,d,C,b,T),

    K1:
    (h-e+3R-2d, 2s, R-d, C-2L+2T-b,
     2h-e+4R-d, h-e+R,
     4C-6L+10T-4b, 2C-2L+2T, 3C-2L+6T-b),

    K2:
    (3h-2e+6R-2d, 2s, h+R-d, 3C-2L+4T-2b,
     3h-e+8R-d, h+2R,
     12C-12L+28T-8b, 4C-2L+6T, 7C-4L+18T-2b). (4)

K1's singleton s exception from Section 2 still applies. Negative matrix entries do not license nonphysical inputs: every allowed grammar step is the literal graph gate. No cone-relaxation membership rule is used.

## 5. Complete normal form for start two

**Theorem.** Every rooted tournament with s=2 is isomorphic, as a rooted graph, to one obtained from the rooted edge by

    u, then one gate Ki, then w,

where u,w are finite O/Z/J words and i is 1 or 2. Conversely every such recipe has s=2. In particular this is a complete grammar for physical nine-coordinate start-two states.

**Proof.** First suppose the input is strong. Apply the inherited canonical start-two decomposition. If its new root has only one entry a into the initial deletion component, that component is itself a strong start-two source. The new root and the arbitrary transitive tail are J followed by O/Z letters, retaining the new root during the tail. Remove this suffix and continue on the smaller start-two source. Vertex order strictly drops, so eventually the decomposition uses two entry vertices, each with start count one, in a strong unique-cycle tournament V.

For a triangle, any unordered pair can be ordered consecutively on its directed cycle. For any larger V, the inherited unit-start block theorem says its entries either are consecutive within that two/three-vertex block, or are the first and third vertices of its three-vertex block.

In the consecutive case, let entries be a,v in that order. Vertex a has outdegree one, with sole out-neighbor v. Deleting a leaves a nontrivial rooted unit-start source U at v: every v-start path extends uniquely with a at its end, so s_v(U)=s_v(V)=1. The old V is literally J(U). Adding the new root beating a,v is exactly K1(U).

In the skipped case, write the block a,c,v. Both a and c have outdegree one, with sole out-neighbors c and v. Delete a,c. Every v-start path in the remaining nontrivial U extends uniquely with a,c at its end, and every v-start path in V has that terminal pair. Hence s_v(U)=1, and V is literally J(J(U)). The new two-entry root gives exactly K2(U). The remaining arbitrary tail is O/Z. U is nontrivial because this skipped case is used only beyond the triangle. The complete inherited unit-start source grammar supplies u. Restoring the removed one-entry layers and tails supplies w.

Now let the start-two source be weak. Its root belongs to its initial strong component V0. Component concatenation gives s=s_root(V0) times the product of Hamilton counts of the later components. Every Hamilton-count factor in the tail product is odd, so that product divides 2 and must be 1. Every later strong component is consequently a singleton, since a nontrivial strong tournament has at least three Hamilton paths. The initial component has start two and cannot be a singleton. It has the strong normal form just proved, and its singleton tail is a final O-run.

Conversely the initial prefix has start one, both allowed nontrivial gates double it, and every remaining O/Z/J letter preserves it. The maps are literal tournaments. QED.

No uniqueness of u,i,w or of a nine-coordinate representation is asserted. The grammar need not generate distinct graphs at distinct words.

## 6. Positive growth and a logarithmic complete target cap

For every nontrivial rooted input with positive s,

    R>=d+s, d>=e, h>=s, C>=2L, b>=L, T>=b.       (5)

The first inequality counts the disjoint root-singleton and spanning root-starting designated covers. The other inequalities are the inherited direct deletion and cover injections. Thus, using (4),

    B(K1)-B
      =3R-e-2d-s +(C-2L)+(b-L)+2(T-b)
      >=2s>0.                                   (6)

For K2 the difference is

    2h-2e+6R-2d +3C-3L+4T-2b-s>0,              (7)

as follows, for example, by writing its ordinary part as
2(h-s)+6(R-d-s)+7s+4d-2e and using (5); its cover part is nonnegative. The inherited O preserves B, while Z and J strictly increase B on these inputs. Hence target pruning B>M is safe at every prefix, even before the gate. Start counts remain positive.

Consider any recipe in the larger grammar allowing any number of K1/K2 gates, with terminal B<=M. Remove its terminal O-run, which preserves B and s. If no letters remain, its order is two. Otherwise let n be the input order of its final remaining operation.

- For final Z or J, inherited cover bounds give B_out>=2^(n-1), hence output order n+1<=floor(log2 M)+2.
- For final K1, L_out=C-2L+2T-b>=b>=2^(n-1), hence output order n+2<=floor(log2 M)+3.
- For final K2, L_out=3C-2L+4T-2b>=2C+2b>=2^(n+1), hence output order n+3<=floor(log2 M)+2.

Here B_out>=L_out because h_out>=s_out, and B_W(2)>=2^|W| follows by assigning each vertex independently to two labeled bins and taking one Hamilton path per bin. Therefore every represented pair with B<=M has a recipe of order at most

    max(2, 3+floor(log2 M)).                     (8)

Combined with Section 5, this is a complete all-order cap for start-two membership in Gamma(B), not merely for a guessed finite seed language. It gives a literal finite algorithm. It does not say a requested B belongs to the fiber.

## 7. Guarded compiler and exact controls

`compiler.py B` accepts a positive odd target B. It enumerates the physical edge grammar with at most one K1/K2 gate and the order cap (8). Stage zero is start one and stage one is start two. It prunes only B above the requested target, uses exact integer arithmetic, and stores the least known order for each complete nine-coordinate state. Because gates have different vertex costs, it uses order-priority search rather than incorrectly treating them as equal-cost letters.

If a state has the requested B and start two, it returns the literal rows and recipe. If the bounded universe is exhausted without a witness, it returns ABSENT_START_TWO_PROFILE. A resource limit returns INCOMPLETE and supplies no absence conclusion. It never labels an absent start-two profile as an absent tournament Hamilton count or as a counterexample to F.

`verify.py` directly counts the fixed basis and all labeled root-zero inputs through order 5, checks the universal gate matrices against their literal children, and handles the singleton exceptions. It checks the inherited O/Z/J maps on the same inputs where applicable, growth and the order-bound inequalities, compares a finite complete start-two adjacency inventory to the independently inherited canonical two-entry generator, and checks compiler positive, negative and resource-limit behavior. Normal and optimized runs must agree.

The all-input gate laws follow from Section 3 and exact nonsingular basis evaluation. The all-order completeness follows from Section 5. The finite computations test these mechanisms and the implementation; they do not establish all-base occupancy. F(8,137) and the unconditional odd-spectrum threshold remain OPEN.
