# The start-three grammar and complete C=3 enumeration

Session197. **All-order rooted classification and unlabelled enumeration.** This packet refines the complete C3 bridge theorem in `C3_COMPLETE.md`. It uses the audited C1 start-one block theorem and refined Fibonacci counts in `START_TWO.md`, and the cycle-preserving two-outdegree-one bijection in `BRIDGE_TWO.md`. No literature-priority claim is made.

Write s_r(T) for the number of Hamilton paths starting at r, and put

    Delta(z)=1-3z+z^2,
    K(z)=z+z^2/(1-2z)=z(1-z)/(1-2z).

K counts adding a new root and a transitive tail: an empty tail has one choice; a nonempty tail has its final vertex pointing to the root and arbitrary other incidences.

The principal conclusions are

    B_3(z)=sum b_(n,3)z^n
          =z^5(8-17z+8z^2)/Delta(z)^2,               (1)

for rooted strong tournaments with s_r=3, and

    C_3(z)=sum c_n z^n
          =z^5(1+5z-17z^2+11z^3-z^4)/Delta(z)^2,   (2)

for all unlabelled strong tournaments with exactly three Hamilton cycles. Rooted and unrooted isomorphism classes are distinguished throughout. These are not labelled counts; nontrivial automorphisms occur.

## 1. A physical clone description of C1 roots starting two paths

Let U be strong with C(U)=1. Its start-one set S consists, for order at least four, of a consecutive block of two or three vertices on its unique cycle. All but the last have outdegree one. At order three the whole triangle is S. The refined class counts are

    a_(n,2)=F_(2n-8),  a_(n,3)=F_(2n-7),  n>=4.

**Lemma1.** A rooted strong C1 graph (T,r) with s_r=2 is obtained uniquely by choosing two cycle-consecutive start-one vertices a,b in a smaller C1 graph U, oriented a->b, and replacing a by a transitive pair r->a. The new root r beats exactly a,b and loses to every other old vertex.

**Proof.** Use the canonical start-two grammar. If deleting r leaves a nonempty transitive tail, every r-start path ends at its final vertex, which points to r. Hence all two paths close and C(T)=2, contrary to the assumption. The deletion U is therefore strong. A single entry would make r outdegree one and again C=s_r=2. Thus r beats two start-one vertices a,b in U.

There are exactly two r-start paths, one entering U at each selected vertex. The unique path starting at a ends at the predecessor of a on the unique cycle, and similarly for b. A path fails to close at r precisely when that predecessor is the other selected vertex. Therefore C(T)=1 precisely when the selected vertices are consecutive on the cycle. Both cannot fail to close, since U has at least three vertices.

By the C1 block theorem, the first vertex a of such a consecutive pair has outdegree one, with its sole neighbor b. Thus the described attachment is exactly transitive cloning of a. Deleting the marked root recovers U and the selected pair, so the description is canonical. QED.

If U has a two-vertex start-one block, its one consecutive pair produces an output with one outdegree-one vertex and hence a two-vertex start-one block. If U has a three-vertex block a,b,c, cloning a preserves two outdegree-one vertices, while cloning b raises a's outdegree and leaves only one. The triangle's three choices form one orbit and give the unique rooted output of order four.

Consequently, for n>=4, the number of C1 graphs with a marked start-two root, classified by the output's start-one set size, is

    |S|=2: F_(2n-8),
    |S|=3: F_(2n-9),                              (3)

where F_(-1)=1 is used only for the order-four boundary. The corresponding generating functions are z^5/Delta and z^4(1-2z)/Delta.

Let M_n count a C1 graph with one marked start-one vertex and one marked start-two vertex. Since C1 graphs of order at least four are rigid, multiplying the two rows of (3) by2 and3 gives

    M(z)=z^4(3-4z)/Delta.                           (4)

Let T_n count a C1 graph with an unordered triple of start-one vertices. The triangle contributes one orbit, while at larger orders the entire three-vertex block is the only possible triple. Hence

    T(z)=z^3+z^4(1-z)/Delta
        =z^3(1-2z)/Delta.                          (5)

## 2. The complete canonical start-three grammar

Let (T,r) be strong with s_r=3. Delete r and write the remaining SCCs in their unique dominance order, with initial component U. Then

    s_r(T)=q * sum_(a in U, r->a) s_a(U),           (6)

where q is the product of Hamilton counts of all later components.

Every nontrivial strong component has H>=3. Moreover H=3 forces that component to be a directed triangle. For an elementary justification, Hamilton-path count is monotone under induced extension by the usual insertion map; every strong graph of order at least four contains the strong four-vertex graph, which has H=5. The strong-four existence argument is in `C3_COMPLETE.md`.

Equation(6) gives exactly two disjoint types.

**Type I: entry sum three, singleton tail.** All components after U are singleton vertices. The entry set inside U is one of:

- one vertex starting three Hamilton paths;
- a start-one vertex and a start-two vertex in a C1 graph;
- three start-one vertices in a C1 graph.

The first case is recursively a rooted start-three input. The latter cases are counted by M and T. Every other U vertex points to r. If the tail is nonempty its last vertex points to r and the other bits are arbitrary. There is exactly one invalid empty-tail case: U is the triangle and r beats all three vertices, making r a source. It would have order four and is excluded.

**Type II: entry sum one, one triangle in the tail.** U is a singleton or a strong rooted start-one graph, and r beats only its marked start-one vertex inside U. The later SCCs consist of exactly one directed triangle and any number of singleton components, in their ordered positions. Incidences from r to these components are arbitrary except that the terminal component must have a vertex pointing to r.

Conversely every allowed description is strong and has start count three. The initial entry is nonempty; a nonempty tail is reachable from U and reaches r through its terminal component. With an empty tail, the only case lacking an inneighbor in U is the explicitly excluded whole triangle. The SCCs, entry marks and incidence pattern are recovered by deleting r, so the grammar is unique on rooted isomorphism classes.

## 3. Rooted enumeration, including triangle symmetries

Let B_1 count rooted strong start-one graphs, with the singleton included. The already proved formula is

    B_1(z)=z+z^3/Delta
          =z(1-z)(1-2z)/Delta.                     (7)

A Type-II tail of total order k>=3 contains i singleton components before its triangle and j after, with i+j=k-3.

- If j=0, the triangle is terminal. Its root-incidence subsets have three rotation orbits: zero, one or two vertices beaten by r. The all-three subset is forbidden. The preceding singleton bits give3*2^i choices.
- If j>=1, all four triangle-subset orbits are allowed. The final singleton bit is fixed toward r; the other singleton bits give4*2^(i+j-1) choices.

Therefore the exact tail count is (2k-3)2^(k-3), with generating function

    G(z)=z^3(3-2z)/(1-2z)^2.                       (8)

These are orbit counts, not counts of three separately labelled triangle vertices. The triangle remains an intrinsic SCC, so no extra identification exchanges it with U or a singleton position.

The canonical grammar gives

    B_3=K(B_3+M+T)-z^4+z B_1 G.                   (9)

Substitution of (4),(5),(7),(8), followed by 1-K=Delta/(1-2z), yields (1). Its initial coefficients at orders5 through12 are

    8,31,106,343,1070,3252,9694,28469.

## 4. Start-three roots without another indegree-one vertex

Call a rooted start-three graph exceptional here if it has no indegree-one vertex other than possibly its root. This is a local definition for this section; it does not mean a bridge-free graph. Let E_3(z) count these rooted classes.

**Lemma2.**

    E_3(z)=z^5(2-z)/Delta.                         (10)

**Proof.** Inspect the canonical grammar, remembering that every old entry vertex gains one incoming arc from the new root, and every non-entry vertex keeps its indegree.

For the one-entry Type-I branch, the input has order at least five. No tail vertex can have indegree one because it receives all arcs from U. The new root has at least four inneighbors in U. Thus the output is exceptional exactly when the marked input is exceptional. This contributes K E_3.

For the mixed-pair Type-I branch, all indegree-one vertices of U must lie in the selected pair. The C1 lemma from `BRIDGE_TWO.md` gives an indegree-one vertex outside the whole start-one set. It must therefore be the selected start-two vertex. By Lemma1 every such vertex has outdegree two; its indegree is |U|-3. Hence |U| must be four. In the unique strong four-vertex graph, there is exactly one valid mixed pair: its start-two vertex and the other indegree-one vertex, which is start-one. This contributes z^4 K.

For the triple Type-I branch, the same C1 lemma excludes every input of order at least four, because an indegree-one vertex lies outside all start-one entries. The triangle input does qualify, with a nonempty tail. All old triangle vertices then have indegree two, all tail vertices have indegree at least three, and the new root is the only possible degree-one exception. This contributes z^3(K-z).

Type II never qualifies. A singleton U gains indegree one from r at its unique vertex. A nontrivial rooted start-one U has an indegree-one vertex different from its marked root, which persists in the output; this holds also for the triangle. That vertex is different from the new root.

Thus E_3=K E_3+z^4K+z^3(K-z). Solving gives (10). QED.

The coefficients for n>=5 are F_(2n-7), beginning2,5,13,34,89,233. This class includes the regular five-tournament marked at a root; it is not restricted to C3 inputs.

## 5. From rooted counts to every unlabelled C3 graph

Let A_n count strong C3 graphs having an outdegree-one vertex, and let P_n count such a graph with one outdegree-one vertex marked. Deleting the marked vertex u leaves an initial component with its sole outneighbor a, followed by an ordered SCC tail. Since s_u=C=3, the possibilities are:

1. (U,a) has start count three and the tail is transitive;
2. (U,a) has start count one and the tail has one triangle plus singleton components.

There are no free incidence bits here: u loses to every vertex except a. For a triangle tail of total size k>=3, there are exactly k-2 possible positions of the triangle. Consequently

    P(z)=z B_3/(1-z)+z^4 B_1/(1-z)^2.              (11)

A strong graph of order at least four has one or two outdegree-one vertices if it has any. Every automorphism fixes these vertices individually, because their mutual arc forbids exchange when there are two. Thus marking them contributes one or two distinct rooted classes, even when the rest of the graph has nontrivial automorphisms.

The audited two-outdegree-one bijection adds u beating only a marked outdegree-one vertex v, and conversely deletes the first vertex of the two-vertex degree-one chain. It preserves strongness and any fixed positive cycle count. Therefore the number of C3 graphs with two outdegree-one vertices at order n equals P_(n-1). It follows that

    A(z)=(1-z)P(z)
        =z B_3+z^4 B_1/(1-z)
        =z^5(1-z)(1+4z-6z^2)/Delta^2.             (12)

Now count C3 graphs with an outdegree-one vertex but no indegree-one vertex. Their marked outdegree-one construction cannot be the triangle-tail case: a singleton initial U supplies another indegree-one vertex, and a nontrivial rooted start-one U supplies one away from a. In the start-three case, absence of indegree-one vertices is exactly the exceptional-input condition of Section4. Hence their marked generating function is z E_3/(1-z).

The same two-outdegree-one bijection preserves the absence of indegree-one vertices. In the forward direction the newborn has indegree at least four and only the marked old vertex gains an incoming arc. In the inverse direction only that marked old vertex loses an incoming arc; it still has outdegree one in a strong C3 graph of order at least five, so its indegree is at least three. No indegree-one vertex is created or lost elsewhere. The unmarked one-sided class therefore has generating function

    O(z)=z E_3(z).                                  (13)

By `C3_COMPLETE.md`, every strong C3 graph has an outdegree-one or an indegree-one vertex. The classes counted by A already include those with both. The remaining indegree-only class is the converse of the outdegree-only class counted by O. Hence the complete unlabelled count is

    C_3(z)=A(z)+O(z),                               (14)

which simplifies to (2). Its coefficients at orders5 through12 are

    1,11,38,124,390,1193,3574,10537.

The independent complete low-cycle census agrees at every checked order, but the recurrence and generating function follow from the canonical grammars and bijections above, not from fitting those numbers.

## 6. Verification and scope

`verify_start_three.py` pins the frozen `verify_start_two.py` helper by SHA256. It physically generates the C1 mixed marks and both start-three branches, canonicalizes rooted start-three graphs by their three actual rooted Hamilton paths, checks the exceptional subclass, and constructs the outdegree-one C3 class and its converse. It verifies both rational identities by a separate coefficient recurrence through order100 and compares the root-three grammar with all33,867 labelled tournaments through order six.

Normal and `-O` runs reproduce the same `start_three_certificate.json`. The independent C3 census separately tests the resulting unlabelled C3 coefficients through order12. The triangle automorphism corrections are essential: simply counting eligible vertices instead of rooted isomorphism classes gives wrong start-three totals.

Together with the earlier C1 and C2 packets, this gives exact unlabelled strong-tournament enumeration at every order for each cycle count1,2,3. It does not enumerate arbitrary cycle counts or resolve the unrestricted Hamiltonian-path spectrum `{7,21}`.
