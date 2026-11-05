# All-order enumeration of the rooted start-one class

Session197, 2026-11-05. This is an exact enumeration of the class classified in `ROOTED.md`. It does not assert that every uniquely Hamiltonian tournament has a strong bridge. The unqualified C=1 census in that earlier packet remains finite evidence.

Let b_n count isomorphism classes of rooted strong tournaments (T,r) on n vertices with exactly one Hamiltonian path starting at r. A singleton is strong by convention. Let a_n count unrooted nontrivial strong tournaments with exactly one Hamiltonian cycle and at least one strong bridge. Write F_0=0,F_1=1,F_(j+2)=F_(j+1)+F_j.

**Theorem.**

    b_1=1, b_2=0, b_n=F_(2n-4) for n>=3;
    a_3=1, a_n=F_(2n-6) for n>=4.

In particular the counts 1,1,3,8,21,55,144 at orders 3 through 9 hold at every order by the corresponding formula for the bridge subclass. No inference that this subclass exhausts all C=1 tournaments is used.

## 1. A canonical decomposition, not an arbitrary construction history

For a nontrivial rooted strong (T,r) with start count one, `ROOTED.md` proves that T-r has a unique initial strong component U, followed by a transitive chain W of singleton components. There is exactly one vertex u in U beaten by r, and (U,u) again has start count one. Every other vertex of U beats r, every U vertex beats W, and if W is nonempty its last vertex beats r. All other arcs between r and W are arbitrary.

Let m=|U| and w=|W|=n-m-1. For w>=1 there are exactly 2^(w-1) choices. For w=0 there is one choice when m>=3. The case m=1,w=0 gives a nonstrong two-vertex tournament and is excluded. The case m=2 is absent because U must be strong.

The source SCC, its distinguished vertex u, and the ordered tail are intrinsically recoverable from (T,r). Consequently two such descriptions yield rooted-isomorphic outputs if and only if their rooted U's are isomorphic and the same tail bit word was chosen. This eliminates possible multiplicities from noncanonical gluing histories.

Thus for n>=3,

    b_n = b_(n-1) + sum_(m=1)^(n-2) 2^(n-m-2) b_m.       (1)

The term b_(n-1) is zero at n=3 and thereafter exactly counts the permitted w=0 case. With b_1=1,b_2=0, the ordinary generating function is

    B(z)=sum_(n>=1)b_n z^n = z + z^3/(1-3z+z^2).          (2)

For example, summing the canonical gadget sizes gives B=z+z(1-z)/(1-2z)*B-z^2; the subtraction removes the sole invalid two-vertex output. This rearranges to (2). The coefficients after b_3=1,b_4=3 satisfy b_n=3b_(n-1)-b_(n-2), exactly the even-index Fibonacci recurrence. This proves the first formula.

## 2. Distinguished outdegree-one vertices

Let c_n count isomorphism classes from b_n whose distinguished root has outdegree one. For n>=3 the canonical decomposition above permits precisely one bit word for each U: every vertex in W must beat the root. Therefore

    c_n = sum_(m=1)^(n-1) b_m.                            (3)

In particular c_3=1,c_4=2,c_5=5,c_6=13.

For any strong tournament of order at least three, an outdegree-one vertex v satisfies

    start_v(T)=C(T),

because every Hamiltonian path from v ends at a vertex pointing back to v, and cutting a Hamiltonian cycle at v gives the inverse map. Thus c_n also counts pairs consisting of a graph from the unrooted bridge class and an outdegree-one vertex, modulo rooted isomorphism. The equivalence between the bridge class and existence of a start-one root is proved in `ROOTED.md`.

A strong tournament of order n>=4 has at most two outdegree-one vertices. Indeed, three such vertices already have three internal outgoing arcs in total, so each must use its sole outgoing arc internally. They would form a directed triangle with no arc leaving it, contrary to strongness.

Every graph in the bridge class at n>=4 is rigid. Its nonempty set of one or two outdegree-one vertices is preserved by every automorphism. If there are two, their mutual arc prevents them from being exchanged, so each is fixed. An automorphism of a unique directed Hamiltonian cycle that fixes one vertex fixes every vertex. Hence the number of rooted classes counted by c_n is exactly the number of those vertices, one or two, per unrooted class.

## 3. The exact two-vertex correction

Let d_n count unrooted bridge-class tournaments having exactly two outdegree-one vertices, for n>=4. There is a bijection

    d_n = c_(n-1).                                      (4)

Given a strong C=1 tournament S with an outdegree-one marked vertex v, add u with u->v and every other vertex of S pointing to u. The result T is strong: S was strong and both directions between u and S are present. Vertices u and v have outdegree one, whereas every other vertex gains one outgoing arc and so has outdegree at least two. Moreover

    C(T)=start_u(T)=start_v(S)=1.

Conversely, in a strong T with exactly two outdegree-one vertices, name them u,v so that u->v. Let the unique out-neighbor of v be w, and set S=T-u. To prove S strong, fix a vertex x outside {u,v,w}. Such a vertex exists because n>=4. A simple path from w to x in T cannot use v or u internally: reaching v would next return to w, while reaching u would next pass through v and return to w. Strongness of T therefore supplies paths in S from w to every such x. Since x->v, appending that arc supplies a path from w to v as well. Every vertex of S other than v and w points to v, while v->w, so all vertices of S reach w. Hence S is strong. Every Hamiltonian path starting at u in T continues through v and then gives a Hamiltonian path starting at v in S, so start_v(S)=C(T)=1. Since v still has outdegree one, C(S)=1. The two operations are inverse on isomorphism classes; the direction u->v determines which vertex to delete.

Finally c_n=a_n+d_n for n>=4, by rigidity and the one-or-two alternative. Equations (3)-(4) imply

    a_n=c_n-c_(n-1)=b_(n-1).

At n=4 this gives a_4=b_3=1; for all n>=4 it gives the stated Fibonacci formula, with the directed triangle handled separately at n=3.

## 4. Verification and scope

`verify_fibonacci_class.py` constructs the canonical rooted class through order nine, checks rooted isomorphism uniqueness, exact start counts and cycle counts, counts marked outdegree-one graphs and underlying unrooted graphs, and checks the two-vertex bijection. A separately generated complete tournament universe through order seven verifies that no rooted start-one graph is missing from the canonical generation. The recurrence itself is checked to order one hundred.

Normal and optimized runs must agree. These are physical controls of the formula and its implementation. The proof above supplies the all-order conclusion. The separate question of a bridge-free C=1 tournament is unchanged.
