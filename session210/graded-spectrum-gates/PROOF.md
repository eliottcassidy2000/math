# A locally defined cycle spectrum and finite closure-response reconstruction

Status: PROVED elementary statements; FINITE-EXACT tests in the explicitly stated universes. No literature-priority claim. Independent review is supplied separately. The unrestricted Chan–Pak Hamiltonian-path spectrum remains OPEN.

## 1. Inheritance, objects, and the connection being made

For a finite tournament T, H(T) counts directed Hamiltonian vertex words. Set H(empty)=1 only for the algebra of covers. Write HC(T) for directed Hamiltonian cycles **modulo cyclic rotation**, not rooted cycle words. Let S_n={H(T): |T|=n}, n>=1.

The closest inherited mechanisms are:

- The standard path-to-cycle closure, including the repository's explicit cone over C3 in THM-052, `01-canon/theorems/THM-052-circulant-scalar-transfer.md`, and the sole-exit path/cycle interface in the audited Session189 packet `05-knowledge/results/session189_response_frontier_20261104/cutarc_incoming/PROOF.md`.
- Ordered path-cover polynomials and their size-weighted odd-cycle form: THM-078, `01-canon/theorems/THM-078-uT-size-weighted-independence.md`. The proof below supplies an independent elementary parity argument, without requiring that file's root assertions or restrictions to odd order.
- Session207's inverse upper-ray and strong-component min-plus reduction, `session207/inverse-size-spectrum/README.md`.
- THM-1370 **with slug** `01-canon/theorems/THM-1370-h-spectrum-omits-7-21-all-n.md` proves the omissions 7 and21. Another unrelated file has the same theorem number.

The canonical hostile is scalar-context blindness: H(A)=H(B) need not imply H(context(A))=H(context(B)). The least-used sidecar here is the complete ordered path-cover polynomial. The concept board is: local degree constraints; path/cycle grading; marked versus unmarked reconstruction; cover polynomials; nested-context observations; dominance gluing.

The map will change a path-counting problem into a cycle-counting problem on a class defined only by vertex degrees. It preserves counts with a two-vertex size shift. At the object level it loses nothing once the canonical gate is identified (except the order-three automorphism boundary); its scalar observations retain the cover polynomial but still lose the tournament. These are explicit transports, not an assertion that a change of vocabulary solves the original existence problem.

## 2. Exact graded spectrum from a local degree condition

Call an arc a->b a **gate** when d+(a)=d-(b)=1. The name abbreviates this exact local condition; it imposes no path-count condition. For a nonempty tournament T define Gamma(T) on V(T) union {a,b} by

    a -> b -> T -> a.

The old arcs of T remain unchanged.

**Theorem 1.**

1. Gamma(T) is strong and has a gate a->b. Conversely, every tournament G of order N>=3 with a gate is precisely Gamma(G-{a,b}).
2. Deleting a,b from a Hamiltonian cycle is a bijection from the cycles of Gamma(T) to the Hamiltonian paths of T. Consequently HC(Gamma(T))=H(T).
3. A tournament of order N>=4 has at most one gate. Thus for n>=2, Gamma gives an unmarked, relabeling-compatible bijection between tournament isomorphism classes of order n and gated tournament classes of order n+2. Their automorphism groups are naturally isomorphic.
4. For every n>=1,

       {HC(G): |G|=n+2 and G has a gate}=S_n.

   If m(h)=min{n>=1:h in S_n}, with infinity for an empty set, the minimum gated order with h Hamiltonian cycles is exactly m(h)+2. Its set of possible orders is the upper ray N>=m(h)+2 when m(h) is finite.

**Proof.** The degree condition forces b to beat every old vertex and every old vertex to beat a, giving the stated converse. Every old vertex reaches b through a, and b reaches every old vertex, so G is strong. Every Hamiltonian cycle must use a->b. Starting that cycle at b gives b, a Hamiltonian word of T, a, followed by the closing arc a->b. This is an exact bijection and counts a cycle once.

For uniqueness, suppose G=Gamma(T) with |T|>=2. Vertex b has outdegree |T|, so cannot be another gate tail; a is already a gate tail with its unique head b. Any possible different tail x in T would have unique out-neighbor a, so a would have to have indegree one to be its head. But d-(a)=|T|>=2. There is no second gate. An isomorphism therefore fixes the gate's two roles, restricts to a core isomorphism, and conversely every core isomorphism extends. The same restriction-extension proves the automorphism statement. If |T|=1 the target is C3, which has three gates and automorphism group C3: the unmarked class correspondence still holds there, but automorphism preservation does not. A marked gate removes this exception.

The spectrum and minimum-order statements now follow in both directions, not just by a construction inclusion. Adjoining a universal source to T preserves H(T), so all larger gated orders occur once the minimum does. QED.

The class therefore has odd HC and never HC=7 or21, using the inherited tournament parity and omission theorems. These are conclusions for the local gate class, not for all strong tournaments. No new global lower bound or eventual threshold follows.

**Exact labeled count.** For N>=4 the number of gated tournaments on a fixed N-element set is

    N(N-1) * 2^binom(N-2,2).

Choose the unique ordered gate, then orient the core freely. At N=3 there are two such labeled tournaments, not six: each of the two oriented triangles has three gates. This is also a direct hostile to forgetting the marked/unmarked boundary.

## 3. The full simple-cycle deck identity and gluing

Let c_k(T) count simple directed k-cycles modulo rotation, and put

    C_T(z)=sum_(k>=3) c_k(T) z^k,
    A_T(z)=sum_(nonempty S subset V(T)) H(T[S]) z^|S|.

**Theorem 2.**

    C_Gamma(T)(z)=C_T(z)+z^2 A_T(z).                 (1)

More precisely the cycles using a or b use both, and for every nonempty subset S of the core their set is in bijection with the paths of T[S].

**Proof.** A cycle containing a must next use b, and one containing b must previously use a. The segment between b and a is exactly a directed path on a nonempty core subset. Conversely every such path closes uniquely. The other cycles lie wholly inside T. QED.

The subset version is lossless as a family of marked counts; the univariate polynomial (1) sums over all subsets of a given cardinality and forgets their labels. It must not be called the complete labeled deletion deck. Deleting core vertices commutes with Gamma as long as at least one remains. Deleting the whole core produces a transitive pair, outside the strong gate class.

For disjoint cores A and B let A->B denote their dominance join. Then

    HC(Gamma(A->B))=H(A)H(B),
    |Gamma(A->B)|=|Gamma(A)|+|Gamma(B)|-2.            (2)

This gives a physical multiplicative gluing law with grade N-2: take the two cores in dominance order and use one outer gate. It is associative on core isomorphism classes and not generally commutative. Formula (2) follows by concatenating the two paths; no scalar factorization is claimed to identify the core uniquely. The same grade connects directly to Session207's inherited additive order cost for strong factors. The retained sidecars multiply as well: P_(A->B)(q)=P_A(q)P_B(q), because each labeled bin concatenates its A-path and B-path; and 1+A_(A->B)(z)=(1+A_A(z))(1+A_B(z)), by choosing the induced subsets independently. Thus the graded construction is compatible with explicit combinatorial products, not merely equal scalar supports.

## 4. Ordered covers and their elementary parity

For q a nonnegative integer let P_T(q) be the number of ordered q-tuples of possibly empty directed paths whose vertex sets partition V(T). In particular P_T(1)=H(T), P_T(0)=0 for nonempty T, and P_empty(q)=1. This differs from counting unordered path covers.

**Lemma 3.** P_T(q) is a monic polynomial of degree n=|T|. For all q,

    P_T(-q)=(-1)^n P_T(q).                           (3)

For n>=1 it has constant term zero.

**Proof.** Work over Q in the commutative square-zero ring R=Q[x_v:v in V]/(x_v^2). Put

    F_T(x)=sum_(S subset V) H(T[S]) product_(v in S) x_v.

Then P_T(q)=[x_V]F_T(x)^q for every nonnegative integer q. Because F=1+N with N nilpotent, the finite binomial expansion defines a polynomial in q of degree at most n. The top-degree term comes only from the n singleton monomials: [x_V]N^n=n!, and binom(q,n) has leading coefficient1/n!. Thus the polynomial is monic of degree n. At q=0 the coefficient is zero for nonempty V.

We show F_T(x)F_T(-x)=1. At a nonempty subset S its coefficient is the signed count of ordered pairs of directed paths whose vertex sets partition S, with sign (-1) to the size of the second path. If both are nonempty, compare their first vertices. Move whichever first vertex beats the other to the front of the other path. This preserves both directed-path conditions. In the new pair, the moved vertex beats its original successor (if any), so applying the rule again moves it back. If one path is empty, move the first vertex of the nonempty path into it; the same argument gives the inverse. This fixed-point-free involution reverses the sign, proving the identity.

Hence F_T(x)^(-q)=F_T(-x)^q. Taking [x_V] proves (3), first at integers and then as polynomials. QED.

This parity is a known part of tournament path-cover theory. The point here is an independent proof that avoids importing an odd-order restriction or any real-root claim.

## 5. A universal local-closure operator

**Theorem 4.** For every tournament T, including the empty tournament if Gamma(empty) is interpreted as the transitive pair,

    P_Gamma(T)(q)=q[P_T(q+1)+qP_T(q)-P_T(q-1)].       (4)

Thus on polynomials the closure acts by

    (Lf)(q)=q^2 f(q)+q f(q+1)-q f(q-1).

**Proof.** For positive integer q partition the covers of Gamma(T) into three disjoint cases.

- a,b are in different bins. Then a is last in its bin and b first in its bin. There are q(q-1) choices of these bins, followed by any q-bin core cover: q(q-1)P_T(q).
- a,b are in the same bin with a before b. They are consecutive. Choose their bin in q ways and split its core portion into the path before a and the path after b. This is an arbitrary (q+1)-bin core cover: qP_T(q+1).
- a,b are in the same bin with b before a. The portion between b and a must be a nonempty core path. There are q choices of this bin and P_T(q)-P_T(q-1) core covers with that bin nonempty.

Adding gives (4). Polynomial equality follows from the infinitely many positive integer q. This proof includes the empty core, where the third case is empty and (4) gives q². QED.

For nonempty T, setting q=1 gives

    H(Gamma(T))=H(T)+P_T(2).                         (5)

The operator raises degree by two and preserves parity. It is injective on polynomial space because its top term is q²f; but that alone would not prove the observation theorem below.

## 6. Finite reconstruction from nested closure counts

For r>=0 define a_r(T)=H(Gamma^r(T)). The r=0 observation is the original count. These are counts of actual tournaments; their sizes are n+2r.

**Theorem 5.** Let n>=1 be known, and put d=floor((n-1)/2). The d observations

    a_0(T),...,a_(d-1)(T)

uniquely determine the entire polynomial P_T. For n=1,2 no observations are needed: the polynomials are q and q². Therefore, for two tournaments of the same order, the following are equivalent:

- Their P polynomials agree.
- All their nested-closure Hamiltonian counts agree.
- Their first d nested-closure Hamiltonian counts agree.

**Proof.** Write f=P_T. By Theorem4, a_r=(L^r f)(1). On evaluations at positive integers,

    (Lf)(k)=k²f(k)+k f(k+1)-k f(k-1).

After r steps the largest possible argument is r+1. There is exactly one sequence of shifts attaining it: every shift is upward. Its coefficient is 1*2*...*r=r!. Thus

    a_r=r! f(r+1)+a known integer linear combination of f(0),...,f(r). (6)

With f(0)=0, the first d responses recover f(1),...,f(d) successively. By Lemma3, f is monic of degree n, has parity n, and has zero constant coefficient. Its unknown coefficients are those at exponents1,3,...,n-2 if n is odd, or2,4,...,n-2 if n is even. There are exactly d of them. After subtracting the known leading term and dividing each evaluation by k or k² respectively, evaluations k=1,...,d give a Vandermonde system at the distinct nodes1²,...,d². It is invertible. The equivalences follow. QED.

This is exact reconstruction of an invariant from a specified family of contexts. It is not reconstruction of T, not an efficient algorithm for evaluating the counts, and not a minimum-query theorem for the discrete set of actually attainable polynomials. The largest queried tournament order is n+2(d-1) when d>=1. It is2n-3 for odd n>=3 and2n-4 for even n>=4. No claim of global spectrum completeness is made.

The local recurrence begins

    a_0=f(1),
    a_1=f(2)+f(1),
    a_2=2f(3)+5f(2)-f(1),

for nonempty cores. The factorial coefficient in (6), not triangularity in the monomial basis, is what makes finite observation sufficient.

### 6.1 The same observations reconstruct the descent distribution

For a permutation word pi of V(T), let f_T(pi) be the number of its adjacent pairs directed forward in T, and let m_f count words with f forward pairs. A binomial-transform identity is

    P_T(q)=sum_pi binom(q+f_T(pi),n)
          =sum_(f=0)^(n-1) m_f binom(q+f,n).         (7)

**Proof.** Concatenate the nonempty bins of a q-bin path cover to get a word pi. Every backward adjacency must be cut. If there are k nonempty bins, choose the other cuts among the f forward adjacencies in binom(f,n-k) ways, and their increasing bin labels in binom(q,k) ways. Vandermonde's identity gives sum_k binom(f,n-k)binom(q,k)=binom(q+f,n). QED.

The n polynomials binom(q+f,n), f=0,...,n-1, are linearly independent: their evaluations at q=1,...,n form a triangular matrix in reversed f order, with diagonal1. Thus P_T and the complete forward-adjacency distribution determine one another. The finite closure observations in Theorem5 therefore recover that distribution too. Word reversal gives m_f=m_(n-1-f), consistent with the parity used in the shorter reconstruction. This is the familiar Eulerian/binomial-transform mechanism, applied here with explicit conventions; no literature-priority claim is made.

At q=2, only the top two f values contribute. For n>=2,

    P_T(2)=(n+1)H(T)+m_(n-2),
    H(Gamma(T))=(n+2)H(T)+m_(n-2).                  (8)

The added term counts words with exactly one backward adjacency. It also has a direct geometric-free bijection: if the unique defect is between v_i and v_(i+1), send

    v_1,...,v_i | v_(i+1),...,v_n
      to v_(i+1),...,v_n,a,b,v_1,...,v_i.

The result is a Hamiltonian path of Gamma(T) whose endpoints do not close to a directed cycle. Conversely every nonclosing Hamiltonian path in Gamma(T) has nonempty core segments on both sides of its consecutive a,b, and recovers exactly this defective word. The other Hamiltonian paths are the n+2 rotations of its H(T) Hamiltonian cycles. At n=1 the defect count is zero and the same scalar formula holds with that explicit convention.

Identity (7) is the exact interface to any further geometric model whose face incidences are determined by forward-adjacency distributions. Such a model needs its own definitions and proof; no unprovided geometric equivalence is required here.

## 7. Sharp boundaries and explicit hostile controls

All rows below are outgoing-neighbor bitmasks, with vertices numbered from0. Bit j in row i means i->j.

**Same scalar, different closure.** At order5, the cores

    A=(10,16,3,6,13),
    B=(24,17,3,6,12)

both have H=15, but

    P_A(q)=q^5+8q^3+6q,       P_A(2)=108,
    P_B(q)=q^5+10q^3+4q,      P_B(2)=120.

Consequently H(Gamma(A))=123 and H(Gamma(B))=135. This is the first possible order for this obstruction: at orders1,2 the polynomial is fixed; at orders3,4 Lemma3 shows the sole free coefficient is determined by H. The explicit order5 controls and the all-order degree argument establish the minimal order, not merely a bounded search.

**The connecting arc matters.** The strong order5 tournament

    U=(18,8,3,5,14)

has outdegrees(2,1,2,2,3), indegrees(2,3,2,2,1), and exactly two Hamiltonian cycles. Thus it has a sole-exit vertex1 and a sole-entry vertex4, but the arc is 4->1, not1->4. There is no gate. Having both extreme-degree vertices does not imply the odd cycle-count law. This is the least possible order for an even positive HC with those conditions: the exhaustive target check covers all labeled tournaments through order4.

**The cover invariant still loses the object.** A source dominating C3 and its converse C3 dominating a sink are nonisomorphic (their score multisets are(1,1,1,3) and(0,2,2,2)), but both have P(q)=q^4+2q². Reversing every path gives cover-count invariance under converse, or dominance multiplication gives q(q³+2q). Every nested-closure scalar therefore agrees, while the cores and, by gate uniqueness, all their closures remain nonisomorphic.

**The empty-core boundary matters.** Gamma(empty) is transitive of order2. It has a degree-defined gate but is not strong and has no Hamiltonian cycle; H(empty)=1 cannot be used in Theorem1. It is included only in the algebraic operator formula, where it gives P=q² correctly.

## 8. Reproduction and scope

Run from this directory using only standard Python:

    python3 verify.py > /tmp/gate-result.json
    cmp RESULT.json /tmp/gate-result.json
    python3 -O verify.py > /tmp/gate-result-O.json
    cmp RESULT.json /tmp/gate-result-O.json

The primary verifier exhausts every labeled core of orders1..5 (1,099 cores), compares its path DP against direct permutations, checks all simple cycles of every closure, all induced-subset alternating-convolution identities, and the operator at q=1,2,3,4 (4,396 checks). It reconstructs every core's cover polynomial from the stated number of actual closure counts and compares against independent subset-convolution evaluation at q=0..5. Direct permutation distributions independently check (7) and the one-defect formula (8) for every core.

A separate exhaustive target pass over all labeled tournaments of orders3..6 recognizes gates from degrees, verifies labeled counts2,24,160,1920, and compares the cycle spectra to S_1,...,S_4. It also finds the unlinked-degree hostile. Monic parity-polynomial reconstruction is tested through degree20 on non-tournament controls, to exercise signs and the general algebra rather than only physically nonnegative examples. These finite tests check implementations; the infinite assertions follow from the proofs.

The equality of constrained spectra is a two-way structural correspondence, not a count-preserving encoding defined by H itself. Nevertheless it transfers the original problem rather than resolving it: no construction of a gate tournament with an arbitrary prescribed admissible count has been supplied. The new useful coordinate is that one simple local operation has a closed cover-polynomial action and finitely many repetitions recover precisely that stronger invariant.
