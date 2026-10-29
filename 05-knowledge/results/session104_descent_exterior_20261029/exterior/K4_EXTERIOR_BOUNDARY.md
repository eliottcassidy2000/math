# K4, bivectors, and the actual 3/4/6 grammar: a real bridge and a decisive obstruction

**Results first.** The six oriented edge coordinates of K4 are naturally the basis of `Λ² Q⁴`. Their boundary gives an exact sequence with three-dimensional cycle kernel and three-dimensional potential-difference target. However, the actual six-state legality/singularity machine is **not** an exterior square of any four-dimensional letter representation, even after an arbitrary invertible basis change. Its `ee` transition has rank two, impossible for an exterior square in dimension four. Thus the grammar's 3/4/6 numbers are not an automaton-level tetrahedral exterior identity.

The question and cross-versus-wedge direction are user-supplied. Grammar, safe transfer, and minimal state counts are inherited from Session103 README and `languages/ACTUAL_GRAMMAR_REFINEMENTS.md`, ultimately the pinned `operation_closed_synthesis_20261029` source at `d3132eafbfa9c43c6d9f1110566b1e9c3a3286e1`. Session103's audited bridge/triangularity proofs supply the all-order singularity interpretation. The rank obstruction, precise K4 contract, and synthetic lift below are this packet's elementary deductions, not claims of historical novelty for exterior algebra. Relevant prior K4 work is Session87 `K4_BRIDGE.md`, inspected locally at `session87-four-vertices/connections/K4_BRIDGE.md`.

## 1. Exact six-state machine and universal obstruction

Use states `A,B0,B1,FA,FB,X`. The first three are safe, with B1 meaning previous b or e; F states have seen a forbidden pair and retain legality type; X is illegal. Start at B0. Outputs are safe, singular, illegal respectively.

|state|a|b|c|d|e|
|---|---|---|---|---|---|
|A|A|B1|X|X|X|
|B0|X|X|A|B0|B1|
|B1|X|X|A|B0|FB|
|FA|FA|FB|X|X|X|
|FB|X|X|FA|FB|FB|
|X|X|X|X|X|X|

All states are reachable: empty word, c, e, eec, ee, a reach B0,A,B1,FA,FB,X. Safe states are separated by one-letter continuations: A versus either B state by a, B0 versus B1 by e. FA versus FB are separated by a; other output groups already differ. This proves six-state Moore minimality.

Linearize a state map by sending its basis vector to the target basis vector (columns are sources). A deterministic linearized map has rank equal to the number of distinct image states. Every letter has rank three. But `ee` sends A,FA,X to X and B0,B1,FB to FB; its rank is exactly two.

**Theorem.** There do not exist four-dimensional rational letter operators L_a,...,L_e and an invertible S with `M_l = S (Λ² L_l) S⁻¹` for every letter l.

**Proof.** Exterior powers are functorial, so `M_e² = S Λ²(L_e²) S⁻¹`. For any linear operator T, `image Λ²T = Λ²(image T)`, giving `rank Λ²T = binomial(rank T,2)`. In dimension four the possible ranks are 0,1,3,6, never 2. Contradiction. This proof works over any field. ∎

The obstruction is stronger than a search over edge labels: it rules out permutation conjugacy, signed-permutation conjugacy, and arbitrary rational conjugacy. In particular the six states cannot carry a faithful transition-compatible wedge interpretation of the four acceptance states. Those four states are obtained by merging FA,FB,X into one rejecting sink. Each of their letter operators already has rank two, so its exterior square has rank one, not the six-state operator's rank three.

**Scope.** This does not rule out nonlinear encodings, larger representations, quotients of other constructions, partial operators with different illegal-word semantics, or a selected word-independent vector-space identification. A bijection of six named states with six edges always exists but does not intertwine transitions. Counting dimensions and deterministic state cardinalities are different invariants.

## 2. The genuine tetrahedral exact sequence

Let V=Q⁴ with vertex basis e0,...,e3, and ε(ei)=1. Identify the edge i<j with ei∧ej. Define

    ∂(u∧v)=ε(u)v−ε(v)u.

Thus `∂(ei∧ej)=ej−ei`, the ordinary oriented incidence boundary. Put H=ker ε. Then

    0 → Λ²H → Λ²V --∂--> H → 0

is exact. To see this, choose t with ε(t)=1. Then `V=Qt⊕H` and

    Λ²V = (t∧H) ⊕ Λ²H,
    ∂(t∧h)=h,       ∂(Λ²H)=0.

Both summands have dimension three. With the standard Euclidean metric the canonical orthogonal choice is t=(e0+e1+e2+e3)/4. The first summand is the cut/gradient space `image ∂ᵗ`; the second is the cycle space. The graph Laplacian is `∂∂ᵗ=4I−J`, hence 4I on H. A triangle cycle e01+e12−e02 is an explicit kernel vector.

This is a **quotient/projection principle**, not a claim that the grammar's counting quotient is this incidence map. No transition-compatible map from grammar states to these vertex potentials has been supplied; §1 excludes the most direct exterior interpretation.

Session87's binary boundary-plus-parity model is related by reduction modulo two, not by identical decomposition. There its cut and cycle spaces intersect in a dimension-two bicycle space; over Q they are orthogonal with zero intersection. In characteristic two, division by four and the Euclidean direct-sum argument fail. Session87 also explicitly distinguished its restricted nine-port decoder from tournament classification. Its word “wedge scaffold” denotes adjacent-edge pairings, not automatically an exterior product functor.

## 3. Cross product, wedge product, and two different 3+3 splittings

The wedge `u∧v` exists canonically in every dimension. In three dimensions a metric and orientation identify Λ²H with H through the Hodge star, giving the familiar cross product. In four dimensions Λ²V has dimension six: no vector-valued cross product arises from this identification.

With Euclidean metric and orientation e0e1e2e3, four-dimensional Hodge star satisfies `*²=1` on Λ²V. It sends

    e01 ↔ e23,    e02 ↔ −e13,    e03 ↔ e12.

Its +1 and −1 eigenspaces each have dimension three. They are **not** the cut and cycle spaces: star exchanges cut and cycle, and its eigenvectors are their diagonal mixtures. Metric and orientation are part of the definition; reversing orientation swaps the eigenspaces. An arbitrary grammar matrix is not thereby a Hodge-compatible operator. In particular none of these dimensional splittings repairs §1.

A differential-geometric “bivector field” would require a base space and a varying bivector at its points. The six machine states alone do not provide that structure. There is nevertheless a precise growing bivector carrier for the native observer, described next.

### A genuine tournament wedge carrier: the native observer’s antisymmetric sector

There is an actual tournament-level six-dimensional wedge carrier that must not be confused with the six legality states. At n=4 the native ordered-two-path observer M2 has 12 ordered-pair start ports and 12 ordered-pair end ports. Let R swap the two run labels, so `R e_(i,j)=e_(j,i)`. Swapping the two paths is a bijection of the counted path covers, hence

    M2_( (i,j),(k,l) ) = M2_( (j,i),(l,k) ),
    M2 R = R M2.

Over Q, the projectors `(I+R)/2` and `(I−R)/2` split the port space into two invariant sectors, each of dimension six. The antisymmetric sector has basis `e_(i,j)−e_(j,i)` for i<j and is naturally identified with `Λ²Q⁴` by

    ei∧ej ↦ e_(i,j)−e_(j,i).

A vertex re-marking π genuinely acts on this sector as `Λ²Pπ`. This is transition compatibility for vertex permutations, not for the construction letters. The latter add a vertex and change the native port-space dimension, whereas the six-state machine retains a fixed decision alphabet. The rank-two ee obstruction in §1 therefore remains untouched.

The integral antisymmetric lattice is generated by the same difference vectors. However, the rational sector projectors divide by two: the sum of the integral symmetric and antisymmetric lattices has index 2⁶ in the full 12-port lattice. One must not discard this integral gluing when making determinant or Smith-form claims. Nor does commutation with run-swap say that the restricted M2 operator is the exterior square of some four-dimensional vertex operator.

This gives a concrete connection to Session103’s inherited bridge kernel: `(ea−es)∧er` means exactly the native integer vector `e_(a,r)−e_(r,a)−e_(s,r)+e_(r,s)`. Its wedge form lives in the observer’s antisymmetric sector, independently of the automaton state count. At general n that sector has dimension binomial(n,2), so its dimension is six only at n=4; the forbidden bridge example itself occurs at larger orders.

More positively, labelled vertex injections V→V′ induce genuine functorial maps `Λ²Q^V→Λ²Q^V′`; each construction insertion thus embeds the old wedge coordinates in a growing carrier. Merely embedding a vector does **not** show that a kernel vector for the old M2 remains in the kernel of the new M2. Session103’s last-forbidden-pair bridge/sink guard supplies precisely the additional kernel-transport argument for its particular certificate. This is a well-defined growing bivector system, not a constant six-dimensional automaton representation and not, without additional geometric data, a differential bivector field.

## 4. Four potentials give a real Plücker/additive relation

For x=(x0,x1,x2,x3), set d_ij=xj−xi. These six numbers are coefficients of the decomposable bivector `(e0+e1+e2+e3)∧x`. Consequently

    d01 d23 − d02 d13 + d03 d12 = 0.

This follows equally by direct expansion. With A=d01d23, B=d03d12, C=d02d13 it becomes A+B=C: the three opposite-edge matchings supply an honest additive relation among products. For strictly increasing x0<x1<x2<x3 all three are positive. Integrality follows from integer potentials; coprimality need not hold before gcd reduction. In fact this carrier imposes **no restriction at all on primitive positive abc triples**, as the following constructive theorem shows. No abc conjecture consequence is claimed. Most six-coordinate bivectors are not decomposable: e01+e23 violates this quadratic relation. Arbitrary six-state probability/count vectors therefore do not automatically describe a single wedge of two vectors.

### Universality after gcd reduction

**Theorem (elementary constructive reduction).** Every primitive positive integer triple A+B=C, up to exchanging A and B, is the gcd reduction of the three opposite-edge difference products of four strictly increasing integer potentials.

**Proof.** Suppose A<B and gcd(A,B)=1. Put d=B−A>0 and choose

    (x0,x1,x2,x3)=(0,d,2d,2B).

These are strictly increasing because 2B−2d=2A>0. Direct calculation gives

    d01*d23 = d*(2A) = 2dA,
    d03*d12 = (2B)*d = 2dB,
    d02*d13 = (2d)*(A+B) = 2dC.

Their common gcd is exactly 2d, since gcd(A,B,C)=1. Dividing by it recovers the prescribed triple. If A=B, primitivity forces (A,B,C)=(1,1,2); potentials (0,2,3,6) give products (6,6,12). This handles the only equal boundary case. ∎

Consequently, any radical/height inequality proved for **all gcd-reduced triples arising from this unrestricted four-potential carrier** would apply to all primitive positive abc triples. The carrier alone gives no smaller arithmetic family and no hidden radical saving. Additional restrictions on the potentials or additional arithmetic estimates would be substantive extra hypotheses. The construction is a session deduction independently checked here, not a historical novelty claim; the proposed formula was supplied by the coordinating agent.

## 5. An honest synthetic 4→6 lift of the three-state counter

One can always build a geometric realization once its non-canonicity is declared. For the inherited

    C = [[1,0,1],[1,1,1],[1,1,0]],

take V=Qt⊕H and L=1⊕C. Then

    Λ²L = C ⊕ Λ²C

on `(t∧H)⊕Λ²H`. This preserves the incidence exact sequence and acts as C on its three-dimensional quotient. It gives a legitimate six-dimensional extension of the counting dynamics, **not** the legality machine and **not** the four-state acceptance transition.

The second block is the contragredient twist `Λ²C ≅ (det C) C⁻ᵗ` after choosing a volume identification. Since det C=−1,

    char_C(z)=z³−2z²−z+1,
    char_(Λ²C)(z)=z³+z²−2z−1.

Thus even the natural synthetic extension contains a different dual three-dimensional block, not two copies of the counter. This construction shows exactly what extra choice is needed to turn the 3/4/6 analogy into a valid linear model, while the rank test proves it is not the actual output automaton.

### The dual block is integrally the shifted counter: an exact seven connection

There is more structure than distinct characteristic polynomials. In the wedge basis (h1∧h2,h1∧h3,h2∧h3), put

    W = Λ²C = [[1,0,−1],[1,−1,−1],[0,−1,−1]],
    S = [[0,−1,2],[−1,−2,4],[−1,−3,5]].

Direct multiplication gives `det S=1` and `W S = S(C−I)`. Therefore the synthetic six-dimensional extension is **integrally similar to C⊕(C−I)**. It is still not the legality machine. One way to find S is to take the cyclic basis (v,Dv,D²v) for D=C−I and the cyclic basis (v,Wv,W²v) for W, with v=(−1,−1,−1); both have determinant −1 and their change-of-basis matrix is S.

Write q(z)=z³+z²−2z−1. Then q(z)=char_C(z+1). For a primitive seventh root ζ and x=ζ+ζ⁻¹,

    q(x)=1+Σ_(k=1)^3(ζ^k+ζ^(−k))=0.

The three distinct real roots are `2cos(2πk/7)`, k=1,2,3. Thus the already inherited counting growth is `1+2cos(2π/7)`, while the exterior dual block has the unshifted three cyclotomic traces as eigenvalues. The known seven-related growth is not a new counting result here; the new assembly identifies its integral exterior-dual block.

For D=C−I, q(D)=0 yields

    (D²−2I)(D+I)=−I,
    D²−2I=−C⁻¹.

On the roots of q the Möbius transformation `φ(x)=−1/(x+1)` therefore equals `x²−2`. The latter sends ζ^k+ζ^(−k) to ζ^(2k)+ζ^(−2k), cycling the three real traces with order three. Directly, `φ²(x)=−(x+1)/x` and `φ³(x)=x` wherever defined; none of q's roots is 0 or −1. This is a genuine algebraic seven/counter/exterior connection, not an identification with a seven-state machine, a field of characteristic seven, or K4 vertex permutation symmetry.

## 6. Positive flows versus signed kernel

The complementary [positive-incidence note](../positive/POSITIVE_INCIDENCE_AND_ATOMS.md) uses this same boundary/potential pairing: a finite oriented graph has no nonzero nonnegative circulation exactly when it admits a strictly increasing vertex potential. This explains how a rank-three incidence map can still detect every nonzero positive flow in an acyclic K4 orientation despite its three-dimensional signed kernel. It links the linear projection to Session103’s D_C potential argument without identifying that graph with the six-state machine. The proof and its finite-versus-infinite descent boundary belong to that note.

## Verification

`python3 check_exterior.py` and `python3 -O check_exterior.py` both pass. The stdlib exact-rational checker verifies the transition table's reachability/minimality and quotient; all letter ranks and the rank-two ee obstruction; incidence exactness and orthogonality; Hodge star and its exchange of cuts/cycles; the synthetic block lift, characteristic polynomials, explicit unimodular shifted-duality intertwiner, and order-three root transformation; and the difference Plücker relation on 625 integer quadruples; and constructive recovery of every coprime pair 1≤A≤B≤100 (3,044 pairs), including the equal boundary case. Finite Plücker tests corroborate the algebraic proof, not replace it. Explicit exception guards survive optimization. Receipt: `results.json`.
