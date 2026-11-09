# Exact path-orbit transfer and terminal stabilizer death

**Status:** PROVED by elementary finite group actions, with exhaustive checks through order five. The recurrence is a concrete state quotient, not a claimed new general theory of group actions. External priority is not claimed. Chan–Pak completeness remains OPEN. No algebraic-geometric theorem is imported into tournament theory.

## Inheritance and scope

Closest proved mechanism: `01-canon/theorems/LEM-003-aut-acts-freely-on-ham-paths.md`, which proves order-rigidity of spanning directed paths. The endpoint/subset Hamilton recurrence is classical. The refinement below retains the pointwise stabilizer of the used subset, which is absent from a quotient recording only state-orbit adjacency.

The inherited all-order absence of 7 and 21 in H/|Aut| is `01-canon/theorems/THM-4812-presentation-fibres-7-and-21-are-absent-at-every-order-orbit-product-lemma-vertex-transitive-lemma-closed-window.md`, including its appended independent repairs. We do not re-prove that theorem here. The SCC/OCF site and relative-color torsors of `01-canon/theorems/THM-4892-tournament-sites-plus-construction.md` are different carriers: this note uses directed path prefixes and their actual automorphism stabilizers.

Anchor: compute Hamilton-path orbits exactly. Niche: identify when residual symmetry vanishes. Wildcard: explain the symmetry content of Session231 moving-core switches. Concept board: prefix histories, endpoint-state quotients, restriction/extension, stabilizer mass, cut-state composition. Canonical hostile: a source followed by a cyclic triangle. Corrected near miss: using one initial coarse orbit and then raw extension multiplicities mixes two different counting measures. Least-used sidecar: pointwise, rather than setwise, stabilizers of used vertices.

## 1. An integral orbit recurrence

Let T be an n-vertex tournament and G any subgroup of Aut(T). For 1 <= r <= n, let P_r be the ordered directed paths on r distinct vertices. The endpoint state of p=(v_1,...,v_r) is (S,v_r), where S is its used vertex set. Write B_r for the G-orbits of all states (S,u) with |S|=r and u in S. States not supporting a directed path are allowed and have count zero.

For each b in B_r, choose a representative (S,u), and put

    K_b = {g in G : g fixes every vertex of S},     k_b = |K_b|.

The subgroup itself depends on the representative by conjugation; k_b does not. Let c_b be the number of G-orbits of paths whose endpoint state lies in b. For a in B_(r+1), put

    N(b,a) = #{v outside S : u -> v and (S union {v},v) belongs to a}.

These numbers do not depend on the representative of b. Define

    W(b,a) = N(b,a) k_a / k_b.                         (1)

**Theorem 1.** Every W(b,a) is a nonnegative integer, and

    c_a = sum_(b in B_r) c_b W(b,a).                  (2)

At level one, c_b=1 for every vertex orbit. At level n,

    sum_(a in B_n) c_a = H(T)/|G|.                   (3)

**Proof.** The stabilizer in G of an ordered prefix p with used set S is exactly K_b: fixing the sequence is equivalent to fixing every visited vertex. To extend a fixed p into endpoint type a, choose v in the set defining N(b,a), modulo K_b. This set is K_b-invariant. For each such v the stabilizer in K_b is the pointwise stabilizer of S union {v}, whose size is k_a. Thus every K_b-orbit in this set has size k_b/k_a. Its number of orbits is (1), proving integrality. It depends on (S,u), not the earlier order in p, so every path orbit of type b has exactly W(b,a) extension orbits of type a. Restriction deletes the last vertex and sends an extension orbit to a unique prefix orbit; hence the contributions in (2) are disjoint and exhaustive. Single-vertex paths have one orbit per vertex orbit. A full path fixes all vertices, so its stabilizer is trivial; (3) follows. QED.

The same proof holds for any finite digraph. Only the terminal bound in section 3 uses the tournament hypothesis.

## 2. The two correct measures and the incorrect hybrid

The action groupoid [P_r/G] assigns each path orbit the reciprocal of its stabilizer order. The contribution of endpoint type b is

    m_b = c_b/k_b.

Dividing (2) by k_a gives the simpler mass recurrence

    m_a = sum_b N(b,a) m_b.                           (4)

Thus N is not an intrinsically invalid matrix. It transports groupoid mass; W transports coarse path-orbit counts. They have different initial conditions: c_b=1 and m_b=1/k_b at level one. At the spanning layer k_b=1 and the two measures agree. Moreover sum_b m_b=|P_r|/|G| at every level, by orbit-stabilizer.

**Minimal tournament hostile.** Let T=1 => C3 with the source labeled 0 and the triangle labeled 1,2,3; adjacency bitmasks are [14,4,8,2]. Here G=C3, H=3. The singleton source state has k_b=3. Its three possible next vertices all belong to one target state orbit, with k_a=1. Therefore N(b,a)=3 but W(b,a)=1. The integral recurrence gives one full path orbit, as it should. Initializing singleton orbit counts at one and then using N gives three instead. For n<=3, no tournament permits such a nontrivial correction between nonempty layers: a prefix with a nontrivial pointwise stabilizer would leave at least three unused vertices. The hostile is therefore minimal in this precise sense.

**Positive control.** C3 => 1 => 1 has H=3 and |Aut|=3, but no transfer correction on a nonempty prefix; both integer recurrences happen to give one. The direction of the ordinal padding matters.

Loss ledger: source = actual directed prefixes; target = G-orbits of (used set,last vertex); map = forget prefix order then divide by G; preserved = valid next-vertex relation and, with W, exact orbit count; lost = prefix multiplicity and unused-vertex symmetry; restoration sidecar = c_b and k_b; decisive test = source => C3. No unproved path-gluing assumption is needed.

## 3. Uniform terminal stabilization, with a sharp boundary

**Theorem 2.** For a tournament and any G <= Aut(T), k_b=1 whenever |S|>=n-2. Thus W(b,a)=N(b,a) for every transition whose source uses at least n-2 vertices. Every prefix action P_r is free in this range.

**Proof.** An automorphism fixing S pointwise is supported on its complement. If the complement has at most two vertices, a nonidentity automorphism could only transpose two vertices. That reverses the unique arc between them, which an automorphism cannot do. Hence it is the identity. QED.

The bound is uniformly sharp for n>=4: take a transitive chain on n-3 vertices followed by C3, and take the prefix consisting of that entire chain in its unique order. The remaining triangle can rotate, so k_b=3 at r=n-3; its next extension has k_a=1. For n=3 the analogous nontrivial group is at the empty-prefix level. No inverse-limit theorem is claimed: this is an exact finite terminal-layer statement.

More generally, if every nonidentity element of G moves at least mu vertices, then k_b=1 for |S|>n-mu. This is the same proof with the actual minimal motion substituted for the universal tournament bound three.

## 4. Seven-coordinate cut states and moving-core gains

For a nonempty tournament C and cut M subset V(C), retain Session231's convention

    state(C,M)=(h,q,A,G,T,S,E).

Here h counts Hamilton paths; q counts ordered pairs of disjoint directed paths covering all vertices, allowing either path to be empty. In such a pair (P,Q), A requires the first vertex of Q to be in M, G requires the last vertex of P to be in M, and T requires both conditions. S and E count full paths whose first and last vertex, respectively, lies in M. Empty distinguished endpoints never meet a condition.

**Proposition 3.** If K <= Aut(C) preserves M setwise, then every one of these seven integers is divisible by |K|. Dividing the complete state by |K| is compatible with ordinal composition, provided the ordinal pair is quotiented by the explicitly specified product action K x L.

**Proof.** An automorphism fixing an ordered spanning pair (P,Q) fixes each path sequence pointwise and therefore every vertex. Thus the action on ordered spanning path pairs is free, just as it is on spanning paths. All cut conditions are K-invariant. Orbit-stabilizer proves divisibility for every coordinate. Session231's ordinal state formula is bilinear in the two complete states:

    (h,q,A,G,T,S,E) * (j,k,B,J,U,V,W)
      =(hj,qk,Ak+hB,Jq+jG,AJ+Tj+hU+EV,Sj,hW).

Its combinatorial decomposition is equivariant for K x L, and division by |K||L| is precisely composition of the divided states. This claims nothing about whether K x L is the full automorphism group. QED.

**Corollary 4.** Suppose two star-closed singleton-padded constructions use the same marked core (C,M), possibly at different positions. Every cut-preserving automorphism of C extends by the identity on padding and closure vertices. If their Hamilton counts differ by delta, then

    2 |Aut(C,M)| divides delta.                       (5)

Indeed both counts are divisible by |Aut(C,M)|, and delta is even by Rédei; the automorphism order is odd. In particular a unit step delta=+/-2 requires a rigid marked core.

This is an explanatory consequence, not a stronger new unit-gain obstruction: Session231's inherited content filter gcd(h,A,G)=1 already implies this rigidity by Proposition 3. The converse fails. On C3 with a singleton cut, Aut(C,M)=1 but (h,q,A,G,T,S,E)=(3,12,3,3,0,1,1), so gcd(h,A,G)=3. Rigid cut does not mean rigid underlying tournament and does not establish a unit loop.

The Session231 twelve-bit witness has adjacency [4,1,18,7,75,31,175,63], cut 124, state (45,1440,878,577,332,25,14). Direct checking finds its underlying core already rigid, and gcd(45,878,577)=1. Its previously proved moving gain -2 is consistent with (5); neither its interval size nor its baseline-selection scope is enlarged here.

## 5. Relation to 7 and 21 and to the source papers

Equation (3) computes the orbit count on the exact carrier used in the inherited quotient-spectrum question. For G=Aut(T), THM-4812 excludes 7 and 21. The recurrence alone does not prove that exclusion and does not prove that every other odd value occurs. It does explain how stabilizer information must pass through partial states before the final free-orbit quotient is reached.

The precise source pins and limited source summaries are in SOURCE_NOTES.md. Our use of them is methodological: preserve transport before compressing, and prove an explicit bound for disappearing extra data. There is no functor between these tournament states and either paper's geometric objects, no imported cohomology theorem, and no p-adic explanation of the forbidden values.

## Reproduction and finite verification

Run from this directory:

    python3 -B verify.py > CHECK.json
    python3 -B -O verify.py > CHECK.optimized.json
    cmp CHECK.json CHECK.optimized.json

The checker has no repository imports and uses explicit exceptions, not assertions. It enumerates all 1,099 labeled tournaments on orders 1 through 5; independently lists every directed prefix by literal permutations; canonizes actual path orbits under the enumerated automorphism group; compares every one of 74,521 endpoint-state orbit counts against (2); and checks all terminal stabilizers. There are 328 transfer edges with W != N. A proper-subgroup control uses the diagonal C3 inside Aut(C3=>C3) and gives H/|G|=9/3=3. The seven-coordinate divisibility check covers all 1,098 (tournament,cut) cases through order four. The short12 core state is independently recomputed by literal ordered path-cover enumeration. This is a finite verification, not an exhaustive large-order spectrum census or an independent audit of the source papers.
