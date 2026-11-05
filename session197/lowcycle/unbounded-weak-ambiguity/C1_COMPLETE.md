# Every uniquely Hamiltonian strong tournament has a strong bridge

Session 197. **All-order structural completion of the C=1 class.** This separate proof closes the qualification deliberately retained in `ROOTED.md` and `FIBONACCI_CLASS.md`; those frozen packets are unchanged. No literature-priority claim is made.

Let C(T) count directed Hamiltonian cycles modulo cyclic rotation. A strong bridge is an arc whose deletion destroys strong connectivity. A strong tournament is 2-strong if deleting any one vertex leaves a strong tournament.

## 1. A unique cycle forces a separating vertex

**Lemma 1 (deletion and insertion).** Suppose T is strong with unique Hamiltonian cycle C, and T-v is strong. Then deleting v from the cyclic vertex word C leaves a directed Hamiltonian cycle of T-v. Moreover that smaller cycle is unique.

**Proof.** By the standard Hamiltonian-cycle theorem for strong tournaments, take any Hamiltonian cycle D of T-v. Strongness of T gives both an in-neighbor and an out-neighbor of v among the vertices of D. On the cyclic word D there is consequently an adjacent pair x,y with x->v and v->y: follow the cycle from an in-neighbor until the first out-neighbor. Replace x->y by x->v->y. The resulting Hamiltonian cycle of T must be C. Deleting v recovers D, so D=C-v and is independent of the chosen cycle. QED.

**Lemma 2.** A 2-strong tournament has at least two Hamiltonian cycles.

**Proof by contradiction.** Suppose its sole cycle is indexed

    0->1->...->n-1->0.

Deleting any vertex and applying Lemma 1 shows that every second-step chord

    i -> i+2   (indices modulo n)                      (1)

is present. A 2-strong tournament has minimum indegree and outdegree at least two, so n>=5.

If n is odd, repeated steps of length two visit every vertex. They form a different Hamiltonian cycle, contradicting uniqueness.

Now let n be even, so n>=6, and examine the arc between vertices 1 and 4.

- If 4->1, the cyclic word

      0,2,4,1,3,5,6,...,n-1,0

  is a Hamiltonian cycle. Its changed arcs are second-step arcs from (1), the given 4->1, and ordinary cycle arcs. It differs from the original cycle.

- If 1->4, the step-two arcs form one cycle on the even vertices and one on the odd vertices. Delete the arcs 1->3 and 2->4, and add 1->4 and 2->3. These reconnect the two cycles into a single directed Hamiltonian cycle. Explicitly it is

      3,5,...,n-1,1,4,6,...,n-2,0,2,3.

  The two added arcs exist by the case hypothesis and the original cycle, respectively. It too differs from the original.

Both orientations are impossible. QED.

Thus every strong T with C(T)=1 has a vertex r for which T-r is not strong. The order-three and order-four cases are included: their minimum semidegrees cannot both be at least two, so they cannot be 2-strong either.

## 2. The separating vertex exposes a strong bridge

**Lemma 3 (cycle product across a separating vertex).** Let T be strong and let the ordered strong components of T-r be D_1,...,D_k, where k>=2 and all intercomponent arcs point forward. Then

    C(T) = A * (product_(i=2)^(k-1) H(D_i)) * B,

    A = sum_(a in D_1, r->a) s_a(D_1),
    B = sum_(b in D_k, b->r) e_b(D_k).                (2)

Here a singleton has one rooted Hamiltonian path, and the middle product is one when k=2.

**Proof.** Cutting any Hamiltonian cycle just after r yields a spanning path in T-r. The strong-component order forces it to traverse D_1,...,D_k consecutively. It starts at a vertex of D_1 beaten by r and ends at a vertex of D_k beating r. Conversely every such choice of component paths closes through r to a directed Hamiltonian cycle. This is a bijection, with r fixing the cyclic starting point. QED.

Strongness of T ensures that r has at least one out-neighbor in D_1 and at least one in-neighbor in D_k. Every vertex in a strong component starts and ends a Hamiltonian path, by cutting a Hamiltonian cycle; the singleton case holds separately. Thus all summands in A and B are positive integers.

**Theorem 4.** Every strong tournament with exactly one Hamiltonian cycle has a strong bridge.

**Proof.** By Lemma 2 choose a separating vertex r, and use (2). If C(T)=1 then A=1. Positivity forces exactly one vertex a in D_1 with r->a. Every later component loses to D_1, and every other vertex of D_1 beats r. Therefore r->a is the sole arc from V(T)\D_1 into D_1. Deleting it leaves no directed path from that nonempty complement into D_1. It is a strong bridge. QED.

This proof uses neither a tournament census nor a conjectured all-order gluing law.

## 3. Consequences for the already proved rooted grammar

The equivalences in `ROOTED.md` now apply to every strong C=1 tournament, not merely the bridge subclass. In particular each such tournament:

- has a vertex of outdegree one;
- has a root from which exactly one Hamiltonian path starts;
- is generated from a singleton by the exact one-seam grammar with G-b transitive and e_b(G)=1.

Conversely the nontrivial strong outputs of that grammar have C=1, as proved there. This is therefore a complete all-order classification of the C=1 class by that recursive grammar.

The exact enumeration in `FIBONACCI_CLASS.md` consequently becomes the count of **all** unlabelled strong uniquely Hamiltonian tournaments:

    a_3=1,       a_n=F_(2n-6) for n>=4.

Every such graph of order n>=4 is rigid by the proof in that packet, so the labelled count is `n! F_(2n-6)`. At order three the directed triangle has three automorphisms and there are two labelled tournaments, not six.

The previously certified arbitrary-depth, arbitrary-gadget bounded image in `ROOTED.md` now also gives the complete scalar C=1 spectrum through 5,001 at arbitrary graph order. Its missing odd values from 3 through 5,001 are exactly

    7,13,21,41,43,55,67,77,83,163,173,231,233,247,251,
    277,299,421,449,503,533,569,583,637,821,1067,1339,
    1861,2647.

The finite closed-state certificate and its growth inequalities, rather than the all-order classification alone, justify that target-bounded statement. No eventual conductor is claimed.

This settles the C=1 structural qualification. It does not classify the C=2 class, prove the full `{7,21}` spectrum, or change the separate H67 classification and its minimum C=3 conclusion.

## 4. Reproduction and scope

`verify_c1_complete.py` independently enumerates literal Hamiltonian cycles on every labelled tournament through order five, checks the deletion/insertion conclusion and separating-vertex product, and physically deletes the exposed bridge. It also checks both even-order cycle surgeries and the odd-order step-two construction on a long finite range, using only the forced arcs.

These checks are implementation and edge-case controls. Lemmas 1–3 and Theorem 4 are the all-order proof. The only standard imported existence fact is that a nontrivial strong tournament has a Hamiltonian cycle, already used with explicit conventions throughout `ROOTED.md`.
