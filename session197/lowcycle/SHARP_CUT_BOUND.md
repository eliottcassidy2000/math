# The sharp nine-cycle bound at a separating vertex

Session197. **All-order theorem, independent of the start-three enumeration.** Dependencies are the audited start-two indegree lemma and C1 pair lemma in `BRIDGE_TWO.md`, together with the separating-vertex cycle product in `C2_CUT_REDUCTION.md`. No literature-priority claim is made.

**Theorem.** Let T be a strong tournament with no strong bridge. If deleting some vertex r destroys strong connectivity, then T has at least nine Hamiltonian cycles and at least seven vertices. Both bounds are sharp. In particular, every bridge-free strong tournament with at most eight Hamilton cycles is 2-strong.

Here a strong bridge is an arc whose deletion destroys strong connectivity. Hamilton cycles are directed cyclic words modulo rotation. A separating vertex is one whose deletion leaves a non-strong tournament.

## Proof of the lower bound

Write the SCCs of T-r in their unique dominance order as D_1,...,D_k, with k>=2. Let

    A=sum_(a in D_1, r->a) s_a(D_1),
    B=sum_(b in D_k, b->r) e_b(D_k),

where s and e count Hamilton paths starting and ending at the indicated vertex. Cutting a Hamilton cycle at r gives the exact product

    C(T)=A * product_(1<i<k) H(D_i) * B.             (1)

All factors are positive. We prove A>=3; reversal gives B>=3.

If A=1, there is exactly one entry a into D_1. The arc r->a is the sole arc from the rest of T into D_1, so deleting it destroys strong connectivity. This contradicts the no-bridge hypothesis.

If A=2, either there is one entry a with s_a(D_1)=2, or two entries a,b with both start counts one. In the first case the audited start-two theorem supplies an indegree-one vertex u of D_1 different from a. In the second case D_1 is C1, and the audited C1 pair lemma supplies an indegree-one vertex u outside {a,b}. In either case r does not beat u, and every later component loses to u. Thus u still has indegree one in T. Its sole incoming arc is a strong bridge, again a contradiction.

Therefore A>=3 and B>=3. Equation(1) gives C(T)>=9. Neither D_1 nor D_k can be a singleton, since then its endpoint sum is at most one. A nontrivial strong tournament has at least three vertices, so |T|>=1+3+3=7. QED.

## Explicit sharp witness

Take two directed triangles U,V and a singleton r, with all arcs

    r -> U -> V -> r.

In the order r,U_0,U_1,U_2,V_0,V_1,V_2, one adjacency presentation is

    [14,116,120,114,33,65,17].

Deleting r leaves the two SCCs U->V. Every Hamilton cycle is r followed by a Hamilton path of U, then one of V, then r. Each directed triangle has three Hamilton paths, so C=3*3=9.

There is no strong bridge. Every arc has an alternative directed path with the same endpoints avoiding it:

- inside U: route through V and r;
- inside V: route through r and U;
- r->x in U: route through x's predecessor in U;
- x in V ->r: route through x's successor in V;
- x in U ->y in V: choose z in V other than y and w in U other than x, and use x,z,r,w,y.

Since the original graph is strong, replacing a deleted arc by its alternative path preserves every reachability relation. Thus all21 arc deletions remain strong.

This is also the only seven-vertex bridge-free strong tournament with a separating vertex, up to isomorphism. The order bound forces precisely two triangle SCCs and no other vertices; A,B>=3 then force r to beat all of the first triangle and lose to all of the last. No uniqueness claim is made for larger orders or for all C9 tournaments.

## Checks and scope

`verify_sharp_cut_bound.py` literally counts the witness's nine cycles, checks every arc deletion and every vertex deletion, and confirms its two SCCs after deleting r. It separately exhausts all33,867 labelled tournaments through order six to check the absence of a smaller bridge-free separating-vertex example and the cut-product identities on actual separating vertices.

Normal and optimized runs agree. These finite tests corroborate the proof; the nine-cycle lower bound at arbitrary order follows from the endpoint lemmas and equation(1). The earlier four-cycle bound remains valid but is sharpened here.
