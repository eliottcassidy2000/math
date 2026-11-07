# Exact endpoint grading behind the binary cycle construction

**Status: PROVED elementary statements below; FINITE-EXACT implementation checks are separate. No literature-priority claim. The unrestricted Hamiltonian-path spectrum remains OPEN.**

## 1. Inheritance and exact conventions

For a finite tournament T, H(T) is its number of directed Hamiltonian vertex words. HC(T) counts directed Hamiltonian cycles modulo cyclic rotation, not rooted cycle words. Set HC(T)=0 for orders below three. For a marked vertex r, s_r(T) counts Hamiltonian paths starting at r.

The source snapshot is `c137c698e632cb76e232a79b48133987755e080a` on `eliottcassidy2000/math`; its only change after Session211 is a hypothesis-index reservation. This work does not modify the stale dirty checkout.

Closest inherited mechanisms:

- Session210, `session210/graded-spectrum-gates/PROOF.md`, Theorem 1: the two-vertex degree gate Gamma(T), with a->b->T->a, has HC(Gamma(T))=H(T). It proves the exact size-indexed equality between the gated cycle spectrum at n+2 and the path spectrum at n. It already supplies gate uniqueness, cycle decks, cover-polynomial transport, and finite invariant reconstruction.
- THM-4862 with slug `01-canon/theorems/THM-4862-one-exit-extension-law-the-readiness-digit-is-the-cycle-count-of-the-entry-attachment-fibonacci-bijection-ready-floors-through-order-10-and-the-full-digit-c2-sweep-realises-every-odd-number-in-307-26551.md`, part (b): adding a vertex whose unique exit is r converts s_r(T) into a Hamiltonian-cycle count. This endpoint-to-cycle interface is inherited, not new.
- Session211, `session211/score-permutohedron/PROOF.md`, Theorems 3-4: gate erasure exposes deletion and two-path responses, with explicit ungated HC=7 and HC=21 witnesses. Those examples already refute transporting the path omissions to all cycles.
- THM-1370 with slug `01-canon/theorems/THM-1370-h-spectrum-omits-7-21-all-n.md`: the omissions are about total Hamiltonian-path counts, not endpoint counts or unrestricted cycle counts.

Anchor: the unrestricted cycle-count support. Niche: a complete logarithmic-size endpoint code. Wildcard: compare one local forced exit with a degree gate. Concept board: path totals; marked endpoints; sole exits; binary subset counts; exact vertex grading; gate obstruction.

The canonical hostile is conflating H with s_r. The corrected near miss is attempting to realize forbidden total-path values 7 and 21 by a cycle construction. The least-used sidecar is a marked starting vertex together with a transitive chain and one apex.

## 2. The inherited endpoint-to-cycle correspondence

For a tournament T of order n>=2 and r in V(T), define F(T,r) by adding x with x->r and v->x for every v!=r.

**Lemma 1 (recovered).** HC(F(T,r))=s_r(T). Moreover, F(T,r) is strong exactly when s_r(T)>0.

**Proof.** Every Hamiltonian cycle must use x->r. Anchor it at x and remove x to get a Hamiltonian path of T starting at r. Conversely, such a path ends at a vertex other than r because n>=2, and that endpoint beats x, so the path closes uniquely. Positivity supplies a spanning directed cycle, hence strong connectivity. If F is strong, a standard elementary Hamiltonicity theorem for strong tournaments supplies a cycle, but that converse is not needed below: equivalently one may use the tournament theorem that a starting vertex has a spanning path exactly when it reaches every vertex. No claim about new provenance is made for either standard statement. QED.

For exact support equality the strong assertion can be omitted entirely: with E_n={s_r(T): |T|=n, r in T},

    {HC(G): |G|=n+1, G has a marked vertex x of outdegree one}=E_n,

because deleting x and recording its unique out-neighbor gives the inverse rooted construction. This is a marked object correspondence and shifts size by one, unlike the unmarked degree-gate correspondence shifting size by two. Zero is allowed. At n=1 the extension is a transitive pair with no cycle while s_r=1, so that case is explicitly excluded.

## 3. Binary endpoint coding on a transitive chain plus one apex

Fix m>=0 and an integer 1<=h<=2^m. Write

    h-1 = sum_(k=0)^(m-1) b_k 2^k,  b_k in {0,1}.

For m=0 this is the empty sum and h=1. Let T(h,m) have chain vertices 0,1,...,m+1, ordered by i->j iff i<j, and one apex z=m+2. Its apex arcs are

    z -> 0,       1 -> z,
    j -> z iff b_(j-2)=1, for 2<=j<=m+1.

All other apex arcs point from z to j. Mark r=1.

**Theorem 2 (complete binary endpoint code).**

    s_1(T(h,m)) = h,       |T(h,m)|=m+3.

**Proof.** A Hamiltonian path beginning at 1 visits z exactly once. The parts before and after z are increasing chain subsequences. Vertex 0 cannot occur before z, since the path starts at 1 and all chain arcs are increasing. Thus the part after z is nonempty and begins at 0. Its entry arc z->0 is always legal.

If the part before z is just vertex 1, there is exactly one path:

    1, z, 0, 2, 3, ..., m+1.

Otherwise, let j>=2 be its last chain vertex. Necessarily b_(j-2)=1, to make j->z legal. The path before z is 1 followed by any increasing subsequence of {2,...,j-1}, followed by j. There are exactly 2^(j-2) choices. Every chain vertex omitted there occurs after z in increasing order, beginning with 0; every such choice gives a valid Hamiltonian path. Conversely the two subsequences and j are recovered uniquely from the path. Therefore

    s_1(T(h,m)) = 1 + sum_(j=2)^(m+1) b_(j-2)2^(j-2) = h.

This is a literal partition of paths, not an asymptotic count or a finite-range observation. QED.

Connection contract: source is the binary expansion of h-1; target is a rooted tournament with a transitive vertex-deletion. The map writes each bit as one apex orientation. Preserved predicate is exact cardinality h of paths starting at the marked vertex. The scalar h alone does not encode arbitrary tournament structure or its total H. Sidecars are chain order, sentinel 0, root 1, and apex z. The smallest decisive tests are h=1 at m=0 and h=2 at m=1.

## 3.1 Positive endpoint support already on strong cores

For h>=2 take the minimal m=bit_length(h-1), so the highest binary digit is one. Every chain vertex reaches the last chain vertex, which beats z; and z reaches every chain vertex via z->0 followed by the increasing chain. Thus the coding core T(h,m) itself is strong. For h=1 the m=0 core is the directed triangle. Consequently every positive integer is a Hamiltonian-path starting count of a **strong** tournament on ceil(log2(h))+3 vertices. This count remains marked; its total H need not equal h.

In the strong sole-exit closure B, cutting the forced arc x->1 also gives a bijection with Hamiltonian paths having **both specified endpoints**: start1 and endx. Thus every positive integer is a fixed-start/fixed-end Hamiltonian-path count of a strong tournament as well. These are endpoint-resolved statements, not statements about total H.

## 4. Relation to the unified family

The main `PROOF.md` unifies two independently developed constructions. Its notation `B-(m,q,tau)` uses a as the sole-exit vertex and b as the coding apex; this note uses x and z respectively. Our original positive family is precisely its `tau=1, q=h-1` slice, up to interchanging the last two vertex labels. Reversing the apex arc produces `B+` and adds one cycle when m>=1. Its arbitrary-tau form supplies the literal marked binary-growth grammar.

The independent verifier `verify_endpoint.py` retains the original labeling and the order-four m=0 boundary. It checks 511 endpoint instances, both apex orientations (1,022 cycle instances), and 62 complete-permutation cycle instances, plus larger controls. The forward family is not always nongated: in this tau=1 slice it has a gate exactly when h=2^m, namely z->0. The original reversed family is nongated for every m>=1. These properties are independently tested, not conflated.

Reproduce: `python3 verify_endpoint.py > /tmp/endpoint.json` and compare with `ENDPOINT_RESULT.json`; repeat using `python3 -O`. The file uses explicit exceptions.
