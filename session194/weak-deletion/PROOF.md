# Weak-deletion closure and the remaining root-deletion sector

Session 194, 5 November 2026. **PROVED, independently reviewed.** The all-order weak-deletion sector is closed below. The unrestricted rooted minimum and the Hamiltonian-path spectrum with proposed exceptions `{7,21}` remain **OPEN**.

## 1. Inheritance and the target

For a tournament T, H(T) counts directed Hamiltonian vertex words; H(empty)=1. Write

    Delta(T,v)=H(T)-H(T-v).

The inherited Grinberg--Stanley odd-cycle collection formula, THM-002 (`01-canon/theorems/THM-002-ocf.md`), gives

    Delta(T,v)=2 sum_(odd directed cycles C containing v) H(T-V(C)).   (1)

A cycle has length at least three, is taken modulo rotation, and retains its directed cyclic order. To see (1), subtract precisely the odd-cycle collections avoiding v, select the unique cycle through v, and sum the independent remaining collections. The weight of that selected cycle is two. This also proves Delta is even.

The closest proved mechanism is Session 193's parity-resolved rooted path identity and extreme-arc closure; its independent review and the exact single-root-triangle minimum g(n) (rooted merge-envelope Theorem 6) are routed through the [published Session 193 guide](https://github.com/eliottcassidy2000/math/blob/169792fb5eeb483e351112b5a09706e6d7bd89f0/00-navigation/S193.md). The latter gives

    g(4..8)=2,4,6,10,16;
    g(9+3j)=26*5^j, g(10+3j)=46*5^j, g(11+3j)=78*5^j.

The canonical hostile is the H33 greedy-arc example from Session 191/193: reversing a reverse cross arc can raise Delta. The new hostile in section 6 shows that the rooted E merge theorem cannot be transferred to Delta. The underused sidecar is the signed *unrooted* two-path-cover convolution within each SCC of T-v.

Concept board: marked odd cycles; ordered SCCs of the actual deletion; signed two-bin path covers; the E-versus-Delta distinction; multiplicative atoms of the desired spectrum. The anchor is a genuine all-order sector closure, the niche is its stronger component-count multiplier, and the wildcard is the direct deletion-merge hostile.

## 2. Two path-parity identities

For a nonempty tournament D and a root a, define

    R_a(D)=sum_(nonempty paths P starting at a) H(D-V(P)),
    L_a(D)=sum_(nonempty paths P ending at a) H(D-V(P)).

A path is a directed vertex word without repetition. Its omitted vertices support one Hamiltonian path, possibly empty. These are equivalently ordered spanning two-path covers with the marked endpoint in the specified bin. Let e_a(D) and s_a(D) count Hamiltonian paths ending and starting at a.

The Session 193 elementary identity, valid for weak tournaments as well, is

    sum_(P starting a) (-1)^(|P|-1) H(D-V(P))=e_a(D),
    sum_(P ending a)   (-1)^(|P|-1) H(D-V(P))=s_a(D).       (2)

For completeness, the first signed sum F_a obeys

    F_a(D)=H(D-a)-sum_(a->b)F_b(D-a).

Induction and the terminal-vertex recurrence give F_a=e_a, starting from the singleton. Reversing all arcs gives the second identity. The unsigned sums are R and L. Also R_a(D)>=H(D) and L_a(D)>=H(D): split each Hamiltonian word immediately before, or immediately after, the named root. These maps are injective into the respective two-bin covers.

Define the unrooted two-bin convolution

    Q(D)=sum_(S subset V(D)) H(D[S]) H(D-S).

The empty/full terms alone show Q(D)>=2H(D) for nonempty D. Its signed analogue satisfies

    sum_S (-1)^|S| H(D[S]) H(D-S)=0                  (3)

for every nonempty D. Here is an OCF proof that fixes all parity conventions. Expand both H factors into odd-cycle collections. Their union is a collection Gamma of disjoint odd cycles in D. Assign every cycle and every uncovered vertex to one of the two sides. Each cycle contributes the signed assignment factor `1+(-1)^odd=0`; each uncovered vertex contributes `1+(-1)=0`. Because D is nonempty, at least one such factor occurs. The total for every fixed Gamma is zero. Thus (3) holds. Equivalently, the total weights of even and odd vertex-count paths P in D, allowing the empty path once and weighting P by H(D-V(P)), are both Q(D)/2.

Every root of a strong tournament lies at the end and at the start of a Hamiltonian path. One elementary proof for the end version takes a longest directed path ending at the root. Any outside vertex beating a path vertex inserts immediately before the first path vertex it beats. If no outside vertex can be inserted, the entire path set dominates its nonempty complement, contradicting strongness. Reversal gives the start version. Hence e_a,s_a>=1 on every strong SCC, including a singleton.

## 3. Ordered-deletion multiplier

**Theorem 1.** Let T be strong, with at least three vertices, and suppose S=T-v is not strong. Write its ordered SCCs as

    C1 -> C2 -> ... -> Ck,       k>=2.

Choose any a in C1 with v->a and b in Ck with b->v; strongness guarantees both. Put

    P=product_(1<i<k) Q(Ci),     with an empty product equal to one.

Then

    Delta(T,v) >= R_a(C1)L_b(Ck)+e_a(C1)s_b(Ck),             k=2;
    Delta(T,v) >= R_a(C1)L_b(Ck) P,                         k>=3.   (4)

In particular,

    Delta(T,v) >= H(S)+1,                                  k=2;
    Delta(T,v) >= 2^(k-2) H(S),                            k>=3.   (5)

Thus always Delta(T,v)>H(S), equivalently H(T)>2H(T-v), and Delta(T,v)>H(T)/2.

**Proof.** Retain the following cycles through v in (1). Choose a nonempty path in C1 starting at a; a nonempty path in Ck ending at b; and in each intermediate Ci choose a possibly empty path. Concatenate them in SCC order between v and its return to v. All intercomponent arcs are forward. The first and last arcs are valid by choice of a,b. The two endpoint paths guarantee at least two nonroot vertices. Retain the concatenation exactly when the total nonroot vertex count is even, making a directed odd cycle of length at least three.

For the other cycles, choose arbitrary odd-cycle collections independently within each unused part of Ci. There are no crossing cycles between SCCs, so their total weight is exactly the product of the omitted-part H values. The root cycle determines all its paths by intersection with the fixed components, including all empty intermediate paths. Its first and last nonroot vertices recover a and b. Hence this construction is injective and counts a genuine subfamily of the marked collections in (1).

If U is the total unsigned path-product weight and Z is its signed weight with sign (-1) to the total nonroot vertex count, the retained weight is U+Z: half the sum selects even count, and the root-cycle weight is two. Formula (2) makes the endpoint signed factors `-e_a(C1)` and `-s_b(Ck)`. For k=2 their product is positive `e_a s_b`. For k>=3 every intermediate signed factor is zero by (3). This proves (4).

Now use R_a>=H(C1), L_b>=H(Ck), e_a s_b>=1, and Q(Ci)>=2H(Ci). Ordered SCCs multiply H, so the products of component H values equal H(S). This gives (5). For k>=3 its multiplier is at least two and H(S)>=1; for k=2 the explicit +1 makes the inequality strict. Finally H(T)=Delta+H(S). No assertion about intermediate arc reversals is used. QED.

The two-SCC singleton boundary is included: if both SCCs are singleton then T is C3, H(S)=1 and Delta=2, and the first bound is equality.

### Independent cut-capacity derivation

The independent review supplies a second proof of (4) from the inherited physical-cut and ordered-join identities, routed through the published Session 193 guide above. For an attachment set X=N^+(v) in S,

    Delta(T,v)=sum_(i in X,j outside X,i->j) c_ij(S),

where every directed-arc capacity is positive. Strongness supplies a in C1 intersect X and b in Ck outside X, and a->b, so Delta>=c_ab. The ordered-join boundary formula gives exactly

    c_ab=R_a(C1)L_b(Ck) product_(1<i<k)Q(Ci)
             +[k=2]e_a(C1)s_b(Ck).

This is the right side of (4). The independent audit also computes c_ab directly as `sum_U[e_a(U)s_b(S-U)+e_b(U)s_a(S-U)]`, with a missing endpoint interpreted as zero, and verifies 500 capacities on 80 explicit ordered joins. Thus the same retained quantity has both a marked-cycle derivation and an exposed-path-cut derivation. The complete audit is in the companion structural-cut-audit packet; the original cycle proof above remains self-contained relative to OCF and does not need the capacity identity.

## 4. All-order exclusion from the sharp-minimum problem

The cited Busch strong-tournament minimum, inherited and source-audited in Session 193, is

    B(3q)=3*5^(q-1), B(3q+1)=5^q, B(3q+2)=9*5^(q-1)

on its positive strong-order domains, with B(1)=1 separately and no strong order two. Thus every strong T of order n>=4 with weak deletion satisfies

    Delta(T,v)>B(n)/2>g(n).                              (6)

The second comparison is immediate at n=4..8. At n=9,10,11 it reads `26<75/2`, `46<125/2`, `78<225/2`; increasing n by three scales both sides by five. This is an all-order residue argument, not extrapolation.

Because Delta is even, a convenient integer version of the first bound is

    Delta(T,v)>=2 ceil((B(n)+1)/4).

Combining (6) with Session 193's extreme-arc theorem leaves only the following sector for either a counterexample to the proposed global minimum or an additional equality case:

* T and T-v are both strong;
* A=N^-(v) and B=N^+(v) are both weak;
* no reverse cross arc goes from the sink SCC of B to the source SCC of A.

The single-triangle minimum is attained at every order, so an unrestricted minimizing rooted tournament, if below or at g(n), must have strong deletion. This strengthens the earlier both-weak/no-extreme reduction. It does not claim that this remaining sector is empty or that its increments have the desired all-order lower bound.

## 5. Relation to the {7,21} spectrum

The source object here is a strong tournament with a specified weak deletion. The map is a positive subfamily of root-covering odd-cycle collections, resolved by the ordered SCCs of that deletion. It preserves an exact contribution to Delta and discards other root cycles. The sidecar is the parity of each path segment; discarding it would give a false product count.

This is a lower-bound theorem, not an existence operation for a prescribed count. Session 193's prime-atom reduction still requires actual tournament witnesses for every odd prime other than 7 and every 7p with p>=5 prime, plus its finite exceptional atoms. Nothing in (4)--(6) constructs those witnesses or supplies the `3h+2` and `3h+4` induction steps. The unrestricted spectrum therefore remains OPEN.

## 6. A precise failed extension of the rooted merge theorem

The tempting inequality

    Delta(T,v) >= Delta(T[A],v) H(T[B])

for strong induced A,B covering T with singleton intersection w!=v is false. A minimum-order witness has outgoing rows

    (2,4,9,19,7),   v=0,
    A={0,1,2}, B={2,3,4}, w=2.

Both factors are directed triangles. Here H(T)=9, H(T-v)=5, and Delta(T,v)=4, whereas Delta(T[A],v)H(T[B])=2*3=6. The only root-containing odd cycles are `(0,1,2)` and `(0,1,2,3,4)`. Each contributes two. This is the smallest possible union order for two strong factors of order at least three sharing one vertex.

The first failed implication is replacing the tagged rooted response E=R+e by the positive deletion statistic Delta in Session 193's path-merge injection. Although both are positive rooted counts, their object types and gluing behavior differ. The strongest inherited survivor is the E merge itself; the new survivor for Delta is the explicitly SCC-resolved cycle construction in Theorem 1. No known theorem is retracted by this hostile.

## 7. Reproduction and finite-control scope

See `README.md` and `verify.py`. The verifier independently computes Hamiltonian subset counts in both endpoint directions, checks literal words and literal root-cycle subfamilies at small orders, exhausts the inherited complete bank through eight, checks the two-SCC singleton boundary and the gluing hostile, and checks the all-order residue comparison by its three prototypes plus a larger finite arithmetic control. Ordinary and optimized runs use explicit exception gates and must agree byte for byte. The proof, not the census, supplies the all-order conclusion.
