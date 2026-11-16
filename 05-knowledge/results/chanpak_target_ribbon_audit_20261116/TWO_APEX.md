# Two-apex bridge surgery and its exact conflict polynomial

Status: proved below; finite controls are separate. This packet makes no priority claim and no claim of closing the full Hamiltonian-path spectrum. No remote writes are performed.

## 1. The objects and the exact new construction

For a tournament T let H(T) count directed Hamiltonian paths. Its odd-cycle conflict graph Ω(T) has one vertex for each directed simple odd cycle (cyclic rotations identified), and adjacency means intersecting vertex sets. We use the standard tournament identity H(T)=I(Ω(T),2), with I the independence polynomial.

A binary word w=b_0...b_(m-1) represents a tournament U(w): its old vertices are the chain 0→1→...→m−1, with apex x pointing to i exactly when b_i=1. A **strong word** starts with 1 and ends with 0, so m≥2. Deleting x is transitive, and U(w) is strong. Define

    h(w)=H(U(w)),
    A(w)=sum_i b_i 2^(m-1-i),
    B(w)=sum_i (1-b_i)2^i.

A and B are positive even integers for a strong word. H is obtained exactly from (h,A,p)=(1,0,0) by

    append 1: (h,A,p) -> (h,2A+1,1),
    append 0: (h,A,p) -> (h+A+p,2A,0).

This elementary one-apex counting recurrence is inherited from Session25, not a new theorem here. The standalone checker also computes H directly on the literal adjacency matrix.

Take two strong words u,v with distinct apices x,y. Form U(u)→U(v), then reverse its **single bridge arc** x→y to y→x. Call the result R(u,v). Put

    h=h(u), j=h(v), P=A(u)+1, Q=B(v)+1,
    a=(h-1)/2, b=(j-1)/2, c=(PQ-1)/2.

**Bridge theorem.** The result is strong, has feedback vertex number exactly 2, and

    Ω(R(u,v)) = (K_a ⊔ K_b) ∨ K_c,                       (1)
    I(Ω,z) = 1+(a+b+c)z+ab z²,
    H(R(u,v)) = hj+PQ-1.                                 (2)

Here ⊔ denotes disjoint union and ∨ denotes graph join. This is a full isomorphism description of the conflict graph, not just a numerical coincidence at z=2.

### Proof

Every old odd cycle in U uses x, and every old odd cycle in V uses y. The old cycle families are cliques of sizes a and b, with no conflicts between them. The sizes follow either directly from one-apex cycle counting or from H=1+2a.

A new cycle crossing the two modules must use the only arc pointing backwards across the cut, y→x. Conversely it is obtained uniquely by an increasing rooted path in U beginning at x, an increasing rooted path in V ending at y, and their intervening forward cut arc. A rooted path is allowed to consist of its apex only. The pair of singleton paths is not admissible, since x→y has been reversed; but its total vertex count is even and it does not enter the odd-cycle count below.

There are P rooted paths in U beginning at x: the singleton, and paths whose first chain vertex is i with b_i=1, followed by any subset of the m−1−i later vertices. Their signed even-minus-odd count is

    -1 + sum_i b_i(1-1)^(m-1-i) = -1+b_(m-1) = -1.

The corresponding count for paths in V ending at y is Q, again with even-minus-odd count −1 since v starts in 1. Therefore the number of pairs with odd total vertex count is (PQ−1)/2=c. Each mixed odd cycle contains both x and y; it intersects every old odd cycle and every other mixed odd cycle. This proves (1) and (2).

For strongness, x reaches every vertex in U, every vertex in V reaches y, all other cross-cut arcs point forwards, and y→x closes the reachability circuit. For feedback number, deleting x and y leaves a transitive chain. But the old cyclic modules contain two vertex-disjoint triangles, so a one-vertex deletion cannot make the result acyclic. Thus the number is exactly 2. ∎

## 2. Neutral middle padding and a genuine coefficient-preserving doubling law

Insert k new chain vertices between the old U-chain and V-chain. Every old U vertex (including x) beats them, and they beat every old V vertex (including y); retain the reversed bridge y→x. Equivalently append k sink letters 1 to u before applying the bridge operation. Call the result R_k; k=0 is R above.

For every k≥1,

    Ω(R_k) = (K_a ⊔ K_b) ∨ K_(2^(k-1)PQ),
    I(Ω(R_k),z) = (1+az)(1+bz) + 2^(k-1)PQ z,
    H(R_k) = hj+2^k PQ.                                  (3)

The feedback number and strongness stay exactly as above. The appended sinks create no new cycles within U, multiply the total rooted-path count P by 2^k, and change its even-minus-odd count to 0. Thus the mixed odd-cycle count is 2^(k-1)PQ. This proves (3). Placing some of the new vertices as sources of V instead makes the same tournament after the evident relabelling: all occupy the same middle chain interval with the same ports.

In particular the quadratic coefficient ab stays fixed for all k, while the mixed universal clique doubles. Writing F=hj gives the exact state law

    H(R_1)=2H(R_0)-F+2,
    H(R_(k+1))=2H(R_k)-F  (k≥1).                         (4)

The exceptional +2 at the first step is essential. A single unqualified doubling recurrence would be false.

For fixed u,v this entire infinite family occupies O(log X) values below X. It does not fill an interval or all admissible residues. Indeed H(R_k)≡hj (mod 2^k) for k≥1.

## 3. Sharp endpoint bounds and limitations of this bridge construction

For every strong one-apex word,

    h ≤ A+1 ≤ 2h-3,   h ≤ B+1 ≤ 2h-3.                   (5)

The left equality for A holds exactly for 10^r, r≥1; the right equality exactly for 1^r0. The equalities for B are the reverse-complements of these.

To prove the left bound, track D=A+1-h. After appending 1, D'=D+A+1; after appending 0, D'=D-p. Starting with the first 1 gives D=1. The first ensuing zero can bring D to 0 iff there has been just one 1; every later 1 increases D by at least 3 and the next zero subtracts only 1. Thus D≥0 at every strong endpoint, with the stated equality case. For the right bound, write w=s0. If s ends in 1, h(w)=h(s)+A(s)+1≥A(s)+2, with equality iff s has h(s)=1, hence s consists entirely of 1s (it starts in 1 and a later zero would form a directed triangle). If s ends in 0, it is already strong and h(s)≥3, so h(w)≥A(s)+3. Since A(w)=2A(s), this proves A+1≤2h−3 and its equality case. Reverse-complement preserves h and exchanges A,B.

Consequently the unpadded bridge obeys the sharp sandwich

    2hj-1 ≤ H(R_0) ≤ hj+(2h-3)(2j-3)-1.                 (6)

Both endpoints are attainable with appropriate dyadic-fan modules. Since h,j≥3, every bridge output is at least 17. This construction therefore cannot supply 13 (which is physically realizable by other two-apex tournaments), and it supplies neither 7 nor 21. For 21 specifically, (6) forces h=j=3, whose unique strong modules and P=Q=3 give 17; padding gives 9+9·2^k≥27.

For u=10, so U=C3, the unpadded formula specializes to

    H(R_0)=3(h(v)+B(v))+2 ≡2 (mod 3).                    (7)

Thus this useful repair subfamily misses both other residue classes. Its padded outputs are 3h(v)+3·2^k(B(v)+1), divisible by 3 for every k≥1. None is a prime. These are limitations of the stated construction, not impossibility assertions about arbitrary two-apex tournaments.

A completely explicit four-parameter family uses u=1^r0^s and v=1^t0^w, all exponents positive. Substitute

    h=(2^r-1)(2^s-1)+2, P=2^s(2^r-1)+1,
    j=(2^t-1)(2^w-1)+2, Q=2^t(2^w-1)+1

into (2) and (3). In particular u=10^s, v=1^t0 attains the lower bound

    H=2(2^s+1)(2^t+1)-1.

## 4. Exact repairs, and what “missing” means

Literal examples, all checked on adjacency matrices:

| u | v | H(R_0) | I(Ω,z) |
|---|---|---:|---|
| 10 | 110 | 29 | 1+10z+2z² |
| 10 | 100 | 35 | 1+13z+2z² |
| 100 | 110 | 49 | 1+16z+4z² |
| 10 | 1010 | 59 | 1+21z+4z² |
| 10 | 1000 | 71 | 1+27z+4z² |

These five counts are absent from the **entire one-apex spectrum**, not just a fixed-order search: the one-apex lower bound h≥2^(m−2)+1 for strong words makes the word enumeration finite for each target. Every nontransitive one-apex tournament reduces to a strong one-apex core by stripping source/sink padding, which preserves h. The checker enumerates every eligible word, with this proved length cutoff. The examples are exact repairs of a structural restriction; they were already known to lie in the unrestricted tournament spectrum. No priority claim is made for their existence.

In contrast, there are **no global missing congruence classes** for one-apex tournaments beyond parity. Session25 already proves every admissible residue modulo every modulus, even in the reset-word family (1^L0^L)^t. Its count of values below X is only O(log X log log X). Thus modular completeness is not integer coverage. Pinned prior source:

https://github.com/eliottcassidy2000/math/blob/1b91742ce55a202b73c80d56561ba56d00b503a6/05-knowledge/results/session25_spectrum_wiring_20261021/session25-spectrum/PROOF.md

The one-apex length bound, classification and positive-density missing-value result are pinned at the same commit:

https://github.com/eliottcassidy2000/math/blob/1b91742ce55a202b73c80d56561ba56d00b503a6/05-knowledge/results/session26_random_bridges_20261022/session26-tournament/APEX_CLASSIFICATION_AND_SPECTRAL_LAW.md

## 5. General two-apex guardrail

For any tournament with a distinguished feedback set {x,y}, partition odd cycles into x-only, y-only, and both-root families of sizes a,b,c. Then

    I(Ω,z)=1+(a+b+c)z+d z²,   0≤d≤ab.                   (8)

Here d is the number of vertex-disjoint pairs consisting of an x-only and a y-only cycle. All other pairs intersect. In particular the discriminant is nonnegative:

    (a+b+c)²-4d = (a-b)²+2c(a+b)+c²+4(ab-d) ≥0.

Thus every such quadratic has real negative roots when d>0 (and a single negative root when d=0 and a+b+c>0). An abstract pair (a+b+c,d) satisfying these inequalities is NOT automatically tournament-realizable. The bridge construction realizes the extremal d=ab together with the specific, physically forced c above; treating c as an independently tunable nonnegative integer would discard the main realization constraint.

## 6. Complete finite compiler and verification scope

For an output ≤X in the bridge family, m=|u|, n=|v| satisfy

    2^(m+n-2) < (A(u)+1)(B(v)+1) < H(R_0)+1 ≤X+1.

The checker safely enumerates m+n≤floor(log2 X)+2, m,n≥2, all strong words at those lengths, and every admissible k using (3). This covers every bridge output ≤X at every order, including neutral-padding outputs. The finite receipt records the exact set at X=10000 and which of those values are absent from the one-apex spectrum, whose own proved cutoff is m≤2+floor(log2(X−1)). At X=10000 the exact result is 4004 bridge values, coming from 20181 word/padding presentations. Of these, 1941 are absent from the entire one-apex spectrum. The one-apex set has 2595 values, so the union has 4536 of the 5000 positive odd integers in this range. These are **finite exact** certificates, not a theorem of asymptotic or all-target coverage. The repaired values at most 100 are exactly 29,35,49,57,59,71,75,89,93.

Independent controls enumerate literal directed odd cycles and their intersections; count Hamiltonian paths by subset DP; compute strongness and a two-vertex feedback deletion; test all eligible short module pairs and padding; and check the endpoint inequalities and their equality characterizations. All gates use explicit exceptions, not Python assertions. Run `python3 check_two_apex.py` and `python3 -O check_two_apex.py`; receipts must be byte-identical.
