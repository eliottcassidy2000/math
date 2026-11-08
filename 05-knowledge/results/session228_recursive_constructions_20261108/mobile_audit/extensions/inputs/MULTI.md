# Exponentially large strong scalar fibers by separated triangles

Status: PROVED by an exact transfer identity and physical graph construction; independent audit requested. This addendum does not claim scalar interval coverage or Chan–Pak completeness.

## 1. Construction and theorem

Take r≥3, let J_r={1,3,5,...}∩[1,r−2], and let N=|J_r|=floor((r−1)/2). Fix 1≤t≤N. For each t-element subset P⊆J_r, start from the expanded binary word(01)^r and replace the pair at positions2j,2j+1 by a cyclic triangle for every j∈P. Positions are indexed from0. The root beats exactly one vertex of each replacement triangle because every replaced pair is01. Other bits specify its singleton adjacencies. Deletion blocks point left to right, as in PROOF.md.

Let H_r=(4^r+6r−1)/9. The resulting tournament W_(r,P) has2r+t+1 vertices and

    d=3^t, E=0, U=3^(t−1)(4^r+2),
    D=3^t H_r +3^(t−1) sum_(j∈P)(4^j−2).             (1)

Remove its initial singleton0 source and final singleton1 sink to obtain V_(r,P). This is strong, has2r+t−1 vertices, and

    d=E=3^t, U=3^(t−1)(2^(2r−1)+1),
    D=3^t H_r +3^(t−1) sum_(j∈P)(4^j−2).             (2)

For fixed r,t, all binomial(N,t) scalar values in either family are distinct. Indeed equality of the sums in(1) is equality of sums of distinct powers4^j, because t is fixed. Base4 digit uniqueness forces the same subset P. This is a statement about different scalar counts, not only nonisomorphic graphs.

Taking t=floor(N/2) for N≥2 gives at least2^N/(N+1) distinct scalar values in a single physical fiber, witnessed by strong tournaments of order2r+t−1. Thus the fiber cardinality grows exponentially with graph order along this explicit sequence. The fixed sidecar depends on r,t; no claim of a single fixed sidecar holding infinitely many integer D values is made.

## 2. Why separated defects add exactly

Use the six-state transfer matrices S0,S1,T1 and v,ell of PROOF.md. Put A=S1S0 and F=T1−3A. A normal01 pair contributes A. A replaced pair contributes3A+F. The key finite identities are

    A(A−4I)(A−I)^2=0,
    F A^h F=0 for h=1,2,3,4.

The first identity and an induction imply F A^h F=0 for every h≥1. These identities are obtained by exact integer matrix multiplication; the checker verifies them explicitly.

Because distinct j in J_r differ by at least2, any two replaced pair blocks have at least one unreplaced A between them. Expand the product of all r pair matrices. Every term with at least two F factors contains a factor F A^h F with h≥1 and therefore vanishes identically. Only the all-A term and single-F terms remain. For a fixed t-element set P this gives

    D(W)=3^t ell A^r v
          +3^(t−1) sum_(j∈P) ell A^(r−j−1) F A^j v.  (3)

The transitive-deletion base has

    ell A^r v=H_r=(4^r+6r−1)/9 for r≥1.

To verify this at every order, apply the cubic recurrence with polynomial(z−4)(z−1)^2 to the sequence ell A^r v for r≥1; the matrix identity supplies it after multiplication by A. The initial values at r=1,2,3 are1,3,9, respectively, and the displayed formula satisfies the same recurrence.

For j∈J_r, the single replacement at position p=2j has j≥1 and p≤2r−4. PROOF.md gives its scalar as B_r+4^j−2, where B_r=3H_r. Therefore

    ell A^(r−j−1) F A^j v=4^j−2.

Substitution in(3) proves the scalar formula(1).

## 3. Sidecars and strongness

The multiple-triangle sidecar identity in PROOF.md gives d=3^t and U=3^t(1+binary_value((01)^r))=3^(t−1)(4^r+2). The unreplaced final singleton1 forces E=0. Distinct pair replacements each replace two singleton vertices by three triangle vertices, explaining the vertex count.

The initial singleton0 and final singleton1 are untouched because every j lies in[1,r−2]. They are a global source and sink. Removing them preserves D,d; removing the source preserves U, while removing the sink halves U by the two append choices for path covers. The new final deletion block is the unreplaced singleton0 at position2r−2, so E=d. The new first deletion block is the unreplaced singleton1 at position1. The root beats that first block and loses to the last block. Every deletion block is strong and all between-block edges point forward. Hence the root reaches all vertices through the first block, and every vertex reaches the root through the last block. This proves strongness and(2).

## 4. Necessary separation and the remaining boundary

Separation is load-bearing: on word0101, replacing both adjacent01 pairs gives D=41. The independent single replacements have D=9 and13, while the unreplaced graph has D=3. A false additive formula would give3·9+3·13−9·3=39. The interaction is nonzero. Thus no unrestricted independent-digit law is claimed. The zero product F A^h F needs h≥1.

The theorem supplies exponentially many scalar choices with one sidecar, but the counts have large common arithmetic structure and widely spaced base4 contributions. It does not supply a prescribed odd target, a contiguous arithmetic progression of unbounded length, or a physical inverse descent covering all odds. Those remain separate research obligations.
