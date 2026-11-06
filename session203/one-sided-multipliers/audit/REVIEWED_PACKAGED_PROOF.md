# Exact one-sided multiplier contexts and the finite ternary tower

Status: **PROOF CANDIDATE, sent for independent review.** This is a constructive context classification, not eventual odd-spectrum coverage. The root-amalgamation formula and zero-slack classification are inherited; no priority claim is made for them.

## 1. The actual context

Let (U,v) be any nonempty rooted tournament. Write h=H(U), s=s_v(U), d=H(U-v), and L for the number of ordered spanning two-path covers whose first path is nonempty and ends at v and whose second path may be empty. Put H(empty)=1.

For any rooted core (G,r), identify r with v and orient every cross edge from G-r to U-v. Denote this one-sided root amalgamation by G*U. For the core use h_G,e_G,R_G,d_G, where R_G counts ordered spanning two-path covers with a nonempty second path starting at r and an optionally empty first path. The inherited exact formula is

    H(G*U) = d_G h + (R_G-d_G)L + (e_G-d_G)s
             + (h_G-e_G-R_G+d_G)d.                 (1)

This is a literal graph construction, not a signed count being treated as a tournament.

Let P(U) use the two-vertex core r->a. Its profile is

    (h_G,e_G,R_G,d_G)=(1,0,2,1),

so H(P(U))=h+L-s.

For every integer m>=1, choose the core G_m=E(2^m+1) with its inherited marked root r. Write a=2^m+1. Its profile is

    (h_G,e_G,R_G,d_G)=(2a-1,a-1,2a,a).

Consequently

    H(G_m*U)=a(h+L)-s
             =a H(P(U))+(a-1)s.                   (2)

On the physically restricted input class s=1 this becomes

    n -> (2^m+1)n+2^m,    n=H(P(U)).               (3)

For m=1, an explicit four-vertex core is (2,8,3,5), rooted at0. It gives3n+2. For m=3 the core gives9n+8, numerically the same as two3n+2 steps. This does not assert that the first output is a new valid P(U') input or that the operation can be repeated indefinitely.

## 2. Complete classification in this one-sided context class

**Claim.** For a positive integer a, a one-sided rooted core G satisfying

    H(G*U)=a(h_U+L_U)-s_U for EVERY rooted tournament U   (4)

exists if and only if

    a=1, or a=2^m+1 for some m>=1.                 (5)

The positive direction is the edge core for a=1 and the explicitly constructed cores E(2^m+1) above.

For necessity, the four functions h,L,s,d are linearly independent on rooted tournaments. The singleton, edge with root last, three-vertex transitive tournament with root last, and directed triangle give, respectively,

    (1,1,1,1), (1,2,0,1), (1,4,0,1), (3,3,1,1).

These four vectors are linearly independent. Thus comparing (1) and (4) forces

    d_G=a, R_G=2a, e_G=a-1, h_G=2a-1.             (6)

If a=1 there is nothing more to show. Suppose a>1. Its merge slack is

    R_G+e_G-h_G-d_G=0.

By the all-order zero-slack theorem, the root r has outdegree at most one. It cannot be a sink, since a sink would give h_G=d_G, contrary to (6). Let its sole out-neighbor be c, and put A=G-r. Every Hamilton path of A except those ending at c can be followed by r, hence

    e_c(A)=H(A)-e_r(G)=a-(a-1)=1.

The exact near-sink extension formula gives

    H(G)=H(A)+R_c(A)-e_c(A),

so R_c(A)=a. Applying nonnegative merge slack again, this time to (A,c), gives

    0 <= R_c(A)+e_c(A)-H(A)-H(A-c)=1-H(A-c).

Every tournament, including the empty tournament, has at least one Hamilton path. Therefore H(A-c)=1, and the slack at c is zero. The inherited complete unit-deletion, zero-slack spectrum is {1} union {2^m+1:m>=1}. Since H(A)=a>1, conclusion(5) follows.

This uses the universal theorem in `session195/zero-slack/PROOF.md`, Sections2--4, not a bounded gadget scan.

## 3. The exact ternary stopping point

An indefinite one-sided family with readout3^k(h+L)-s would require a=3^k in (5). For k>=1, the equation3^k-1=2^m has only k=1,2:

- If k is odd, reduction modulo8 forces m=1, then k=1.
- If k=2j is even, (3^j-1)(3^j+1)=2^m. Both factors must be powers of two differing by2, so they are2 and4. Thus j=1 and k=2.

The observed cores at a=3 and9 cannot be extended to a=27 by ANY one-sided core of ANY order. This does not exclude mixed-mask contexts, two sensitive interfaces, a different encoding, or a different output profile. Nor does it exclude the scalar count53, which is known to occur; it excludes the required profile (h,e,R,d)=(53,26,54,27).

The unit-start restriction in (3) is real. Even when U has s=1, the constructed output need not be representable as P(U') with s(U')=1. No threshold or closed ternary induction follows from (2).

## Dependencies and local evidence

- Root-amalgamation: `session197/path-polynomials/ROOTED_ORDERED_AMALGAMATION.md`.
- Nonnegative/zero merge slack and exact unit-deletion spectrum: `session195/zero-slack/PROOF.md`, Sections2--4.
- The endpoint-gate worker independently found the a=3 and a=9 tensors by exact core enumeration through six vertices. This finite agreement is a control, not the proof of the classification.

This packet is a constructive complement to the fixed-root even-multiplier obstruction. That obstruction is not a dependency here. The verifier checks six explicit core factors and 138 graph pairs on 23 rooted inputs; the all-order classification follows from the proof and its stated inherited dependencies, not from finite testing.

    python3 verify.py > NORMAL.json
    python3 -O verify.py > OPTIMIZED.json
    cmp NORMAL.json OPTIMIZED.json
