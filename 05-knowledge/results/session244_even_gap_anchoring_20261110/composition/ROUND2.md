# Round2: corrected inverse, efficient endpoint input, and primitive-response lifting

Status: exact physical identities and conditional target descent, not complete target selection.

## 1. Correction and independent algebra audit

The first draft reversed both signs in the F− inverse. The correct formula is

    h=(h′−3e′)/2, e=3h′−10e′.

It returns (9,3) from (171,51). `verify.py` now explicitly checks this inverse in all31375 integer cases; normal and optimized replay are unchanged because the forward identities were correct. The neighbouring `logic/ROUND2.md` inverse differences are also correct: fixed output endpoint and gap2N require (Δh,Δe)=(N,6N) for F− and(N/7,6N/7) for F+. Its determinant/content conclusions do not require correction.

The structural classification previously described as an advance is **inherited**, from Session231 `descent/ROUND2.md`, Sections1–2, with the degree obstruction already in THM-4873/4875. `STRUCTURE.md` now explicitly credits this inheritance and labels the independent reachability proof as an alternative proof. The small-count and prime gcd statements are scoped corollaries.

## 2. Linear-size synchronized gap compiler

Use the efficient endpoint compiler proved by this session's selector lane. For N≥1 it supplies a strong source

    (D,E,U,d)=(g+2N,N,2g+2N−ε,g+N−ε),

where1≤g≤N and ε∈{0,1}. The N=1 boundary uses the singleton response source; for N≥2 the inherited inverse O/Z response word supplies g=g(N).

Substituting into the new five-letter compiler gives

    L_N=31g+37N−11ε,
    R_N=31g+39N−11ε,
    E_out=d_out=9g+11N−3ε,
    U_out=H_out.

Thus both outputs are strong, have the same order, and

    R_N−L_N=2N,  R_N≤70N.

Their order is at most ceil(log2 N)+8 for N≥2 and8 for N=1. This is a quantitative improvement for this synchronized-interface construction; it does not improve the inherited compact one-vertex gap pair's5N bound. It still does not place the gap at an arbitrary lower endpoint.

## 3. A mixed-tail endpoint-one source, then minimal positive response slack

Begin with the directed edge0→1 and append its first Z vertex2, obtaining a triangle. Append any finite O/Z tail at the fixed root0, either empty or ending in Z. Reroot at vertex2.

Its endpoint count is exactly1. Vertex2 has only vertex1 as an incoming neighbor, and vertex1 has only vertex0 as an incoming neighbor. Every path ending2 therefore has the forced suffix0,1,2. The remaining tail vertices form a transitive chronological chain, so their order is forced. If the tail is nonempty, its last vertex points to0 exactly when the final letter is Z. This constructs the unique path. An empty tail is the triangle case. The source is strong because its final Z closes all old vertices through root0.

Write this source state as (D,1,U,d). Apply J at vertex2. The result has

    (H,e,R,d′)=(D+U−1,D−1,D+U,D).                 (1)

Consequently R=H+1, and gcd(H,R)=1. This is exact control of a response coordinate, stronger than merely finding a coprime pair of total counts. The newly rooted output is strong because it points to an old root that reaches all vertices, and old nonroot vertices point to it.

For the selector lane's special nonempty tail O^k Z, k≥1, its proved state is

    D=3·2^k+3, U=4·3^k+6·2^k−1,

so (1) gives

    H=4·3^k+9·2^k+1, R=H+1.

The parameter restriction k≥1 is essential. At k=0 the separate graph edge.ZZ has D=5,U=9 and yields H=13,R=14; the displayed exponential formula does not apply. For k=1,2 the new H values are31,73. The general mixed-tail source is not restricted to this one-parameter family.

## 4. Exact affine target descent through an arbitrary suffix tournament

Let A be any strong rooted tournament with response R_A=H(A)+1. Put h=H(A). Choose any nonempty tournament B, writing b=H(B), q=Q(B). Form A⇒B and append a new vertex that beats exactly the marked root of A, marking that new vertex. The inherited saturation factorization gives

    H_out=(b+q)h+q, E_out=bh, R_out=H_out.

This is an exact affine map with known intercept, not a fitted scalar model. For a prescribed odd target M, this construction reduces the problem to

    h=(M−q)/(b+q),

provided the quotient is a positive odd integer with an actual strong realization satisfying R=h+1. Since b+q≥3, the target rank strictly decreases. The necessary coprimality test for a prime target is gcd(b,q)=1, agreeing with `STRUCTURE.md`.

For B transitive of order j≥1, b=1,q=2^j, giving

    H_out=(2^j+1)h+2^j, E_out=h,
    gcd(H_out,E_out)=1.

Thus every mixed-tail source in Section3 produces an explicit primitive saturated family. Taking B to be the four-vertex count5 tournament gives b=5,q=32 and H_out=37h+32. That is another exact residue branch with gcd(b,q)=1, rather than a prime guarantee.

What remains is substantial: the descended h must carry the response-h-plus-one sidecar. Ordinary attainment of h is insufficient, and no complete invariant guaranteeing this condition has been proved. Nor do the transitive-B residue classes alone cover all targets: M=2^a−1 has M+1 a power of2, so no odd divisor2^j+1 can satisfy the required congruence M≡−1 mod(2^j+1). General B or a different physical source class is necessary for those targets. This is a precise positive descent with an unclosed physical membership condition, not a solution of Chan–Pak.
