# Induced one-subdivisions from a dyadic guard and a mod-seven split

2026-10-19. Independently [reviewed](../INDEPENDENT_REVIEW.md); retained as an alternative construction and provenance. Its height constant is superseded by the [main pair-extension proof](../SCALED_PAIR_EXTENSION_AND_TWO_DEGENERATE.md). No external-priority claim.

## The result

For every finite simple graph G on n>=1 vertices, its exact one-subdivision (replace each edge by a path of length two, using a private new vertex) is an induced subgraph of the positive-integer square-sum graph. All labels are below 50176 n^4. This is an explicit deterministic formula; no search, prime distribution, or factorization is required.

The entire incidence graph of K_n is constructed at once, and arbitrary G is obtained by deleting the new vertices corresponding to unwanted edges. The stronger construction therefore has n+binom(n,2) vertices, while the restriction to G has n+|E(G)| vertices.

## Construction

Let a be the least integer with 2^a>=64n, and put lambda=7*2^a. Then

    448n <= lambda < 896n,   a>=6.

Give original vertex i in {1,...,n} label A_i=lambda*i+3. For every i<j, put

    h=j-i,
    r_ij=lambda*h/4-1,
    X_ij=r_ij^2-A_i.

These are integers, and r_ij is odd. The two prescribed incident sums satisfy

    A_i+X_ij=r_ij^2,
    A_j+X_ij=(r_ij+2)^2,

because 4r_ij+4=lambda(j-i).

## Positivity and distinctness

Write H=lambda*n+3. For every pair,

    X_ij >= (lambda/4-1)^2-H > H.

Indeed lambda>=448n gives lambda²/16 >=28lambda*n, so the claimed inequality follows after subtracting lambda/2+2H-1. Every new label is positive and larger than every original label.

For a fixed h, distinct i give distinct X_ij because the term subtracted is lambda*i+3. For two gaps h'>h, their squared-root terms differ by at least

    3lambda²/16-lambda/2 > lambda(n-1).

This exceeds the largest possible difference between two original labels. Hence labels from distinct gap classes cannot coincide. Thus all original and all new labels are distinct.

## Remove within-family edges modulo seven

All originals are3 mod7. Since lambda is divisible by7, every r_ij is−1 mod7, and every X_ij is5 mod7. Sums within the two families are respectively6 and3 mod7. Neither is a square residue modulo7. Thus no original-original or new-new edge occurs.

This part does not distinguish which original vertices meet which new vertex. That exact incidence is enforced by a different modulus, below.

## Remove unintended incidences by a dyadic root guard

Suppose A_k+X_ij is a square s², with s>=0. Set r=r_ij and t=s-r. Then

    lambda(k-i)=2rt+t².                                  (1)

Since r is odd and lambda is divisible by8, s is odd, so t is even. Substituting 2r=lambda*h/2-2 into (1) shows

    t(t-2) is divisible by 2^a.

Write t=2z. Since z and z-1 are relatively prime, one of them is divisible by2^(a-2). Therefore

    t = 0 or 2 modulo 2^(a-1).                           (2)

On the other hand, rationalizing the square-root difference gives

    |t| = |lambda(k-i)|/(s+r)
         <= lambda(n-1)/(lambda/4-1) < 5n.

The modulus in (2) is at least32n. The only representatives of its two residue classes in (-5n,5n) are t=0 and t=2. Equation (1) now forces k=i or k=j respectively. Hence X_ij has exactly its intended two neighbors among all original vertices.

Together with the mod-seven nonedges this proves the induced one-subdivision assertion.

## Height and computational size

Each new label is less than r_ij² < lambda²n²/16 <50176 n^4. Original labels obey the same loose bound. Computing all labels takes O(n²) arithmetic operations on O(log n)-bit integers. For a graph supplied by an edge list, only its n+|E(G)| labels need be computed.

## Connection and exact losses

The two-vertex bridge used earlier created a length-three path with three residue families. Here a difference-of-squares factorization with root gap2 creates a length-two path. The new dyadic side condition forces the two allowed roots and excludes every third incidence. The role of the finite moduli is divided: modulo7 prohibits whole families of edges, while a power of2 selects the exact endpoints.

This keeps adjacency of the subdivision, not adjacency of G on its original vertex set. A one-subdivision is always bipartite, so chromatic information and odd-cycle parity of G are destroyed. The earlier two-subdivision construction remains independently useful: replacing edges by length-three paths preserves cycle parity. Neither construction produces Euclidean unit-distance realizations or proves the Hadwiger–Nelson bounds.

In particular the square-sum graph contains induced one-subdivisions of K5, K3,3 and the Petersen graph by the same formula. This demonstrates incidence universality with prescribed subdivision depth, rather than a coincidence among their vertex counts.
