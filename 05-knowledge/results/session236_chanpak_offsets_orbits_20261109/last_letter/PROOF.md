# All-word tower-ray rigidity for the complete 57-matrix alphabet

**Status:** proof candidate, independently audited before final promotion.
**Scope:** the exact ternary alphabet of all marked tournament blocks of orders1 through4, not arbitrary blocks and not orbit cofiniteness.

Let G be the57 distinct nonnegative integer matrices in `matrices_base3_order4.txt`, ordered lexicographically. They are exactly the THM-5070 matrices N_3(B,w,u') for all1105 marked labelled tournament blocks of orders1 through4. The file also records one block representative per matrix in `.reps`. Define

    T = [[1,3,0],[0,0,1],[1,1,0]],     v=(3,1,2),     s0=(1,0,1),
    C = {(h,E,O): h>0, E>=0, O>=0, h<=3O}.

## LL1. All-word rigidity

For every finite word P=M_r ... M_1 over G, the following are equivalent:

1. Pv is a positive multiple of v.
2. Every letter M_i is a power T^k with k in {1,2,3,4}.
3. P=T^K for an integer K>=0 (the empty word gives K=0).

Consequently no word over G sends s0 to the ray R_{>0}v.

This proves an **all-word** replacement for THM-5070(5)'s finite check, whose inherited scope was length at most4. It does not prove or refute HYP-9640(A)'s first-coordinate coverage.

### Certificate 1: an invariant cone

For every M in G, its (1,1) entry is positive, all its entries are nonnegative, and the row

    3 row_O(M) - row_h(M)

is entrywise nonnegative. Thus M maps the entire nonnegative h>0 cone into C. Both v and s0 lie in C.

This inequality also has a direct physical explanation, valid for every tournament and every real x>=1. Put D=I(Omega(T-u),x). Classifying odd-cycle packings according to whether they cover u gives

    h = D + x * sum_{odd cycles C through u} I(Omega(T-V(C)),x).

Break each such cycle immediately before u. It becomes a distinct odd-order path starting at u, hence contributes to O; the singleton path contributes D. Therefore

    h <= x O - (x-1)D <= x O.

The displayed matrix-row certificate is stronger for the specific alphabet: it holds even on abstract nonnegative inputs, without assuming physical realizability.

### Certificate 2: classify the possible last letters

If Ms lies on the target ray, then

    A_M s = 0, where A_M has rows row_h(M)-3row_E(M), row_O(M)-2row_E(M).

The57 matrices split exactly as follows.

- Four matrices are T,T^2,T^3,T^4, at zero-based indexes1,9,29,50. They are invertible and Tv=2v.
- Fifty-two matrices have an integer linear combination of the two rows of A_M of the form (a,b,c), with a>0 and b,c>=0. Since h>0 and E,O>=0, that combination cannot annihilate s. The full52 certificates are recorded in CHECK.json. Their coefficients have absolute value at most3.
- One exceptional matrix is

        M_* = [[7,18,12],[2,6,7],[5,7,6]].

  Its equations are h-9O=0 and h-5E-8O=0. Their nonnegative solutions with h>0 are precisely positive multiples of (45,1,5). These are outside C because45>3*5.

Thus, for s in C, Ms on the target ray forces M=T^k. Inverting gives s on that ray as well.

### Induction on the last letter

Start with v. Every intermediate vector remains in C. If the final vector lies on its ray, Certificate2 forces the last letter to be a tower power and the previous vector to lie on the ray. Repeat backwards through the entire word. This proves1 implies2. Multiplying powers gives2 implies3, and Tv=2v proves3 implies1.

The same backwards argument from s0 would force s0 onto the ray, which is false. So the ray is never hit by the physical singleton-seeded orbit. Avoiding that ray alone does not constrain the orbit's possible first coordinates.

## Why the invariant cone matters

The unqualified claim "every non-tower matrix has no nonnegative preimage of the target ray" is false for the57-matrix alphabet:

    M_*(45,1,5)=131(3,1,2).

The cone eliminates precisely this surviving exceptional preimage. This is a positive hostile control, not a failed computation.

## Verification and mathematical status

Run:

    python3 05-knowledge/results/session236_chanpak_offsets_orbits_20261109/last_letter/verify.py --rebuild
    python3 -O 05-knowledge/results/session236_chanpak_offsets_orbits_20261109/last_letter/verify.py --rebuild

The rebuild enumerates all1105 marked labelled blocks using the inherited THM-5070 transfer. The checker then uses only integer arithmetic to verify all57 cone rows, all52 separating forms, all four tower identities and the unique exceptional solution. Normal and optimized runs agree. `--mutate` corrupts a generator and must fail.

The certificate list is finite-exact; **the theorem quantifies over every word length** because the finite local statements feed a mathematical induction, not a longer brute-force search. The proof assumes only the frozen alphabet and the cited transfer law. A separate independent reconstruction checks the block interpretation at audit.

## Precise consequence for the requested offset strategy

Graph offsets preserve d-2^(n-1). A tournament digit proof that insists all its block products preserve the ternary tower direction cannot introduce any non-tower word in this alphabet. Only the old tower remains. Adaptive states, other directions, larger blocks and other tournament operations are not excluded. The binary56 problem is a different orbit and is not decided by this ternary theorem.
