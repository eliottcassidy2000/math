# Signed motion with an internally growing core

**Proved on the stated physical families.** The goal is not to preserve each large coordinate. It is to preserve the combination that enters the chosen count, making the high drift cancel and leaving a small signed residual. This does not establish an absolute-target selector.

## One contact mask, exact readout

Use the [previous strong head construction](https://github.com/eliottcassidy2000/math/blob/6d38e8ceb0ccc8a061270e52507820bd4ecf53e6/05-knowledge/results/chanpak_computational_tapes_20261115/HEAD.md). The old tape consists of TT_a → D → TT_b followed by a sentinel t. The new control x beats the tape except for V→x inside D, and t→x. Both a,b≥1, so the output is strong. Write (h,q,E,L,J) for H(D), the optional two-slot count, the count with slot one ending in V, the count with slot two beginning in V, and the common intersection, respectively. With X=2^a,Y=2^b,

    H = qXY + EX − LY + h−J.

Therefore with guards a=k+u,b=k+v and k≥1, the two large coefficients are only q·2^(u+v) and E·2^u−L·2^v. Equality of those combinations leaves the offset h−J entirely available. When comparing states, it is essential not to discard h or J individually, or to demand that E and L separately match.

### All five raw coordinates can move in a three-count window

Here labels use lexicographic forward-arc code and masks are sums of powers of two of vertex labels. The following states produce the three consecutive odds `128·4^k−63·2^k−3, −1, +1` (at k=1: 383,385,387).

| n | code | mask | u | v | h | q | E | L | J |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
|4|0|9|0|3|1|16|9|9|4|
|5|4|24|0|1|5|64|15|39|6|
|4|63|1|0|3|1|16|1|8|0|

Every pair of consecutive states changes all five raw coordinates. The outer control and sentinel use the same rule in every state. Holding total order fixed is not needed here; the separate one-core packet supplies fixed-order examples.

### A five-count window changes the quadratic response in mid-route

These five states have common high coefficients (1344,−318); the offsets are −59,−57,−55,−53,−51. The same assertion holds at every k≥1. The middle state changes q, E and L together, then the route returns to the original q without losing a step.

| n | code | mask | u | v | h | q | E | L | J |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
|6|10|45|1|2|9|168|99|129|68|
|6|10|43|1|2|9|168|95|127|66|
|7|2097141|82|0|2|9|336|198|129|64|
|6|10|51|1|2|9|168|91|125|62|
|6|8|53|1|2|9|168|89|124|60|

A “window” here is a list of physical witnesses, not a claim that each transition is one arc reversal, or that consecutive windows at different k overlap. For a fixed finite core bank, even allowing all a,b does not fix absolute placement: the t→x case alone contains at least half of q(D)2^(a+b) by placing the first guard vertex in slot two, so there are only O((log N)^2) candidate witnesses with H≤N per fixed core and mask. Letting the growing module inside the core change takes us beyond finite fixed-core banks, but the explicit one-parameter families below are still sparse.

## The high drift can grow *inside* the same strong core

Let the template D have order 6, code 22. At vertex 2 substitute any nonempty tournament P, denoting the strong result D(P). The vertex labels external to P are renumbered 0 through 4 in original order (old vertices 0,1,3,4,5); vertices inside P are appended. Use masks `4` and `17` on those outside labels, i.e. old vertex {3} and old vertices {0,5}. Use the left/right guard lengths (1,2) and (2,1), respectively. The **same** strong core D(P) occurs in both witnesses, with the same total order. No additional nontransitive block is involved.

Let (p1,p2,p3)=(Q1(P),Q2(P),Q3(P)). The actual profile coordinates, written as coefficient triples on (p1,p2,p3), are:

| coordinate | mask 4 | mask 17 |
|---|---|---|
| h | (3,6,0) | (3,6,0) |
| q | (−6,48,42) | (−6,48,42) |
| E | (12,14,0) | (6,14,0) |
| L | (−1,6,6) | (−4,24,12) |
| J | (0,0,0) | (2,4,0) |

Substituting gives the pair

    H_− = 312 Q3(P) + 394 Q2(P) − 17 H(P),
    H_+ = 312 Q3(P) + 394 Q2(P) − 15 H(P).

Thus moving between masks changes L and J by terms that grow without bound on transitive P, but they cancel *exactly* in the final difference. For P=TT_m, both outputs are strong (the core itself grows), and H_−=312·3^m+394·2^m−17 while H_+=H_−+2. The first four pairs are 1707/1709, 4367/4369, 11559/11561, 31559/31561. This is not uniform multiplication by a spectator block: the high terms read higher slot covers of the module inserted **inside** the strongly connected template.

### Why the finite check proves the polynomial identity for all P

For each of h,q,E,L,J, partition the one- or two-slot cover by the exterior vertex words and the maximal visits to P in each slot. Contacts with P are uniform, and all marked restrictions concern exterior vertices. Every return to P in a slot needs two distinct exterior vertices; with five exterior vertices and at most two slots the total number of visits is at most 4. For any fixed exterior pattern with r≥1 visits its factor is P_r(P), the ordered nonempty r-path cover count. One scalar linear combination of P1,…,P4 therefore suffices, with coefficients independent of P. The change of basis to Q1,…,Q4 is invertible by Q_j=Σ_r binom(j,r)P_r. At P=TT_m, Q_j=j^m. The matrix (j^m) for m,j=1,2,3,4 is a Vandermonde matrix with each column multiplied by j and is nonsingular. Thus four independently enumerated exact transitive substitutions certify the table for every tournament P, not merely the four sampled sizes. For additional controls the replay uses genuinely nontransitive P and counts the final witnesses independently.

The exact issue remains repeatability and absolute selection. When P is arbitrary and H(P)>1 the gap is still amplified by H(P); the two high coordinates of P are common to the pair. A suggested next target is to make two such signed kernels on a single live profile, with the lower term of one canceling the H(P) factor of the other. A second option is to let P vary between the two witnesses and solve for its actual Q2,Q3,H coordinates jointly, rather than holding all three fixed.
