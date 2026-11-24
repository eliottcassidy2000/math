# The first blue line, and actual routes to compact Moore classes

## Work in the native metagraph, not the edit cube

A class node is an ordinary, unlabelled tournament isomorphism class. An incidence is represented by a chosen directed Hamilton path P, relabelling its vertices in path order. A native step complements **all** nonspine tiles in that chart. It is blue if the tile mask is invariant under reflecting the spine. Different paths in the same tournament may give different lines, and a class can have many representatives. These are not one-coordinate edits of a single rooted mask.

The transitive class has only one rooted mask, zero. It has exactly one native line, blue, to the all-ones mask. On the chosen spine that neighbor has arcs i→i+1, all other arcs going backward. Its directed Hamilton paths correspond to compositions of n with no adjacent singleton parts: write the maximal increasing consecutive runs as intervals, which must occur in decreasing interval order. Two neighboring singleton intervals are the only obstructed boundary. Thus the Hamilton counts have generating function

    (1−z²)/(1−z−z²−z³).

For n=3,...,10 they are 3,5,9,17,31,57,105,193. Everything after that first line may branch and even have cycles or coincidences; this is not a tree assumption.

## A lower bound on path length, valid even with black lines

Let R_P reverse the arcs on the path only. A native move is the converse of R_P. Under R_P only two outdegrees change, by −1 at the start and +1 at the end. If sorted scores are s_0,...,s_(n−1), then

    D(T) = 1/2 Σ_i |s_i−i|

is unchanged by converse and changes by at most one per native line (relabeling and sorting do not enlarge matching L1 distance). It is zero at the transitive class. Thus D is a lower bound on native distance, unrestricted and blue-only.

For the original doubled observer at degree d, m=d²+1, the score lists are m+d−i and m−d−1+i, i=0,...,m−1. For d≥3 the shared score interval in the two lists has length 2d+2; each tail moves by d+1 and the lower half of the overlap contributes (d+1)². The result is D=(d+1)(m−d−1), also directly true at d=2. The old observer lower bounds are therefore **6,24,336,185136** at n=10,20,100,6500. Degree 57 remains conditional. The compact HS tournament has lower bound **106** at n=50.

## Exact shortest all-blue paths

Identifiers below are the minimum mask over **all** Hamiltonian relabelings of a class, not only its blue presentations. Bits are pairs (i,j), j≥i+2, in lexicographic order, with 1 indicating the arc j→i. Each line lists successive *source* masks after the required reroot; complement against all ones at the given size. The actual Hamiltonian reroot orders are in [native-replay.py](native-replay.py) and are independently checked against the previous output, not assumed from the height.

| Target | n | Class ID | Hamilton paths | Shortest native distance | Source-mask witness |
|---|---:|---:|---:|---:|---|
| Regular pentagon | 5 | 22 | 15 | 3 | 0,26,18 |
| Original doubled pentagon, order 01423 | 10 | 27588026 | 6473 | 6 | 0,35437707265,1325379632,24640,41970477573,9515068760 |
| Compact Petersen | 10 | 368700356 | 8265 | 7 | 0,5369004036,7520420484,12313403264,24768,10072225120,19847149298 |

Every incidence displayed is blue; the score bound proves each distance is shortest even among all lines. For compact HS, the [packed replay](hs-prefix-check.py) gives 91 consecutive blue moves *out of* the target reducing its debt from 106 to 15. It does **not** reach transitive. No complete native path to compact HS or the larger 20/100/6500 observer charts is proved here.

## Many classes already share the first neighbor's scores

For n≥4, look at the precise sorted score multiset of the first blue neighbor: 1,1,2,...,n−3,n−2,n−2 (for n=4 take 1,1,2,2). Order vertices by increasing score. For every proper prefix of length k there is exactly one arc from the prefix to its complement, since its total wins exceed k(k−1)/2 by one. Consequently every lower-to-higher arc determines an interval of cuts, and these intervals partition all n−1 cuts. The exceptional arcs therefore form a single increasing path on some breakpoints including both extremes. All other arcs go from higher to lower index.

The two tied pairs at the ends can be swapped independently, toggling the first and last interior breakpoints; all strictly intermediate score positions are intrinsic. Hence exactly **2^(n−4)** ordinary classes share that score multiset. The transitive class is adjacent to just one of them. Converse reflects the central breakpoint set, giving **2^ceil((n−4)/2)** self-converse classes in this slice. These are precisely backward-key strips of the [Chan–Pak transfer framework](../chanpak_tournament_sl2_20261124/): the spacers between successive breakpoints are the strip transfer inputs.

At n=5 and n=6 the strip with only the two extreme keys has the same H as the first blue neighbor, 9 and 17 respectively, as well as the same scores. Its H is 1+2^(n−2) by choosing which interior vertices precede the unique ascent. Yet its native distance is **two**, not one (finite replay). Even the pair of observables (Hamilton count, sorted scores) can therefore hide which class the first blue line actually enters.
