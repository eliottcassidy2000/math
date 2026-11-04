# Exact full-response separation in the three-state gate's edit neighborhoods

**FINITE-EXACT:** no homogeneous-pair-free tournament obtained by reversing one, two, or three edges of any of the nine specified S185 gate inflations has the same named full two-port response as its source. Three prime scalar collisions are exhibited and rejected by a single literal response query. This is not universal reconstruction, an unbounded neighborhood theorem, or a full-response collision.

## Observer and the quartic boundary

All vertex names stay fixed. Attach x,y with x->y, with a_v=1 meaning v->x and c_v=1 meaning v->y. Let F_T(a,c) be the actual Hamilton-path count. The unique Boolean multilinear polynomial has bidegree at most (2,2), hence total degree at most four. No quartic coefficient may be discarded when comparing arbitrary T and U. The degree-three comparison battery of `two_port_twelve_principles_20261103/PROOF.md`, T4, applies to T versus its literal converse, not to arbitrary pairs with different module structures.

The known current-main reconstruction census already reaches order ten. Our independent labelled orders2–6 check is only a control, not a new reconstruction bound.

## Defined search class

Use S185's seeds

    Q0=(16,5,33,39,46,3),
    Q1=(16,5,49,39,42,3),
    Q2=(18,4,49,39,42,3).

Replace vertex5 by the increasing transitive module on names5,...,4+r, for r=6,7,8. Let the result be T_j(r), of order11,12,13 respectively. S185 proves

    H(T_j(r))=4*3^r+9*2^r-3+2j.

For each of these nine *named* sources, enumerate every subset of1,2,or3 distinct old edges and reverse precisely those edges. Retain targets with no homogeneous pair. This necessary prime-tournament filter is broader than being prime: no assumption identifies pair-free with prime.

There are464709 source/edit-set cases. Exactly28293 resulting targets have no homogeneous pair. The exact Hamilton scalar differs from its source for28290 of those. Since F_T(0,0)=H(T), these are already separated by a full-response query. The remaining three targets are below; all three are prime and strong, and all have a different singleton response. Therefore **every target in this defined class is full-response separated from its source**, without assuming that H,S,E,C determine full F.

## Three prime scalar hostiles

All three occur at r=6,j=0. The source is

    T=(16,5,2017,2023,2030,1987,1923,1795,1539,1027,3),
    H(T)=3489.

It has the five homogeneous pairs {5,6},{6,7},{7,8},{8,9},{9,10}.
Reverse the following three-edge sets:

1. {0,6},{0,9},{6,8}, giving
   U1=(592,5,2017,2023,2030,1987,1666,1795,1603,1026,3).
2. {0,7},{0,9},{6,8}, giving
   U2=(656,5,2017,2023,2030,1987,1667,1794,1603,1026,3).
3. {0,8},{0,9},{6,8}, giving
   U3=(784,5,2017,2023,2030,1987,1667,1795,1602,1026,3).

Every Ui is strong, has no proper module of size2,...,10, and satisfies H(Ui)=3489. Thus these are genuine homogeneous-pair/prime scalar collisions, with fixed physical labels. Nevertheless, the actual attachment query a={0},c=empty gives

    F_T({0},empty)=17055,
    F_Ui({0},empty)=15579.

Since F(0,0)=3489 on both sides, the named linear coefficient of a_0 is13566 for T and12090 for every Ui. The difference1476 is already visible before the higher coefficients are examined. No profile sorting, name permutation, response normalization, or dropped quartic coefficient is involved.

## Exact algorithms and checks

`perturb.cpp` builds each adjacency matrix literally, enumerates edge subsets in increasing index order, tests all vertex pairs for homogeneity, and uses subset endpoint DP for H. Its extra E,S filters are diagnostics only: the proof above needs just its complete scalar-match list and the independent literal queries. Unsigned64-bit arithmetic suffices because every computed old-graph count is at most13!<2^64.

`verify.py` independently rebuilds the named source, checks the three exact edge differences, enumerates every proper vertex subset to verify primality, tests directed reachability, and computes H and F by predecessor-based subset DP. Its checks use exceptions rather than removable assertions. Normal and optimized receipts agree. The largest literal attachment has13 vertices, so there is no floating point arithmetic.

`search.cpp` is a separate small-universe complete coefficient control. For every permutation of the n+2 named vertices, it rejects backward old-old steps and the forbidden y->x step. Every remaining old-port edge contributes its literal Boolean factor a_v,1-a_v,c_v,or1-c_v. Expanding and summing gives every coefficient, including quartics. All (1+n+binomial(n,2))^2 possible bidegree-(2,2) coordinates are retained. Equality is tested by entire integer vectors, not hash equality. The numbers of distinct responses are2,8,64,1024,32768 at n=2,...,6, exactly the numbers of labelled tournaments. At n=6 there are484 retained coordinates.

`degrees.cpp` additionally observes that truncation through degree two already distinguishes this small universe. This is a bounded diagnostic only, not permission to discard cubic or quartic terms elsewhere.

## Reproduction and limitations

    g++ -std=c++17 -O3 -Wall -Wextra -Werror perturb.cpp -o perturb
    ./perturb
    python verify.py
    python -O verify.py
    g++ -std=c++17 -O3 -Wall -Wextra -Werror search.cpp -o search
    ./search

This packet does not decide four-edge neighborhoods, arbitrary names/isomorphism alignments of the order11–13 sources, arbitrary order-eleven reconstruction, or hidden module-converse cores of size at least twelve. It does not claim that the scalar hostiles are new in the literature. The useful exact boundary is: the recent scalar three-state gate produces actual prime scalar collisions under very small perturbations, but none in this prescribed search survives the full named response.
