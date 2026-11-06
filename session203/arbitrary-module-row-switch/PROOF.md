# An arbitrary-module pair with the same padding row

Status: **PROVED, independently audited in the stated scope.** This is a positive all-order equal-profile construction. It does not yet widen a fiber beyond two choices or prove eventual odd-count coverage.

## Physical construction and statement

All rows below encode outgoing edges, with vertices numbered from 0. Start from the rooted tournaments

    A = (102,0,10,83,135,30,182,47), root 0,
    C = (142,0,34,22,71,155,47,94), root 0.

For any tournament X, replace vertex 5 of A, or vertex 3 of C, by a copy of X. Every outside vertex has the same direction toward all vertices of X that it had toward the replaced vertex. Retain all other edges. Write A[X], C[X] for these actual tournaments. Empty X means deleting the designated vertex. The designated root is never replaced. Both contexts retain a global sink, so no strong-output claim is made.

For a rooted tournament U with root r, define:

- h=H(U), s=the number of Hamilton paths starting at r, e=the number ending at r;
- d=H(U-r);
- R=the number of ordered spanning two-path covers (P,Q), where P may be empty and Q starts at r;
- L=the number of ordered spanning two-path covers (P,Q), where P ends at r and Q may be empty;
- beta=the number of ordered spanning two-path covers of U-r, with either path allowed empty;
- B=h+L-s and V=R+beta.

The padding row is (B,e,V,d). For j>=1 let p_j(X) count ordered covers of V(X) by j nonempty directed paths. Put p_0(X)=1 for empty X and 0 otherwise. Equal path slots remain labeled, so p_j is divisible by j!.

**Theorem.** For every finite tournament X, including the empty tournament, A[X] and C[X] have exactly the same padding row. Their complete profile differences are

    Delta(h,s,e,L,R,d,beta)
       = (2(p_1+p_2), 2(p_1+p_2), 0,0,0,0,0).       (1)

The common row is

    B = 93p_0 + 253p_1 + 204p_2 + 48p_3,
    e = 0,
    V = 300p_0 + 1022p_1 + 1178p_2 + 520p_3 + 72p_4,
    d = 11p_0 + 29p_1 + 21p_2 + 4p_3.                (2)

For A[X], the remaining profile coordinates are

    s = 5p_0 + 16p_1 + 15p_2 + 4p_3,
    h = 37p_0 + 111p_1 + 95p_2 + 24p_3,
    L = 61p_0 + 158p_1 + 124p_2 + 28p_3,
    R = 100p_0 + 382p_1 + 486p_2 + 236p_3 + 36p_4,
    beta = 200p_0 + 640p_1 + 692p_2 + 284p_3 + 36p_4. (3)

For C[X], add 2(p_1+p_2) to h and s in (3), leaving the other coordinates unchanged.

For a singleton X this recovers row (253,0,1022,29), with s=16 or 18. For a transitive X of size m>=1, p_1=1 and p_2=2^m-2, so the endpoint gap is 2^(m+1)-2. No transitivity assumption is needed for (1)-(3).

## Exact core-gap proof

For an arbitrary module X, a directed spanning path in a fixed surrounding core has a unique ordering of the core vertices and a unique set of occupied gaps. Each occupied gap contains one nonempty path of X. An empty interior gap requires a forward core edge. An occupied gap requires its left endpoint, if any, to beat X and its right endpoint, if any, to lose to X. Thus the number of skeletons with j occupied gaps is the coefficient of p_j(X). This bijection works for all internal orientations of X.

The coefficients for h come from summing these gap polynomials over all core orders. To compute s, require the first core vertex to be the root and its initial gap to be empty. To compute e, require the last core vertex to be the root and its final gap to be empty. Deleting the root gives d in the same manner.

The cover coordinates admit the same exhaustive gap calculation using auxiliary directed graphs. These graphs are counting devices, not physical tournament outputs:

- For R, merge an auxiliary separator immediately before the root into that root. Make every other core vertex point into the merged root, and make every module vertex point into it. Retain the root's original outgoing edges. A path in this possibly bidirected graph splits uniquely as P followed by the original root-starting path Q, giving R.
- For L, merge a separator immediately after the root. Retain all original incoming edges and allow the merged root to beat every other vertex. Paths split uniquely into the original root-ending path P and Q, giving L.
- For beta, first delete the root and then add an auxiliary separator with edges in both directions to every other vertex, including the module. Splitting at the separator gives all ordered two-path covers, with either side allowed empty.

The same gap bijection applies to these auxiliary directed cores: an empty gap tests the actual directed edge, and an occupied gap independently tests its incoming and outgoing permissions. No tournament parity theorem is applied to an auxiliary graph. All physical outputs A[X] and C[X] remain tournaments.

The standalone module_tensor.py enumerates all core orders for these seven profile coordinates in each context and outputs the complete integer coefficient arrays in MODULE_TENSOR.json. Direct subtraction gives exactly (1), and addition gives (2). The calculation exhausts only fixed finite cores, so these are exact identities for arbitrary X, not interpolation from small module counts.

A separate C++ implementation, verify_module.cpp, computes literal Hamilton-path subset tables and cover convolutions for every labeled X through order 4, including empty X. It checks all coordinates against (1)-(3) on 152 actual context graphs. Optimized and unoptimized builds give identical receipts. These literal tests are controls of the coefficient certificate, not its all-order justification.

## What this does and does not provide

The pair preserves the actual padding state for every X and allows two different marked endpoint values when X is nonempty, with no restriction on its internal structure. Common O/Z suffixes therefore preserve the collision. The standard physical additive context H=3B+2s gives two counts differing by 4(p_1+p_2), but no iteration or interval coverage is inferred.

The row depends only on the first four path-cover moments of X. A substitution preserving all four of those moments also preserves both endpoint values in (3), and therefore does not automatically add a new independent choice. Nesting the pair itself changes the unrooted path-cover moments, so equal-row status is not by itself a closure theorem.

There is a precise weaker moment condition. For two nonempty module inputs, preserving B,V,d imposes three independent linear equations on their four moment differences. Their rational kernel is generated by

    Delta(p_1,p_2,p_3,p_4)=(-1728,3420,-5427,7768).     (4)

Along (4), the A endpoint changes by 1944 and the C endpoint by 5328. Requiring each p_j to be divisible by j! restricts integral multiples in (4) to multiples of 6; the primitive corresponding direction in unordered q_j=p_j/j! is (-10368,10260,-5427,1942). These are arithmetic necessary conditions for actual modules, not a realizability claim. Finding physical module inputs along this kernel, or a richer retained interface with independently composable switches, remains open.

## Reproduction and provenance

    python3 module_tensor.py > tensor_replay.json
    cmp tensor_replay.json MODULE_TENSOR.json
    g++ -O3 -std=c++17 verify_module.cpp -o /tmp/module_pair
    /tmp/module_pair > replay.json
    cmp replay.json MODULE_NORMAL.json
    g++ -O0 -std=c++17 verify_module.cpp -o /tmp/module_pair_debug
    /tmp/module_pair_debug > debug_replay.json
    cmp debug_replay.json MODULE_DEBUG.json

The finite starting collision came from regrouping the existing rooted tournament profiles, not a new tournament enumeration. The discovery probes are separate from this standalone certificate. The proof requires no completeness assumption for the starting bank: the displayed two cores and their exact identities suffice.
