# Primary-source positive K4,4 construction

This is a separate addition; the bounded exclusion through N=10000 remains unchanged.

## Attribution and scope

Ajai Choudhry, *Two sets of integers such that all elements of the sumset of the two sets are perfect squares*, arXiv:2508.07806v1 (11 August 2025), https://arxiv.org/html/2508.07806v1, introduction and equation (1.1), reproduces a numerical example from Guy and a two-parameter family attributed to Jean Lagrange through Guy. Consequently bare integer K4,4 existence is established. The code evaluates that published family; the parameter formulas are not new research here. Our normalization, bounded enumeration, gauge optimization, and exact checks below are independent computations. Neither the source construction nor our witness claims a magic square of squares.

## Exact positive witness and bounds

The published numerical example, shifted to positive disjoint parts, gives

    L = {323,54755,119715,264195}
    R = {1,79201,227206,1258561}.

The accompanying JSON verifies all sixteen sums by their exact integer square roots. Reflection defects, defined for a sorted four-set X by X[0]+X[3]-X[1]-X[2], are 90048 and 952155, respectively. Both are nonzero.

Our complete bounded scan of primitive integer parameters -200<=u,v<=200 in equation (1.1) checks 97856 pairs, of which 97780 yield positive disjoint normalized parts. Its least achievable maximum label is 1112833, attained at (u,v)=(-11,5), with

    L = {8648,44943,242063,725903}
    R = {1,62641,242353,1112833}.
    Roots = [[93,267,501,1059],
             [212,328,536,1076],
             [492,552,696,1164],
             [852,888,984,1356]].

All eight vertices are positive and distinct. The reflection defects are 447545 and 807840. Thus neither side has the equal opposite-pair sums required for reflection/parallelogram decoration, under any ordering. These defects are unchanged by gauge translation; scaling cannot change whether they vanish. In particular this is a concrete bare biclique, not a certificate of the decorated/MSOS condition.

Combining the independent bounded exclusion and this witness gives

    10001 <= least N containing an ordinary K4,4 <= 1112833.

The upper endpoint is not claimed optimal globally, nor over the entire unbounded parametric family. It is only the optimum in the precisely specified primitive parameter box after the following normalization and gauge procedure.

## Normalization and exact optimal gauge

Compute all sixteen roots of a_i+b_j and their gcd h. Let b0=min(B). Then

    A'_i=(a_i+b0)/h^2, B'_j=(b_j-b0)/h^2.

These quantities are integral: every sum is divisible by h^2, and every difference of two such sums is too. This normalizes the common root gcd to 1 and makes min(B')=0. Reject duplicate labels within either side. A positive realization is L=A'-t, R=B'+t, with integer 1<=t<=min(A')-1. Disjointness forbids at most sixteen integers t, namely those satisfying 2t=a'_i-b'_j. The objective is

    max(max(A')-t,max(B')+t).

It is V-shaped with real minimizer (max(A')-max(B'))/2. Clamp its floor to the feasible interval; testing distance at most 17 around that point, together with the adjacent ceiling, contains an optimal permitted integer because at most sixteen integers are forbidden. If none are permitted, reject that member. The implementation explicitly verifies positivity, pairwise distinctness, and all sixteen square sums.

## Reproduction

`python primary_construction.py` regenerates primary_construction_result.json for parameter bound 200. Explicit exceptions, rather than assertions, guard all checks, so optimization mode does not suppress them. The result records the original numerical example, bounded-family witness, root matrices, normalization, gauge, counts, and nonzero decoration defects. No original repository files or remotes were modified.
