# The globally smallest pure ternary transitive-clone clock

Session 197. **The least positive amplitude of a family obtained by transitive cloning of one fixed tournament vertex and satisfying `H_m=A*3^m-1` is `A=60`.** The conclusion still holds if that identity is required only for all sufficiently large positive integers `m`.

An attaining source has eight vertices, outgoing rows

    [2,12,145,149,99,143,175,19],  root 7.

Replacing its root by a transitive `m`-vertex module gives

    H_m = 60*3^m-1,  m>=1,
    179, 539, 1619, 4859, 14579, ... .

Thus successive clone multiplicities implement the actual returning recurrence `H_(m+1)=3H_m+2`. The attaining construction is from the companion hybrid-construction work. The present theorem adds a complete all-order exclusion of smaller amplitudes. It concerns this precise single-root transitive-cloning mechanism, not every possible graph-changing compiler, and does not solve the unrestricted Hamilton-count spectrum.

## 1. Exact run coordinates and the matching bound

Let `(T,r)` be a rooted tournament. For `k>=1`, let `a_k` count legal words using every vertex of `T-r` once and `k` indistinguishable root symbols, with no two root symbols consecutive. Equivalently there are `k-1` nonempty internal path fragments between roots and two possibly empty exterior fragments. The inherited maximal-run substitution formula, also obtained directly by cutting at clone runs, is

    H_m = sum_k a_k * k! S(m,k).

Indeed, each of the `m` distinct ordered clones is assigned to one of `k` nonempty ordered runs; inside every run the order is forced. There are exactly `k! S(m,k)` such assignments. Every coefficient is nonnegative, and `a_1=H(T)`.

Put `O=N^+(r)` and `I=N^-(r)`. Form the bipartite graph whose edges are the arcs `o -> i` with `o in O, i in I`, and let `nu` be its matching number. Then

    max { k : a_k>0 } = nu+1.                         (1)

For the upper bound, every internal fragment starts in `O` and ends in `I`, so it contains a crossing edge `O -> I`. Different fragments give vertex-disjoint crossing edges.

For the converse, take any matching of size `nu>=1` and use its edges as internal fragments. Extend the first matched edge to a path containing all vertices not used by the other matched edges, preserving the order of its two endpoints. This is possible by repeatedly inserting a new vertex immediately before the first path vertex it beats, or appending it if none exists. In that extended path, cut immediately before the first `O` vertex and immediately after the last `I` vertex. These positions are in the correct order because the original matched edge was preserved as an ordered subsequence. The middle fragment starts in `O` and ends in `I`; the prefix consists only of `I` vertices and the suffix only of `O` vertices. They can be placed before the first root and after the last root. All other matching edges remain nonempty internal fragments. This constructs a legal word with `nu+1` roots. For `nu=0`, no internal fragment is possible, while `a_1>0` by the existence of a Hamilton path. This proves (1).

If `H_m=A*3^m-1` eventually with `A>0`, the largest nonzero run index must therefore be three: the highest exponential term of `k! S(m,k)` is `k^m` with positive coefficient one. Hence `nu=2` and

    H_m = a_3*3^m + (a_2-3a_3)*2^m
          + (a_1-2a_2+3a_3).

Comparison of these distinct exponential terms proves

    a_3=A,  a_2=3A,  a_1=3A-1.                     (2)

In particular the identity then holds for every `m>=1`. Since tournament Hamilton counts are odd, `A` is even. For a source with matching number two, (2) is equivalent to the two finite exact checks

    H_2 = 3H_1+2,    H_3 = 9H_1+8.                (3)

These are a sufficient certificate because the full run bound is proved first; they are not an unsupported fit to three sample values.

## 2. Removal of all irrelevant strong components

Strong components in a tournament form an ordered chain, and Hamilton counts multiply along that chain. If the root is a singleton strong component, transitive cloning leaves the count constant, so it cannot give a positive-amplitude ternary family.

Otherwise, replacing the root by transitive clones keeps its enlarged component strong. Each clone lies on the same root-containing directed triangle with two fixed old vertices. All other strong components remain unchanged. If their product count is `c`, then `c` divides both `H_1=3A-1` and `H_2=9A-1`. But

    gcd(3A-1,9A-1)=gcd(3A-1,2)=1,

because `A` is even. Therefore `c=1`. All those other strong components are singletons, and deleting them preserves the complete clone-count family. Thus a strong rooted source suffices for every possible amplitude.

For a hypothetical smaller amplitude, `A<=58` and the strong source has

    H_1=3A-1<=173.

Busch's strong-order bound gives `5^(n-1)<=H_1^3`. The exact comparison

    5^10=9765625 > 173^3=5177717

therefore forces `n<=10`. This is an all-order reduction; arbitrary singleton padding has already been removed rather than incorrectly bounded by the strong-order theorem.

The external theorem is [Arthur H. Busch, *A Note on the Number of Hamiltonian Paths in Strong Tournaments*, EJC 13 (2006), N3, Theorem 1](https://www.combinatorics.org/ojs/index.php/eljc/article/download/v13i1n3/pdf), DOI `10.37236/1141`. Its precise scope and integer consequence were checked directly from the official paper.

## 3. Complete generation of all bounded strong cores

**Strong deletion lemma.** Every strong tournament of order at least four has a vertex whose deletion is strong.

Take a proper strong vertex subset `S` of maximum cardinality. A directed triangle supplies such a subset. If `|S|=n-1`, there is nothing to prove. Otherwise, a vertex outside `S` with both an in-neighbor and an out-neighbor in `S` would extend it to a larger proper strong subset. Thus every outside vertex either loses to all of `S` or beats all of `S`. Call these two sets `L` and `U` respectively. Strongness forces both sets nonempty and some arc `l -> u` from `L` to `U`. Then `S+l+u` is strong. If there were more than two outside vertices this would contradict maximality. If there are exactly two, deleting any vertex of `S` leaves a strong tournament: its nonempty remaining part satisfies `S-v -> l -> u -> S-v`, irrespective of its internal arcs. This is again a larger proper strong subset, a contradiction. The lemma follows.

It is therefore complete to start from the directed triangle and repeatedly add one vertex in every mixed orientation to the old strong tournament. Every such extension is strong, and every strong tournament arises in this way by reversing a deletion chain.

The Hamilton cap is also safe at every intermediate order. Insertion of a new vertex into each old Hamilton path, preserving its old order, injects old paths into new ones. Hence Hamilton counts are nondecreasing under induced extensions. Every deletion of a final count-at-most-173 graph still has count at most 173.

Representatives are deduplicated by sorting vertices by outdegree and minimizing their adjacency tuples over every permutation within equal-outdegree groups. This is exact isomorphism canonicalization because every isomorphism preserves outdegrees. No approximate graph invariant or random sampling is used.

The resulting complete strong inventories are:

| Order | Extension candidates | Accepted before deduplication | Classes with H<=173 |
|---|---:|---:|---:|
| 3 | initial triangle | initial triangle | 1 |
| 4 | 6 | 6 | 1 |
| 5 | 14 | 14 | 6 |
| 6 | 180 | 180 | 35 |
| 7 | 2170 | 2168 | 351 |
| 8 | 44226 | 3756 | 659 |
| 9 | 167386 | 1833 | 347 |
| 10 | 176970 | 356 | 66 |

Every vertex is tested in every retained graph with `H_1=5 modulo6`, the necessary residue for even `A`. Exactly **1624** such rooted candidates have matching number two. None satisfies both equations (3). By the reductions above, no smaller amplitude is possible at any source order.

## 4. The attaining clock and verification

For the displayed eight-vertex source, direct matching and exact clone counts give

    nu=2,  (a_1,a_2,a_3)=(179,180,60).

The expansion in Section 1 then proves `H_m=60*3^m-1` for all positive `m`, attaining the lower bound. Actual graphs for multiplicities one through five provide additional direct controls.

`verify_minimum_ternary.py` compiles and runs `verify_minimum_ternary.cpp` in a temporary directory. The C++ search exhausts every mixed extension under the safe cap and emits a deterministic complete representative bank. The Python layer independently checks every retained graph's strongness, exact canonical representative and Hamilton count, and independently retests all candidate roots using augmenting-path matching and a backward endpoint recurrence. The C++ search uses subset matching and a forward recurrence. The Python layer also enumerates every labelled tournament through order six to check the small strong banks independently of the growth construction.

Run normally and under `python -O`, and optionally with `TERNARY_CPP_MODE=debug` to use an unoptimized C++ build with library checks. Outputs are deterministic. Explicit exceptions, rather than Python assertions, enforce all checks. The receipt reports the complete bank hash and each control count. A C++17 compiler is required.

This is a globally sharp theorem for one specified physical clone operation. It does not assert that amplitude 60 is minimal for unrelated ternary compilers or that the new family covers every prime or `7p` target.
