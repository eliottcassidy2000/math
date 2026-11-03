# Independent audit of the multi-birth symmetry budget

**PASS.** The proof in `BOUNDARY_BUDGET.md` is valid for its stated directed
and undirected graph scopes. The prime/7p consequence is necessary, and the
initial sharp family is correctly scoped to unrestricted tournaments. The
subsequent B4 refinement proves strong sharpness for r>=1, k>=2r.

The load-bearing action was checked explicitly. A source group supported on
M_i preserves M_i setwise because it fixes its complement pointwise. On a
path stratum with exact touched-block set S, every old endpoint of an edge
incident to a new vertex lies in a block of S or outside all supported
blocks. The product over i outside S therefore fixes those endpoints.
Old-old edges survive because the product consists of source automorphisms.
The touched set remains S, and the action on full vertex orders is free.
Thus the action need not extend to a graph automorphism of the enlarged
source. There is no missing global attachment-invariance assumption.

Prime-wise minimization of the stratum divisors may remove the largest
2r weights independently for each prime. This correctly computes the
valuation of the gcd; it does not assert that one deletion set minimizes
all prime valuations simultaneously. Summing strata preserves their common
divisor even when some strata are empty.

For independently rotating triangle modules this gives
3^max(k-2r,0) | H. Every target p or7p with p!=3 is prime to3, so r>=ceil(k/2)
is necessary. This conclusion depends on retaining the old induced source;
reorienting an old arc may destroy the supplied group actions and must not
be treated as a legal attachment in this theorem.

## Independent literal tests

`boundary_independent.py` imports no producer or repository code. It builds
arcs as ordered pairs, enumerates Hamilton words by direct path recursion,
and applies explicit rotations to every untouched triangle along each
word. The controls comprise:

* all8 attachments to one triangle;
* all64 attachments to two ordered triangles;
* eight declared attachments to the source C3[C3,C3,C3], which is not an
  ordered join of its three modules.

Across these80 sources/cuts it checks69,342 literal Hamilton words and
104,316 generator images, including preservation of the touched stratum,
nontriviality of each generator action, and all stratum/global divisors.
The mixed two-triangle cut{0,3} explicitly breaks the first source rotation
as a graph automorphism of the extension, yet its untouched-stratum actions
remain valid. This is a hostile to the stronger assumption the theorem
correctly avoids.

The H41 seed is counted literally. Its triangle blocks are (0,1,2) and
(3,4,5), and the new vertex6 beats exactly0 and3. The script independently
builds joins of r such seeds and k-2r triangles for
(k,r)=(1,0),(2,1),(3,1),(4,2),(5,2), and checks every retained old arc after
deleting the explicitly named new vertices. The old source is exactly the
ordered join of k triangles. Ordered-SCC multiplication then gives
41^r 3^(k-2r); its3-adic valuation is k-2r. At this stage the independent test uses ordered-join multiplication, not
an independent enumeration of all17-vertex paths. The later B4 controls
below add an independent subset-DP computation at that order.

Run:

```sh
python3 -I -B boundary_independent.py
python3 -I -B -O boundary_independent.py
```

Both runs agree and all checks use explicit exceptions. No all-prime
realization or stronger strong-target sharpness is inferred.

## B4 continuation: higher boundary divisor and strong sharpness

The new-edge congruence is valid: every affected path contains the changed
new-new edge, which consumes two new-vertex incidences and leaves at most
2r-2 old-new incidences. The untouched-group action fixes the new endpoints
and preserves the condition of containing that edge, so each affected-path
set has divisor D_(r-1). Cancellation of unaffected paths is literal.

The strong-sharpness construction keeps all old-new edges fixed and changes
only the first-apex/last-apex edge. The original components remain internally
strong; forward cross connections through old vertices remain available;
the reversed apex edge closes the return route. The strong conclusion is
therefore genuine. For the final one-apex component, the exceptional t=0
value41 is retained, and t>=1 uses3^t(45*4^t-1). Its3-adic valuation is t.
The new-edge difference has valuation at least t+2, so cannot cancel this
leading factor. No r=0 all-strong claim or prime production is inferred.

The independent checker now constructs and counts the full graphs for
(r,t)=(1,0),(1,1),(1,2),(2,0),(2,1), using its own subset endpoint DP. It
reproduces41,537,6471 for the one-apex seeds and the changes1681->7009,
22017->107985 at orders14 and17. It checks the corresponding9 and27
divisibilities, exact remaining3-adic valuations, strong reachability,
and every old source arc after deleting the retained named apices.
Normal and optimized Python outputs agree.
