# Independent extension audit: complete d=1 and d=3 membership

Status: PROVED independently as below; producer extension draft not yet frozen when this file was first written. This theorem classifies low deletion-count sidecars, not unrestricted Hamilton-path counts.

## Deletion structure

Hamilton-path count is nondecreasing under adjoining vertices: any path admits an insertion of the new vertex (use the first vertex it beats, unless there is none, in which case append). Fix one insertion deterministically; deleting the new vertex recovers the original path, proving injectivity. Thus an induced subgraph's count is a lower bound.

A tournament has H=1 iff it is transitive: a directed triangle has H=3, and every nontransitive tournament contains one. A transitive tournament's unique vertex order gives its unique Hamilton path.

If H=3, choose a directed triangle C. Any outside vertex with mixed adjacency to C produces the strong four-vertex tournament, with H=5, impossible. All outside vertices are therefore in L (beating all of C) or R (losing to all of C). Neither L nor R can contain a directed triangle: it and C form an ordered pair of C3 blocks with H=9. A backward edge R→L would induce the five-vertex cyclic composition of L-singleton, C3 and R-singleton, with H=15, also impossible. Thus L and R are transitive and L→C→R. Conversely every such ordered sum has H=3. This is an elementary structural proof, without relying on an unaudited all-order strong-floor theorem. The classification itself is inherited low-spectrum knowledge, not a new priority claim.

## d=1 rigidity

For root-deletion TT_m, its unique ordering and root neighborhood yield a binary word b. The one-block U recurrence gives U=1+binary_value(b). Leading zero vertices are global sources of the entire rooted graph; deleting them preserves D and E as well as U and d. Hence for every integer U≥1, the canonical word bin(U−1) determines exactly one (D,E) pair for d=1. The root-alone case U=1 also has (D,E)=(1,1), agreeing with canonical word0. D is computed by the six-state transfer with singleton matrices only; E=1−last_bit, including canonical word0. Not every arbitrarily specified E is possible.

## Complete d=3 decision procedure

Given candidate integer sidecar (E,U,3), reject unless U is a positive multiple of3 and E∈{0,1,2,3}. Put M=U/3−1≥0 and w=00+bin(M), with bin(0)=0. For each adjacent position p in w, compute D by the six-state transfer and E_p by

    E_p = 3−(2w_p+w_(p+1)) if the triangle is the final block;
    E_p = 3(1−w_last) otherwise.

Then the exact set of realizable D at this sidecar is the set of those computed D for which E_p=E. Every candidate is literally realized by the constructed rooted tournament, so sufficiency is immediate.

For necessity, the deletion classification gives TT_l→C3→TT_s. Compress the root's triangle neighborhood size k∈{0,1,2,3} to its two-bit binary word. The universal U formula forces the resulting effective word to represent M. It therefore differs from bin(M) only by leading zeros. If the triangle is not wholly contained in this leading zero prefix, extra leading zeros are global source vertices and can be removed or prepended without changing the profile; the two zeros of w retain all placements crossing the first significant bit. If the triangle lies wholly within the leading zeros, the entire leading prefix is an ordinal sum of singleton sources and one C3 source component. Its scalar multiplier is3 regardless of its position; all such placements have the same profile as p=0 in w. This proves exhaustive stabilization. The M=0 case gives only (D,E,U,d)=(3,3,3,3), consistent with all-zero words and either terminal-block case.

This provides at most len(bin(M))+1 distinct D values for fixed (E,U,3), not infinitely many values at a single numerical sidecar. The mobile alternating family proves the cardinalities are unbounded as U varies. Both conclusions are compatible.

## Independent finite receipts

`classify.py` imports only this audit's `verify.py`, never producer code. Commands:

    python3 classify.py > classify.normal.json
    python3 -O classify.py > classify.optimized.json
    cmp classify.normal.json classify.optimized.json

All checks PASS with byte-identical output. Of all 1,099 labelled root-at-zero tournaments of orders1–5, 443 have d1 and272 have d3; every one satisfies its stated canonical classification. For every M0–127 and every leading-zero count2–7, all7,308 placements yield exactly the same stabilized set of (D,E). These are finite sanity controls; the all-order necessity rests on the proof above.

## Frozen producer addendum verdict

The subsequently frozen `inputs/MEMBERSHIP.md`, sections1–4, is **PASS with no correction required**. Its smaller backward-edge obstruction using just two C3 vertices is indeed the strong four-tournament with H5; the H15 witness above is a different, compatible five-vertex obstruction. The explicit case split around the first1 correctly proves finite-bank completeness. Its d1 endpoint identity E=U mod2 is equivalent to the final-bit expression above. The bitlength upper bound is asymptotically sharp along U=4^r+2 because the lower and upper bounds differ by only2. The statement that3 is the first positive deletion level with unbounded scalar fibers uses the standard oddness of tournament Hamilton-path counts, together with d1 rigidity. No claim about the smallest-order collision is implied.
