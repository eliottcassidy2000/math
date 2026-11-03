# Independent audit of S173 square-root words

**PASS after one corrected boundary condition.** Audited final NOTE.md, check.py, twins.json and final certificate.json at `/tmp/s173_squares`.

## Proof review

- Theorem 1: alternating affine coordinates give exactly the stated bounds; duplicate offsets within a parity force collision and cross-parity collisions give precisely 2x=o−e. Empty root words cause no exception.
- Theorem 2: the first moment fixes the start for odd n; for even n the second moment has positive quadratic coefficient n and gives at most two starts. This does not assert existence.
- Theorem 3: equality of edge sums forces alternate ±t shifts. Each residue-t chain must permute by adjacent swaps, yielding even chain lengths and exactly the prescribed alternating t-block coloring. Sufficiency is the reverse construction. Both n=36 paths, every vertex/edge, and their exact ±6 differences were independently checked.
- Connector cardinality: deletion of k extra vertices yields at most k+1 old pieces. The three/four old components therefore need two/three additions, attained by the explicit paths.
- Added height result: at equality in the component bound, there are exactly c old pieces, exactly one per old component; added vertices are internal and isolated along the path and connect endpoints of distinct whole components. The component 527–2–7 therefore needs a seam through 527 or 7. Exhausting possible one-vertex seams to the other component endpoints gives minimum 317. The two displayed paths attain this, proving the conditional height optimum.

The initially stated general one-vertex-connector iff omitted exclusion of endpoints. For a=2,b=7, the sole positive candidate is z=2, so it repeats an endpoint. The author corrected the theorem to exclude endpoints and occupied vertices and added a regression. The named connectors and main theorems were unaffected.

## Independent controls

`independent.py` imports no author code. Normal and optimized outputs agree. It verifies:

- 4,095 interval controls with arbitrary, possibly negative edge sums and translated/negative intervals;
- 22,748 alternating-shift partitions, including odd orders absent from the author's even-order controls;
- both explicit square-sum twin placements and retained Hamilton paths on 32,44,72,527;
- anchor components using union-find, independently of the author's DFS;
- all height candidates by direct integer z scanning, independently of factor enumeration. The finite bound z<|a−b|² follows from |a−b|=(v−u)(v+u)>2u, u²=z+min(a,b). These rows agree exactly with the author's table;
- the endpoint-repeat hostile example.

The final author checker was rerun under normal and optimized Python. The two reruns match each other and the final supplied certificate byte-for-byte. Its finite run establishes paths through 527, not an all-n termination theorem.

## Source verification

The Jan 21, 2018 contribution signed by Robert Gerbicz in [OEIS A090461](https://oeis.org/A090461) directly states the all-n path range n≥25 and cycle range n≥32. The packet correctly treats this as cited prior work rather than its own proof. The underlying forum proof was not independently inspected here.

The [primary repository README](https://github.com/devisnik/square-paths) contains the exact 32-cycle used by the checker and the cycle cutting/reinsertion description. The [GerbiczExtender source](https://github.com/devisnik/square-paths/blob/master/src/main/java/extend/GerbiczExtender.kt) and [CutAndInsertExtender source](https://github.com/devisnik/square-paths/blob/master/src/main/java/extend/CutAndInsertExtender.kt) support the packet's brief algorithm descriptions.

No forbidden canon file was read. This is a targeted proof, arithmetic, implementation, and source audit, not a literature novelty audit.
