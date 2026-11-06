# A one-arc edit repairs each induced-preservation hostile

Status: **VERIFIED exact finite witnesses**, independently audited. This is not an all-target constructor.

The all-order induced-containment theorem in `INDUCED_EXTRACTION_BARRIER.md` forces any construction of the five indicated target counts to alter the old induced E(u) somewhere. The obstruction is sharp in the following precise edit sense:

**Among operations that keep all old vertices, reverse old arcs, and add arbitrary new vertices, the minimum number of old-arc reversals is exactly one for each of the five hostile pairs.**

The lower bound is the all-order theorem, not a two-birth-only search. For the upper bound, reverse the marked arc r->c in every case, then use the following outgoing masks for the first and second new vertices:

| u | Target | First mask | Second mask |
|---|---:|---:|---:|
| 9 | 55 | 0 | 82 |
| 17 | 103 | 97 | 229 |
| 19 | 113 | 49 | 113 |
| 35 | 209 | 0 | 335 |
| 35 | 211 | 0 | 409 |

Old names are r=0,c=1,p_j=j+2. New names are appended in order. A new vertex's mask records exactly the old names it beats; all remaining edges point into it. `ONE_ARC_REPAIRS.json` retains every literal row, order and verified count. The five recorded graphs are nonstrong; no strongness claim is made or needed for positive existence.

This operation performs an actual physical reduction before adding vertices: reversing r->c makes r a global sink, and the count drops from u+g(u)-(u mod2) to g(u). It does not treat subtraction of the unwanted sector as a free formal operation. The later incidence masks are genuine target-dependent data, not a common recipe proved to work for every u.

## A separate repair of 391

A separate finite target considered was 391=6*65+1. It is attained with a different old-arc reversal:

1. In E(65), reverse c->p_0. Since the binary word of 65-1 is 1000000, this makes the whole old nine-vertex tournament transitive.
2. Add a universal sink, outgoing mask 0.
3. Add a vertex with outgoing mask 57 on the ten existing names.

The resulting eleven-vertex tournament has outgoing rows

    (514,1536,2043,1011,995,963,1923,1795,1539,1024,57)

and exactly 391 Hamilton paths.

This is an exact finite repair, not a power-of-two formula. Applying the identical named-arc and mask rule to E(129) gives 783 rather than the desired 775. That hostile is explicitly retained in the verifier.

## Reproduction

    python3 verify_one_arc_repairs.py
    python3 -O verify_one_arc_repairs.py

Both produce identical output. The verifier builds the graphs from the physical operations, verifies every tournament edge, compares literal rows, and independently counts Hamilton paths by full vertex-subset DP. It also checks the failed adjacent power-of-two generalization.

The exploratory scan is not included in this release packet. The proof of the five one-arc upper bounds and the 391 repair needs only the explicit witnesses and the included verifier, not completeness of that scan.
