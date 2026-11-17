# The asymmetric port gap is a binary address

**Proved:** after any fixed asymmetric prefix and suffix of the triangle frame, the change in the seam gap under further shared contacts is exactly an affine binary address with one bit correction. This decodes every shared suffix which restores a desired gap. It does not decode the absolute midpoint.

Use the two-port triangle of the [original switch](../three_visit_triangle_switch_20261117/README.md), with the triangle at any position of the transitive exterior order. A contact digit has bit 0 if x beats the vertex and bit 1 if a beats it. Digits 1 and 2 contact the ports unequally. Fix any such exterior, retaining both orientations of x–a. Subsequently append only shared suffix vertices, either S (beaten by both ports) or R (beats both ports); all old nonports beat the new vertex.

Only paths using x–a can change. Contract the adjacent port pair, leaving two increasing external runs around it. Maintain the six states (L,R): L is empty, endpoint incompatible, or endpoint compatible with the incoming arc to the first port; R records only whether the run after the second port is occupied. For each orientation separately, the state counts are computed by routing every exterior vertex into either run, checking entry at the outset of the right run, and inserting the triangle by ordered cuts of its three vertices. At the finish, read the total for states with L empty or compatible. Subtract the minus state vector from the plus state vector.

Write D for this signed finish count, U for the signed mass of states with R occupied, and B for the signed mass of the state with L nonempty incompatible and R empty. Crucially the signed total of all R-empty states is zero: these paths simply put the whole nonport tournament on the left, not yet checking its junction to the first port. Thus its signed mass cancels in the two orientations.

* Appending S preserves D: its left insertion ends incompatibly and its right insertion freely extends every acceptable old state. Appending R adds U+B to D: every state now admits a compatible left insertion, and preexisting occupied right runs may continue.
* Either letter doubles U. For R, only occupied right runs admit the right insertion, and they also admit a left insertion. For S, all states admit the right insertion, but their total signed mass equals U, since the R-empty mass is zero. In either case B becomes zero: only the old R-empty states can contribute to B, always together or not at all.

Recording the subsequent word in chronological order with bit j set for R at step j (starting with j=0), put z=Σ[R_j]2^j. The **signed gap** is therefore

    D_final = D + U z + B (z mod 2).

For U≠0, each prescribed signed gap has at most two suffix addresses, ignoring neutral trailing S's: solve the integer affine equation once for each starting bit, then check divisibility, sign and parity. For U=0 a single S preserves D and makes both future corrections vanish; already stable for all words iff B is also zero. Thus the three probes empty, R, RR alone decide stability under all shared suffixes. They do not test asymmetric letters in the continuation.

## Checks and an exact return

| Fixed prefix | Fixed suffix | Initial signed gap D | B | U | exact unit addresses z |
|---|---|---:|---:|---:|---|
| 21 | empty | −2 | 0 | −12 | 0 |
| 1 | 132 | +2 | 3 | 105 | 0 |
| empty | 21 | +10 | −3 | −3 | 3,4 |

In the last line the fixed pair has counts (103,113). The arithmetic decoder, taking signed target −2, gives addresses 3=RR and 4=SSR. Literal full counts are (1307,1305) and (1379,1377), respectively. The signed cancellation is paid by the address, not by a permanent endpoint balance. This is a different exact choice than choosing a tile with an identically preserved gap.

Scanning all exterior prefix/suffix splits through seven vertices with at least one unequal contact found respectively 1,2,18,62,206,558 unit-gap addresses at lengths 2,...,7, but none was stable under all shared suffixes. This is a finite statement only. Run model.py for the state test and the scan. The useful design equation for future cores is now U=0 together with the required initial D; otherwise a desired correction consumes a uniquely constrained suffix choice per starting bit, which can then be tested for its absolute midpoint.
