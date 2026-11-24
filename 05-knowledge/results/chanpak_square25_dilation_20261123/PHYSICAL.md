# A physical square seam can skip a scale only by reversing its sign

Write \(r_m=2^m+1\). A **marked seam** is a pair of tournaments on the same vertices that differ at a distinguished arc. Its readings are their Hamiltonian path counts. In an induced inclusion, deleting the same new vertices from both completions recovers the old pair and preserves the distinguished unordered vertex pair. The sign is preserved when the lower reading restricts to the lower reading. This keeps exactly the placement data that separate numerical coverage from a recursion on actual tournaments.

**FINITE-EXACT.** For \(2\le m\le9\) there is a strong marked seam with readings \(r_m^2-2,r_m^2\). The minimum common orders are as follows. The supplied outgoing masks and distinguished arcs independently generate all sixteen tournaments. The lower bounds follow because the lower reading exceeds the known maximum at the preceding order in every case. Singleton padding preserves height. [Physical witnesses and subset recurrence](physical.py) · [Exact order maxima](../path_spectrum_two_angles_metagraph_20261106/README.md)

| Root | Counts of the two completions | Minimum vertices |
|---:|---:|---:|
| 5 | 23, 25 | 6 |
| 9 | 79, 81 | 7 |
| 17 | 287, 289 | 8 |
| 33 | 1,087, 1,089 | 9 |
| 65 | 4,223, 4,225 | 10 |
| 129 | 16,639, 16,641 | 11 |
| 257 | 66,047, 66,049 | 11 |
| 513 | 263,167, 263,169 | 12 |

## A classification over all orders for the first three seams

**PROVED by finite exact exhaustion.** Consider the three pairs of readings \(A=(23,25), B=(79,81), C=(287,289)\). Allow arbitrary representatives at arbitrary orders and either sign for the inclusion.

- No marked \(A\) occurs as an induced pair inside any marked \(B\). No marked \(B\) occurs inside any marked \(C\).
- Marked \(A\) does occur inside marked \(C\), always with reversed sign. Up to simultaneous converse there is exactly one inclusion of their strong active cores, with respective orders six and nine.

For an explicit inclusion, use the nine outgoing masks

```text
322 488 451 341 487 461 256 73 128
```

Reversing the arc on vertices 4 and 5 raises the full count from 287 to 289. Restricted to the first six vertices, the same reversal lowers 25 to 23. Both full completions are strong. Among the new vertices, 6 and 8 form a transitive suffix below all six old vertices, and the last inserted vertex 7 beats old vertices 0 and 3. The arc on 4 and 5 has no changed endpoint. Thus a two-step jump preserves its physical location and exchanges which orientation gives the larger reading. [Exact replay](physical.py)

### Why the search answers the arbitrary-order question

In a unit seam \((h,h+2)\) with \(h\ge5\), replacing the marked arc by a digon identifies a common active core. Both orientations of that core are strong. The common exterior has height dividing \(\gcd(h,h+2)=1\), hence is a transitive prefix and suffix. An old strong core embedded as an induced subtournament remains inside the strong component of the new tournament containing its marked vertices. Consequently any induced inclusion at any orders restricts to an induced inclusion of the two active cores. [Strong active-core theorem](../chanpak_arc_fiber_grammar_20261122/README.md)

The Busch lower bound for strong tournaments bounds the core orders by 6, 9, 11 at the lower readings 23, 79, 287. Enumerate the possible source cores by deleting one endpoint of the arc, taking one representative for each possible parent, then restoring the endpoint with each contact mask. The two marked readings for an apex with mask \(u\) and pivot \(i\) are the entries indexed by \(u\) and \(u\mathbin{\mathrm{xor}}2^i\) of the same all-apex table. Canonicalize the strong marked results up to simultaneous converse. For \(A\) this gives two classes at order six. For \(B\) it gives two each at orders seven and eight and none at nine. The only representatives omitted by the converse identification have the same extension possibilities after reversing every arc. [Busch bounds and component grammar](../chanpak_arc_fiber_grammar_20261122/README.md) · [Parent recurrence](../chanpak_arc_fiber_grammar_20261122/BIRTH-CUT-COVER.md)

Starting from each core, add **every** common new-vertex contact mask in every remaining position up to the target cap. If one of the two running counts exceeds its prescribed terminal count, discard that branch. This uses monotonicity under arbitrary vertex insertion. From each old path, insert the new vertex at the beginning if it beats the first vertex, at the end if the last vertex beats it, or between successive vertices at a transition from an inneighbor to an outneighbor of the new vertex. Choosing one such insertion per path is injective by deletion. All intermediate tournaments, including reducible ones, are admitted by the search. Each marked target containing the source has some ordering of its additional vertices and hence is reached. [Exact extension program](bank_scan.cpp)

| Source → target | Expanded parent records, preserving sign | Expanded parent records, reversing sign | Exact hits |
|---|---|---|---|
| \(A\to B\), parent orders 6–8 | 2, 50, 218 | 2, 49, 208 | 0 in either sign |
| \(B\to C\), parent orders 7–10 | 2, 88, 437, 2,080 | 2, 86, 430, 2,044 | 0 in either sign |
| \(A\to C\), parent orders 6–10 | 2, 128, 4,884, 37,298, 245,512 | 2, 128, 4,854, 37,022, 243,352 | Preserving 0. Reversing 6, 48, 360 at orders 9, 10, 11 |

To obtain the uniqueness statement, fix the labelled six-vertex source in each of its two representatives modulo converse. The automorphism group of a marked unit pair acts freely on the Hamilton paths of each completion. Its order divides their coprime counts and the group is trivial. Therefore every target with a specified image of the fixed source core contributes exactly \((n-6)!\) labelled extension histories. At order nine the six histories give one inclusion. At orders ten and eleven, distribute the additional vertices between a prefix and a suffix of that same strong target. These give two and three distinct inclusions and account for all \(48=2\cdot4!\) and \(360=3\cdot5!\) histories. A strong target core of larger order violates the Busch bound. The statement concerns the inclusions of active cores. Extra source singleton positions, when present in an embedding, are forgotten by this reduction. [Marked automorphism lemma](../chanpak_arc_fiber_grammar_20261122/README.md)

## Consequence for a recursive construction

Replacing any collection of vertices away from the arc by arbitrary nonempty tournaments always leaves an induced copy of the base pair. Choose one representative of each replacement block in both completions. This applies to unequal flanks, simultaneous enlargement of any number of ports, and arbitrarily many iterations with fixed arc endpoints. The two consecutive non-inclusions therefore rule out every such recursion starting with \(A\) or \(B\) and moving one dyadic scale at a time. They also rule out arbitrary ordinary vertex extensions retaining that seam. The reversed inclusion gives an explicit construction past this obstruction when one scale may be skipped. A uniform walk across further roots is open. [Related independent lower and upper constructions](SQUARE-RAILS.md)

Run `python3 physical.py` in this directory for the eight minimum-order seams and the reversed inclusion. `python3 physical.py --census` exhausts both adjacent stages. `python3 physical.py --skip-census` also exhausts the classification of the two-scale jump. The optional censuses use nauty's `gentourng` and a C++17 compiler.
