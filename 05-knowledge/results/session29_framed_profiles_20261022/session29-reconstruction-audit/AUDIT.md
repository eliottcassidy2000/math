# Independent audit: intrinsic seams and rooted target reconstruction

Status: **ACCEPTED**, subject only to the source hashes recorded alongside this audit. No mathematical correction was required.

## Scope and proof review

This audit accepts the all-length classification of the strong targets reached from the specific binary A/B join family, not classification of the entire tournament metagraph. Its key mechanism is the intrinsic one-arc-cut chain, not extrapolation from a finite order-10 census.

For incomparable sets A,B with one outgoing arc each, put P=A intersection B, U=A union B, Q=A minus B, R=B minus A. The exact cut identity is

 d+(A)+d+(B)=d+(P)+d+(U)+|Q||R|.

When P is nonempty and U is proper, strong connectivity makes the right side at least three. If P is empty, |A||B|>=4 exhausts the two budgets; if U is the full vertex set, |Q||R|>=4 does so. These cover all incomparable cases under the stated size restrictions. Thus the relevant cuts form a chain and each prescribed size5j has a unique cut.

The seam construction is strong, since its fixed Hamiltonian path closes using the reversed first-last chord. Consequently every isomorphism preserves each cut of size5j, hence every five-vertex block. Unique outgoing seam arcs identify the required roots. Conversely marked-block isomorphisms assemble into a target isomorphism. This proves both necessity and sufficiency, and the direct-product automorphism statement. Converse takes a prefix cut to the complementary suffix cut, which proves the address-reversal rule rather than merely suggesting it.

Every finite rooted letter class is nonempty, so every Cartesian-product address is realizable. Distinct addresses are not silently identified. The counts169*17^(r-2), the self-converse counts, and the merged Burnside GF follow. The ordinary source words remain distinguishable by their ordered strong-component types A and B. The line-incidence product counts all actual tile choices, and there is no factor-of-two quotient on lines: source vertices are nonstrong while target vertices are strong. Converse merging changes endpoint assignments, not the tile/line population.

The odd-length central class {755,949} is genuinely fixed as a doubly rooted isomorphism class but has no fixed tile representative in this alphabet. The resulting no-blue assertion is correctly restricted to incoming lines from this family. It does not say those self-converse targets lack blue lines globally.

## Independent implementation

`independent.py` imports no producer Python module. It rebuilds the18-letter alphabet from Hamiltonian permutations of the two five-vertex seed matrices, then independently performs the24-permutation rooted and6-permutation doubly rooted canonicalizations. It reconstructs every two-block target directly from its arcs.

Verified:

-13 outer rooted classes,17 internal doubly rooted classes,5 fixed internal classes;
- every one of324 full10-vertex subset-DP Hamiltonian counts;
- every supplied isomorphism/converse vertex witness, arc by arc;
- every balanced one-arc cut and nesting of all relevant cuts in these targets;
- all source/target/color incidence products and169 ordinary/91 merged target classes;
- marked-block automorphism factors, which give tile-weighted frequencies289,34,1 at automorphism orders1,3,9;
- connected94-vertex restricted multigraph, exactly8 bridges, all blue;
- all316 supplied alternative paths, hence genuine global nonbridge certificates;
- connected86-vertex black boundary component with306 edges.

For the distinct-class Hamiltonian generating functions, this audit uses the already independently certified three-state letter matrices, verifies that their columns/matrices descend to the newly reconstructed rooted classes, and checks the ordinary rational resolvent directly. The fixed-target series are checked with a separate **nine-coordinate Kronecker transfer**, not the producer's six-monomial implementation. Besides18 coefficient checks, the claimed denominator annihilates the nine-state observation row exactly. Therefore the rational identities are proved at all orders, not inferred from the matching prefix.

The verified ordinary target mass GF is

 (93313z^2+3359832z^3+4057812z^4)
 /(1-823z-37140z^2-45180z^3).

At r=2, ordinary mass93313 and fixed mass7661 give merged mass50487. These are distinct-class sums, not the previous edge-weighted mass195472.

## Universal multigraph bridge statements

The black-pair proof is valid for n>=3. F and R commute. Neither F nor FR has a fixed tile, because the first-last chord is a coordinate fixed by R and toggled by F. A black F-edge therefore has a distinct mirror image with the same converse-merged endpoints. Hence a raw black edge is never a bridge. The assertion is not valid after suppressing parallel edges, and the report explicitly avoids that conflation.

The parity-charge reduction is also valid: automorphism groups of tournaments have odd order, since an involution would swap the endpoints of some arc. Hamiltonian counts are odd, and the automorphism action on ordered Hamiltonian paths is free. A self-converse merged node therefore has odd degree H/Aut, whereas a nonself-converse merged node has even degree2H/Aut. All black boundary counts are even. Any bridge must be blue and separates an odd number of self-converse nodes on either shore. Contracting full black components preserves the blue-edge bridge question; the quotient degree parity equals the parity of the number of self-converse nodes in that component. Loops count twice in degree. None of these full-metagraph degree statements may be applied to truncated target degrees in the restricted boundary graph.

## Limits retained

The eight restricted bridges are not certified global bridges. The other316 are certified global nonbridges because their alternative paths are already present. No claim about the full simple graph's bridges, all strong tournaments, A000568 asymptotics, or the global Hamiltonian-count gap is inferred.

Normal and optimized-mode independent outputs are byte-identical. The source and deliverable manifests specify the audited snapshot.

## Expanded fixed-k theorem and polynomial recognition

The final source sections9-11 were included in this audit. The reconstruction argument extends verbatim to a mirror-closed alphabet of fixed-size k>=2 tiles with unconstrained independent word choices. Source blocks need not be strong: their strong components have a canonical total order, no component crosses a join boundary, and cumulative component sizes recover each prescribed size-k boundary. This establishes source reconstruction without silently requiring strong alphabet letters. It does not extend the special three-state Hamiltonian observer to arbitrary alphabets.

The O(n^3) cut recognizer is correct. A strong bridge must occur in the union of a root's in- and out-arborescences; otherwise their surviving union still gives strong connectivity through the root. Deleting one arc from a tournament leaves SCC reachability totally ordered: every pair of SCCs has an arc unless both are the singleton deleted endpoints, and those endpoints still have a reverse path. All zero-outgoing cuts are therefore suffix unions. Restoring the deleted bridge converts exactly these into one-outgoing-arc cuts. At most2n-2 candidates, each requiring O(n^2) adjacency-matrix work, suffice.

`check_recognition.py` independently verifies this mechanism without importing producer code. It deletes every arc, not just candidate arcs, in all570 strong labeled tournaments of orders3,4,5, uses Boolean transitive closure rather than the producer's SCC implementation, and compares the recovered suffix closures to exhaustive subset cuts. All5590 arc deletions and1438 strong-bridge instances pass, including membership in the independently computed arborescence union. Its normal and optimized outputs are identical.
