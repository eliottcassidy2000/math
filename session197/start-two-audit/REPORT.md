# Independent audit: start-one intervals and the start-two grammar

**PASS**, no correction. This report is bound to the exact three hashes in SOURCE_SHA256SUMS. Frozen sources were copied for replay and left unchanged.

## Proof review

The exact old-root seam equation is valid: paths avoiding the seam contribute g*s_x−t*p_(x,a), and paths using it split into the ordered two A-runs starting at x and a and the gadget cover whose first run ends at b. The empty final gadget run is included in L. For the marked old root itself, seam use is impossible, giving g for nontrivial A and g−t for singleton A.

With a nonempty transitive tail, t is0 or1 and L>=2. If t=1, the root and the first/final tail vertices contain a triangle, so g>=3. Thus no other old vertex retains start count one. The tail formula counts the increasing initial tail fragment up to its final vertex that beats b, then the unique marked A-path and the complementary increasing tail. It gives precisely the claimed two/three-vertex consecutive block and outdegree pattern, including the all-zero triangle exception.

With an empty tail, D>=1 comes from the two nonempty segments of the old unique cycle. Equality s_x−p_(x,a)+D=1 forces s_x=p_(x,a). Every such path closes through a→x, so uniqueness of the cycle forces x to be the cycle successor of a and s_x=1. The inductive interval then identifies this case with outdegree(a)=1. Conversely, in a two-cover starting at x and a, that latter path must be singleton a; the remaining x-start path is unique. Hence D=1. The singleton-gadget step and all small cases are sound.

The refined Fibonacci identities use rigidity only at orders at least4. The triangle's unordered start-one pairs form one orbit, not three. This gives D2(z)=z^3(1−z^2)/(1−3z+z^2).

Deleting a start-two root has the exact SCC start-count product stated in the proof. Any later nontrivial component would contribute a factor at least3; hence the tail is transitive and the initial entry sum2 is either one start-two vertex or two start-one vertices. The entry set cannot be the entire initial component, ensuring an incoming edge when the tail is empty; a nonempty tail's last vertex supplies the return edge otherwise. The root-deletion SCC decomposition and marked unordered entry set are intrinsic, so no construction-history overcount occurs.

The resulting generating function is z^4(1−z)(1−z^2)/(1−3z+z^2)^2. Its numerator ends in degree7, so the displayed fourth-order homogeneous recurrence begins at n=8, as claimed. The class is rooted start count two, not all C=2 graphs; the regular five-tournament is correctly excluded.

## Independent checks

Both source normal and Python -O runs reproduce the frozen certificate byte for byte.

The separate independent.py uses a starting-endpoint subset recurrence, explicit permutations for cycles and paths, and an independent transitive-closure SCC decomposition. It checks all33,867 labelled tournaments through order6. All6,146 unique-cycle graphs satisfy the interval/outdegree law. All13,584 labelled root instances with start count two have exactly the claimed intrinsic SCC grammar. Their rooted class counts are1,5,18 at orders4,5,6. A full factorial root-fixing canonicalization of representatives independently confirms the source's two-Hamilton-path canonical keys do not identify distinct rooted classes. See INDEPENDENT.json.

The separate seam_checks.py enumerates actual one-seam gluings from all start-one sources through order4 and every permitted transitive-tail gadget through four tail vertices. Direct subset counts verify1,263 physical gluings,3,648 old-root seam rows and3,871 tail start-count rows, including the empty-tail and singleton-source boundaries. See SEAM_CHECKS.json.

These bounded checks validate the implementations and controls. The all-order classification and generating functions rest on the induction and canonical SCC arguments reviewed above. No unrestricted spectrum completeness or complete residual C=2 classification is asserted.
