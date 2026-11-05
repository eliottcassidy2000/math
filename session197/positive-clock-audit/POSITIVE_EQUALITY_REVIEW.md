# Independent review: unique positive-clock equality core

Status: **PASS, no correction required.**

Frozen source files in `session197-positive-clock-audit/`:

- `EQUALITY.md`: SHA256 `96c612b1fa50648dfe3d04fb8da5ab35d8d7945650743d026086b0dc36543ee7`.
- `check_positive_clock.cpp`: `cb04f653b223b4152e2b24ecfe4fe1c352816a9946e1ba56b3b4fe318b473f81`.
- `independent-opt.txt`: `6459b8349c896e476620b70a3464e409bd52a1e4e8d6663afea9d8a5964fda2c`.
- `equality-controls.json`: `be0c5306d3229f52b2a4f0e559e14116633d3a4abe18ad764fc6abcfdf71a28c`.

The SCC gcd reduction uses gcd(217,649)=1 correctly. All outside components are singleton factors. The root-fixed anti-involution preserves the strong root component and reverses the SCC order, forcing equal singleton padding on its two sides. At count217, the exact Busch order-eleven minimum225 is needed; the weaker uniform estimate at order11 would not suffice. Together with the order-at-least12 bound and oddness of the core, this reduces the equality case to the complete order-at-most-nine bank.

The source exact isomorphism test preserves the marked root and checks every adjacency relation. Its separate unmarked automorphism test establishes full rigidity, not only a trivial root stabilizer. Its all-root test correctly pads rectangular crossing matrices with zero rows or columns; it does not assume every candidate root is balanced. Thus the root uniqueness claim is valid.

The converse equal-padding construction preserves every clone count and admits the stated anti-involution. Ordered SCCs force every automorphism to fix each singleton and restrict to a core automorphism. Therefore rigidity and unique-root status persist. One unlabelled source at each odd order n>=9 consequently has n! labelled copies, n! labelled rooted pairs, and (n-1)! copies with a prescribed root label. Other orders are excluded. These counts are not divided by an unproved symmetry factor.

## Independent controls

A separate checker, `check_positive_equality.cpp`, uses the complete Hamilton-word-normalized universe rather than the source's O/I matrix universe. It finds exactly eleven amplitude72 presentations and maps each to the displayed rooted source by an independently implemented exact backtracking isomorphism test. It independently checks the full automorphism group is trivial and tests every root using subset matching and physical clone counts. Only root4 has the required rank and counts.

The same checker enumerates all 9! vertex permutations, tests Hamilton adjacency and reversal-anti invariance directly, and finds exactly eleven fixed words. This differs from the source's half-word enumeration. Every one of the eleven explicit certificate words was additionally validated for vertex set, adjacency and symmetry. Hence the normalized mass and physical-presentation multiplicities are consistent: eleven fixed words and 4!=24 O-side labelings.

Optimized and unoptimized library-checking independent runs agree. The frozen source checker was replayed and reproduces its receipt byte for byte.

Independent checker SHA256: `2a46e0a6e5ac18d974ba4d13c9914dc40e73c4c1552e0f6bf8626c05fc18b967`.
Independent receipt `positive-equality-independent.txt` SHA256: `63a71bf1f38369d4947bcfc0f9ee40cb995e304799dac372cbd93de53dc6c982`.

The equality classification and labelled/unlabelled counts pass at the frozen source hashes above.
