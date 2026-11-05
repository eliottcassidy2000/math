# Independent audit: exact replacement fibers and the H=43 boundary

Verdict: PASS on the revised title, which makes no unproved minimality claim.

- PROOF.md SHA256: 8045be10920a9b5a3cdc629586de010cb322bd989b8279469ed9607b39c9beee
- verify.py SHA256: 663cf4c9313ddd6e6605987a191fd9632375f41a517bc8ac773f268872007af3
- positive_model.py SHA256: f7790919712d11692a3ba1f463f272b2193603af61f879148b4f54b3b503c146

The positive model is byte-identical to the previously reviewed finite-grammar
verifier. All four Python dependency files are also unchanged. Isolated normal
and optimized runs pass and reproduce the frozen receipts exactly.

The two replacement rows are exact linear identities on source profiles. The
three infinite source families are physical words in O/Z/J. Their prefix
formulas follow by induction under Z or O, and the displayed suffixes give the
claimed full profiles. The linear response corrections are respectively -1,
+1 and +1. The stated coprime differences prove primitivity for the first two
families; the unit coordinates prove it for the third. Parent/child values and
orders follow as stated. These are genuine specified-fiber constructions,
without an assertion that every desired parent meets such a fiber.

The complete C2 inverse bound is 7 at target43. Its seven source profiles yield
exactly the stated eleven product-admissible tests and the single unordered
profile pair N(t,O(triangle)). The conservative C3 bound14 includes every
permitted singleton boundary. Its twenty profiles and 2,964 admissible tests
have support {15,27,29,31,33} through43, so no C3 presentation supplies a rebase.
These are justified all-order inverse exclusions, not depth-limited searches.

The entrywise comparisons P_F>=P_N and 3 G^t P_N>=4 P_N are correct. Since
physical nontrivial sources have nonnegative positive coordinates, three total
source letters force every final gate above87. All 34 pairs of words of total
length at most two and ten gates give precisely {43,69,77,79,81,83} through87.
Transpose closure handles the reversed ordered parent. Compiling indexed
source moves into O/Z/J shows that no inherited physical source move is omitted.

The graph-class step is also valid. A source with H1 is transitive, and U2
forces the two-vertex source. A start-one H3 source has a triangle initial
component and singleton tail; any strong component with H3 has order at most3
by cutting a Hamilton cycle. Its U6 thus forces exactly one tail vertex.
Rooted source choices are unique up to isomorphism. The complete C2 graph normal
form therefore leaves exactly the two ordered parent constructions, which are
converses. The exhaustive 128 births rule out85 and87 for either class,
regardless of the child's cycle count. This conclusion is only about adding
one vertex to these C2/C3 parent classes, not arbitrary-depth embeddings or
parents with more cycles.

Additional independent controls use no code imported from the packet:
- indexed M_L/O source enumeration, instead of its positive O/Z/J bank;
- all raw-coordinate C2/C3 inverse and 340 descendant evaluations;
- full Hamilton-word counts of all five parent/child graphs;
- enumeration of all eight-vertex Hamilton words, classifying each word by
  its required birth-mask bits, which reproduces all128 individual birth counts;
- exhaustive labelled rooted tournaments through order4 for the two parent
  source profiles.

All pass. See ../target43-independent.py and ../TARGET43_INDEPENDENT.json.
The explicit eight-vertex outputs85(C2) and87(C2 and C3) are verified. They
escape the stated extension architecture by replacing the parent structure.
No general odd coverage or first-failure theorem follows or is claimed.
