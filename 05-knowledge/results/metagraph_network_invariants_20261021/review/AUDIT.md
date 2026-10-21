# Independent audit: reflection-paired metagraph networks

Verdict: PASS for the explicitly scoped statements and finite-exact controls in the final producer PROOF.md. This is an independent mathematical and computational audit, not a claim of literature novelty, all-order connectivity, or metagraph classification.

Final audited proof SHA-256:

    0528b1db1005b29f1db631cda03c8b7b63a5a794b08b60e84c32784c402d136d

The parent requested one terminology-only correction before freezing: the order-five coefficient 20 and its term are now called “lowest-degree,” not “leading.” The original and final hashes and exact scope are recorded in proof_wording_revision.json. The producer's SHA256SUMS was refreshed for that proof alone. No repository or Git operation was performed, and no other producer file was changed.

## 1. Mathematical proof audit

### Carrier and multiplicities: PASS

The C-orbits are literal undirected line instances. Reflection R commutes with C and preserves converse-merged owners. Its only possible fixed action on a nonblue C-line would satisfy Rt=Ct, impossible at a fixed free coordinate. Thus each nonblue line belongs to a free two-element R-orbit with the same endpoint multiset. The integer electrical network is exactly G=B+2Q, not B+Q, and not the simple support graph.

Loops count twice in stub degree and mass, once as literal edges, and zero in the electrical Laplacian and every spanning tree. Paired black loops remain distinct in the literal graph. The scripts and independent reconstruction use these conventions consistently.

An individual literal black edge is not a bridge because its distinct reflected parallel mate remains. This says nothing analogous about bridges of Q. The older black-quotient README already proves the constructive bridge families; this packet correctly treats those as inheritance and does not reclaim them.

### Componentwise colour-tree polynomial: PASS

For a fixed connected component of G, B+Q has the same connected support. Its reduced Laplacian L is positive definite. With S=L^(-1/2)L_Q,red L^(-1/2), both S and I-S are positive semidefinite. The zero eigenvalue multiplicity is c_Q-1, the unit eigenvalue multiplicity is c_B-1, and the remaining v+1-c_B-c_Q eigenvalues lie strictly between zero and one. This proves the exact minimum/maximum powers, positive coefficients throughout the interval, and strictly negative nonzero roots. The one-vertex case is the empty determinant 1 and is consistent.

Cauchy–Binet with literal edge columns establishes the weighted matrix-tree identity; loops have zero columns, and parallel edges are distinct choices. The endpoint contraction formulas correctly retain all nonloop multiplicities. Substitution z=2, not z=1, restores the literal black choices. Newton's normalized log-concavity is applied to the shifted degree-r polynomial, with no extra missing zero roots.

No global connectivity premise is used. The finite connected observations at n=3,...,7 are correctly separated from this all-order componentwise argument.

### Critical group and bicycle space: PASS

Modulo two, the integer Laplacian of G equals that of B. The incidence exact sequence

    ker(A_B A_B^T) -> ker(A_B) intersect im(A_B^T)

has kernel ker(A_B^T), of dimension c_B. It gives full Laplacian nullity c_B+b_B. Deleting a root reduces nullity by exactly one: normalize the kernel function at the root; the omitted equation follows from zero row sum. Hence the critical-group two-rank is c_B-1+b_B.

Summing chip counts on blue components gives a well-defined surjection to even-sum component parities. Its dimension is c_B-1; the bicycle contribution is additional and is not misidentified as black connectivity. An NS node is a singleton blue component. For any distinct target in the same full component, the source-target divisor has nonzero component parity, hence even order. If a component contains an SC node the guaranteed rank is at least the number of NS nodes; if it contains only NS nodes the guaranteed bound is one less. The proof preserves this distinction.

For a fixed grounded integer Laplacian, the least common denominator of the full rational solution equals the divisor order: integral multiples of the rooted solution are exactly integral firing solutions. This is a claim about the whole potential vector, not one terminal resistance. The order-five counterexample to the latter strengthening is exact.

### Kron reduction and abstract rewire: PASS

The actual order-five source 2 and sink 0 give resistance 6092/5137, internal voltage 11327/10274, and full voltage denominator 10274. The response conductance is 5137/6092 and the relevant harmonic-extension coefficient is 11327/12184. Preserving Dirichlet-to-Neumann response does not preserve the original integer firing lattice without the extension data. These statements are proved by direct rational equations, not by a numerical approximation or a metaphorical quotient.

The four-edge Q switch removes present edges, leaves all multiplicities nonnegative, preserves all paired degrees and loops, and keeps B unchanged. Thus literal degrees, masses, colour counts and charges remain unchanged. The tree counts 41096 and 39232 differ. The packet explicitly calls this an abstract admissible network rewire, NOT an alternative realizable tiling-to-owner map; this limitation is necessary and respected.

## 2. Electrical versus random-walk conventions

The audited conclusions concern the physical integer conductance Laplacian. They do not identify Q with the same electrical network as G, or either graph with its simple support. They also make no holding-time or lumped Markov-dynamics claim.

If a future extension uses the uniform-half-edge node Markov surrogate, let A have off-diagonal literal conductances and diagonal twice the literal loop count, and let D contain full stub masses. Then P=D^(-1)A includes loop holding probabilities, while the electrical Laplacian is D-A and the loop contributions cancel. Discarding the loop ledger preserves direct-current energy but need not preserve that clock. Moreover, this node chain is a re-randomized compression, not automatically the projection of deterministic complement on individual tilings. The cited endpoint-transport result already distinguishes them. No additional dynamical theorem is approved here.

## 3. Independent finite-exact checks

independent_audit.py uses only the Python standard library. It imports no producer code and enforces content pins before using either producer outputs or atlas inputs. Its separate paths are:

- Full-adjacency canonicalization under every vertex permutation through n=5, including a separate converse merge and SC-label check. This does not reuse the producer's Hamiltonian-path-filtered fibre algorithm.
- Reconstruction of all C-lines and actual R-pairs through n=7 from the pinned tiling owner map. It checks edge multiplicities, source masses, loop stubs, colours, reflection degrees and charges.
- Integer Bareiss determinants, exact finite-difference tree-polynomial interpolation, both contraction endpoint formulas, and exact Sturm certificates for strictly negative roots.
- Full literal blue incidence matrices with distinct parallel columns and zero loop columns, independently checking the bicycle and critical-group ranks.
- Fraction Gauss–Jordan inverses and direct full Kirchhoff equations for every ordered NS-source pair through n=6.
- Critical-group exponent recovered from the common denominator of the inverse reduced Laplacian. Together with order and two-rank, it independently fixes the claimed Smith invariant lists in these cases: after the largest factor, the product of remaining nonunits is exactly 2^(rho-1), so all are 2.
- Componentwise theorem checks on 888 small generic component cases, including disconnected whole graphs and isolated vertices; explicit parallel-edge, loop and Q-bridge hostiles. This is not an expanded tournament census.

The independent output records 109301 explicit checks. Normal and python -O outputs are byte-identical. All guards use explicit exceptions, not assert statements.

Actual finite records:

- n=3: 2 nodes, c_B=1, bicycle=0, two-rank=0; Psi=1.
- n=4: 3 nodes, c_B=2, bicycle=0, two-rank=1; Psi=z; 2 NS-source ordered pairs checked.
- n=5: 10 nodes, c_B=3, bicycle=1, two-rank=3; Psi=20z^2+197z^3+599z^4+603z^5+165z^6; tree count 41096; exponent 10274; 18 NS-source ordered pairs checked.
- n=6: 34 nodes, c_B=23, bicycle=3, two-rank=25; tree count 9672890651856613234344767407709495492608; exponent 576549211255110098978565180761188; 726 NS-source ordered pairs checked.
- n=7: 272 nodes, c_B=185, bicycle=11, two-rank=195. No integer Smith-form or full voltage computation is claimed at n=7.

For n=6,7 the owner classification and SC labels are inherited frozen inputs, not independently reclassified in this audit. The independent classifier covers the complete n<=5 universes only. All atlas hashes are in frozen_inputs.json and match those reported by the producer.

The pure-blue regular n=5 vertex independently has mass three, three blue stubs, zero paired black degree and three distinct neighbors. This validates the essential hostile against the old unrestricted pure-blue-leaf statement. Even-order purity is cited only via the corrected current theorem.

As an extra audit observation, 10 of the 18 order-five and 174 of the 726 order-six NS-source ordered pairs have odd reduced resistance denominators, while every complete voltage vector has even denominator. These bounded counts are not promoted to an all-order claim.

## 4. Producer replay and pin-failure tests

The three producer scripts were copied unchanged into this audit directory and run here, so their generated JSON did not overwrite producer outputs. For probe, hostile and voltage_audit, ordinary and optimized stdout matched each other and the frozen producer stdout byte for byte; their regenerated JSON matched the producer JSON byte for byte. All ten original producer SHA256SUMS entries verify after the authorized wording correction.

The producer scripts record, but do not themselves enforce, source hashes. The independent checker closes that reproducibility gap by rejecting any changed source before computation. negative_pin_tests.py alters only temporary copies and verifies that both normal and optimized runs reject (1) a changed atlas file and (2) a changed producer proof. All four cases fail closed with the precise expected RuntimeError. No original input is modified by these negative tests.

The audit pins both the final producer files and the reviewed inheritance documents. THM-4739 and THM-4740 are not opened or used as mathematical inputs. Their reserved status is not converted into a dependency. No historical Smith-diagram conductance, black-connectivity assumption, global blue-connectivity assumption or complete metagraph classifier is imported.

## 5. Reproduction and final boundary

From any working directory, using explicit locations if the directory layout changes:

    python3 independent_audit.py --repo /path/to/math --producer /path/to/session20-metagraph-networks
    python3 -O independent_audit.py --repo /path/to/math --producer /path/to/session20-metagraph-networks

The supplied default layout is the two sibling stage directories beside math. Compare stdout; all reported checks are deterministic. negative_pin_tests.py uses that supplied layout. The copied producer replay scripts require SymPy; the independent checker does not.

Approved scope: elementary all-order deductions for the precisely defined componentwise reflection-paired integer network, plus the named finite-exact controls and information-loss hostiles. Open scope remains metagraph classification, all-order full/blue/black connectivity, complete bridge classification, alternative-tiling realizability, and any exact reconstruction from electrical invariants alone.
