# Independent review: self-converse strong C2 classification

**PASS, no correction.** This review is bound to the frozen self-converse-c2 packet below. All eleven files in its RELEASE_FILES.sha256 were checked again after the audit and remained unchanged.

## Binding

- PROOF.md: 3e6d6e7344413504f299329fec05304600deef9432911f7aa9fecbd8fd4cb908
- converse.py: ea9ef882db60d29a4152630c4f5475999eecb65f67afeb8c65241eefb7e35c2f
- verify.py: ffc9235ceee0b80b3d3f04d08157fc093a68191ab5f678926a25df449a4a0b1d
- core.py: dc32dc1ef3647b694015cb4ff1ef99415c162bfbaf58d3b545175b09f0378042
- bilinear.py: 33cfe6cea487b17da44e083d7e81bd860dfeddcab3096deb6d135191e5537ac0
- All three source validation JSON files: 50ad3f92bb5228f015b5d8fd416c57d786603f69dd4809aff51505c6fc3f6472
- BILINEAR_C2.md: 5eab9deb2871819d8a1b5684c4987f2f4da807adc55063a30675d531c5ad3813
- C2_COMPLETE.md: 4cc0bf11efc435d7a2746daf6d8f1a58c49c94e57e608911b24464facdd1057d
- FIBONACCI_CLASS.md: d95a133e360b96589009ebd8a1c8449525c4a2847a88b33a4b3f2b839739e63e

## All-order proof review

The converse maps check every type of incidence. For N the central vertex is fixed under source exchange; for M the two central vertices are exchanged. Within-source arcs reverse, the root-root exceptional seam reverses into the corresponding seam, and the remaining central edges have the stated paired incidences. The displayed index maps therefore prove the identity for arbitrary actual sources, not merely equal response quadruples.

The named bilinear theorem gives an ordered graph-class bijection, not merely a scalar surjection. Its injectivity follows from its explicit surjection and the independent complete class generating function. Applying this exact statement to F(V,W) and F(W,V) proves self-converseness precisely when the rooted source classes agree. Both the regular-five exception and the nontrivial-source restriction are retained. Rigidity of every other C2 graph is explicitly present in the named C2_COMPLETE dependency, so the uniqueness and fixed-point claims for the anti-involution follow.

The diagonal generating function doubles one chosen rooted source, then adds one or two central vertices. It is not the generating function of unordered arbitrary pairs. Its expansion and the exceptional labelled contribution24 are correct. The existential scalar image legitimately uses numerical quadruples: any actual source realizing such a quadruple may be used twice, even though two different source graphs can share that quadruple.

The diagonal lower bound3U²+1 makes the inverse finite. U bounds physical source order; the inherited strictly-U-increasing oracle supplies every required actual state and a recipe. The sharper M bound is valid but unnecessary for completeness. No finite-order extrapolation enters the inverse.

The parity induction is valid for O, for L=1, and for L>=2 separately. D is odd and E+U is even, so the deletion term vanishes modulo4 in both diagonal formulas. The remaining terms give1 modulo4 in both parity cases for E,U. The conclusion excludes only self-converse C2 graphs outside that residue, with R5/H15 expressly excepted; it does not assert an unrestricted Hamilton-count obstruction.

## Independent physical controls

The standalone C++ checker does not use the source source-grammar, cycle-canonicalization, cycle-matching anti-isomorphism routine, or Hamilton-count helper.

- It exhausts all33,867 labelled tournaments of orders1 through6, using a backwards subset path recurrence only to identify every root with start count one. There are26,175 such rooted inputs. Their unique root-start Hamilton words give canonical rooted arrays, yielding exactly1,1,2,5,13,34 classes by order.
- From that independently reconstructed bank it builds686 ordered N/M graphs through order10. The bank covers all presentations through order9; at order10 it deliberately omits sources of order7, so it is not advertised as a complete order10 graph census. All self-converse diagonal presentations through order10 are included.
- Literal path DFS, with no memoization or endpoint-count formula, counts H and cycles in all686 graphs. It agrees with the bilinear scalar formula and C=2. Every source-swap map is checked edge by edge.
- A separate degree-compatible exhaustive bijection search computes the complete automorphism and anti-isomorphism sets. It finds rigidity and the exact diagonal criterion in all686 graphs. Including R5 and the two hostile graphs, there are689 full automorphism and689 full anti-isomorphism checks.
- R5 has H15, C2, and five automorphisms and five anti-isomorphisms. The complete self-converse coefficients at orders5 through10 are2,1,2,2,5,5.
- The two specified order7 sources are independently canonicalized and have distinct rooted arrays but the same response(27,9,27,9). Literal DFS gives H3321 and C2 for both the mixed and diagonal order15 constructions. Exhaustive bijections find zero versus one anti-isomorphism respectively; both graphs are rigid.

The separate algebra checker computes all57 gadget constants by literal path enumeration, rather than the source g/Q recurrence. A strictly-U-layered enumeration then reproduces all159 physical response states with U<=57. This yields the same194 diagonal values through10000, independently checks1000 exact inverse targets, verifies the positive bounds/residue conditions, and checks101 coefficients of the rational generating function.

## Reproduction and receipts

Compile check.cpp with C++17 and run it. Optimized (-O2) and debug (-O0 -D_GLIBCXX_ASSERTIONS) outputs agree byte for byte:

- check.cpp: 848e2165a1634051f1feb16255f14adef7592cd6cd90099b74e05e7a066c2df1
- independent-normal.json and independent-debug.json: 7f478b21fea40f55982a39eb8427fcc7167c7bcc4620e27a83cd81d3132a449c

For check_algebra.py place an unchanged copy of the reviewed source packet in the adjacent source/ directory. Its imports are used only for comparison with the source oracle/inverse; its gadget constants and layered bank are computed independently. Normal and -O outputs agree:

- check_algebra.py: 095455e0ec0b041686a025781f70007f72cc59c7de51ad123d3442ba38ac02d7
- algebra-normal.json and algebra-optimized.json: ac69e3d763c286ac7634e10e754b63e203f4ae859b18c2b04addd5addf3377e9

The original source verifier was run normally and with -O only in that audit copy, because it writes its receipt. Both reproduced the frozen validation bytes exactly (50ad3f92...). All mathematical checks use explicit failures rather than removable Python assertions.

This review confirms the stated classification and inverse conditional on the named already-reviewed full C2 classification and source grammar. It makes no new cofinite-coverage or literature-priority claim.
