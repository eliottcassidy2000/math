# Independent review of four supplied lattice manuscripts

Review date: 2026-10-22. All four PDFs were read in full, including references and appendices; extracted texts are adjacent. These are manuscript claims, not assumed established literature. All authors are listed as SRFP311T1 Collaboration.

## Inventory and decisive outcomes

- `extremal_lattice_96.pdf`: *An Extremal Even Unimodular Lattice in Dimension 96*, 5 pages, dated October 9, 2026. Claims full extremality.
- `dual_horizon_96.pdf`: *Dual-Horizon Spectral Rigidity in Dimension 96: Deterministic Cyclotomic-Symplectic Architecture, Rational Certification, and the Extremal Frontier*, 6 pages, October 2026; PDF metadata October 21. Explicitly leaves extremality open.
- `niemeier_genus5_restriction.pdf`: *On the Restriction of Genus-5 Siegel Theta Differences of Niemeier Lattices to the Jacobian and Hyperelliptic Loci*, 16 pages, October 2026; PDF metadata October 21.
- `algebraic_compression_rigidity.pdf`: *Finite Algebraic Compression in ADE Root-System Moment Problems: An Exact Four-Moment Rigidity Threshold*, 11 pages, October 2026; PDF metadata October 7.

Strongest verified positive: independent exhaustive enumeration confirms the claimed rank-43 threshold for the explicitly supplied four polynomial moments, including unique collision.
Strongest negative: the Niemeier manuscript's hyperelliptic nonvanishing contradicts Codogni's published theorem, applying to exactly the same five pairs. The old extremal96 claim cannot be promoted to a theorem: the later manuscript retreats to a partial-search statement, and the full no-small-vector certificate has not been independently reproduced.

## 1. Extremal96 and dual-horizon96

### Independent exact controls

The old paper p5 links https://srfp311t1.com/G96_reduced.txt . It was downloaded as inert numeric data. File SHA256 matches exactly:
`7a6a3b8f761a9f02e1a687be8040aaf403c7ac75c33d92f010f2dd131f17bc24`.
Our `check_matrix.py` verifies shape96x96, integer symmetric matrix, even diagonal, determinant1 by exact DomainMatrix arithmetic, and all96 positive rational LDL pivots. Thus this file really defines an even unimodular positive-definite rank96 lattice. The initial pivots are28,21,56/3 and last three63/1742,22/567,1/22. These belong to the OLD reduced matrix, not the newer unreduced hash.

No complete shortest-vector enumeration was run. A bounded floating-point-assisted heuristic and an LLL control find no counterexample; this is NOT evidence establishing minimum10. Smallest basis norm20 is only an upper bound on the minimum, not a lower bound. Do not confuse it with certification. Old manuscript p4 item6 asserts qfminim(G,8) returned zero, but no complete search log/tree certificate is supplied in the PDF.

Exact truncated q-series computation (`controls.py`) gives the unique weight48 candidate with q1..q4 zero:
Theta = E4^12 -2880 E4^9 Delta +2324160 E4^6 Delta^2 -465638400 E4^3 Delta^3 +5611550400 Delta^4.
Its q5 coefficient is565866362880, q6 is2972108280960000. Old paper p2 equation2 instead states184139089776640 at q6, which is wrong; later paper p2 equation1.5 is correct. Our finite computation verifies positivity through q12 only, NOT all coefficients. Even all-coefficient positivity would be a necessary modular consistency condition, not existence of a lattice realizing the formal series.

### Proof audit

Old p3 Proposition3.2 proves a bound only on pure tensors and then claims it for all lattice vectors: decomposable-vector estimates do not bound arbitrary sums of tensors. Old p4 Proposition4.2 asserts indecomposable gluing contains no pure tensors, an implication requiring proof and false in general as a tensor-subspace assertion. The forbidden-line count of at most about10^12 has no derivation or reproducible enumeration in the text.

Most importantly, its discriminant bilinear form is symmetric in odd characteristic7, NOT alternating symplectic. In the split orthogonal48-dimensional case the number of maximal totally isotropic24-spaces is 2 product(i=1..23)(7^i+1), not product(i=1..24)(7^i+1). The fixed isotropic-line inclusion probability is1/(7^23+1), not1/(7^24+1). Nonisotropic lines cannot occur in an isotropic gluing space. This invalidates the particular counting argument, though correcting this factor alone would not defeat an independently proved small forbidden-set bound.

New p3 Theorem2.1 odd-prime gluing proof is sound: integrality plus49 norm(x) even forces norm(x) even. Theorem2.3 is a useful evasion equivalence on isotropic lines: their coset norms in an integral overlattice are even; checking complete forbidden classes is essential. It does not enumerate that forbidden set.

New p4 relative integral basis argument and displayed Schur identity7A-RR^T=I are directly checkable and consistent. Positive HJ implies positive Schur complement HJ/7, so the adjective 'indefinite hyperbolic' in the abstract should not be taken as a proven signature statement. HJ's alpha_pos twisting element is not explicitly specified in the PDF's Proposition3.3, nor are all96 pivots or the full new matrix given. Thus the new paper is not self-contained for reconstruction despite its certification language.

New pp5-6 sublattice minima/horizons are asserted computational certificates; they were not independently rerun. The completed-square reduction in Theorem6.4 is mathematically valid conditional on exhaustive Schur-candidate enumeration and exact lifted-slice test. Bridge Participation Theorem7.1 is simply the correct contrapositive: a short vector cannot have all bridge coordinates zero. It does NOT exclude short vectors with bridge support. Remark7.2 explicitly leaves minimum in{2,4,6,8,10}. Its matrix hash differs from the old paper, so do not silently treat the papers as certifying identical coordinate data.

The advertised verify_lambda96.py download returned an empty file in this check; no external script was executed. The old actual matrix remains a useful verified object, not a verified extremal object.

## 2. Niemeier genus5 restriction

### Established literature resolves the main restriction question

Primary source: Giulio Codogni, *Hyperelliptic Schottky Problem and Stable Modular Forms*, Documenta Mathematica21(2016),445-466. Publisher https://ems.press/journals/dm/articles/8965467 ; full paper https://emis.de/ft/9558 .

Theorem6.1, printed p461, says rank24 even unimodular forms with equal root count have theta differences vanishing on the hyperelliptic locus in EVERY genus. The paper lists precisely these five pairs. Since root count=24h, it applies to all h=6,10,12,18,30 here. Therefore attached Theorems4.9/4.11 pp12-14 cannot establish actual hyperelliptic nonvanishing. Table5.1's23 hyperelliptic classes is false: there remain19. Distinct h remain distinct by elliptic boundary specialization.

Codogni Theorem6.4, printed p464, additionally states the equal-root-count difference is zero on Mg for g<=4 and cuts a slope12 divisor on M5. For the distinct Niemeier pairs this supplies the known nonzero genus5 Jacobian restrictions. Thus genus5 Jacobian separation is compatible with established mathematics, but the manuscript's proposed hyperelliptic route is not. This is an excellent clean example of ambient information disappearing on a sublocus even when it survives on the larger Jacobian locus.

### Independently checked finite claims

`root_counts.py` constructs standard E8 roots with integer coordinates scaled by2, restricts to E7/E6 via orthogonality, then independently counts5-cliques with pairwise inner product4 (unscaled1). Dividing unordered star-frame counts by12 yields A5 subsystem multiplicities E6=36,E7=1344,E8=40320. Thus the unusual E7 count in Theorem2.6 is independently verified, not assumed from its normalizer argument.

The classical Dn D5 count binomial(n,5), An A5 count binomial(n+1,6), and Dn A5 count32 binomial(n,6) follow the coordinate descriptions. The nonembedding of A5 inD5 via determinant ratio6/4 is sound. The ambient nonzero Fourier coefficients therefore provide meaningful certificates; ambient nonzero does not imply restricted nonzero.

### Specific defects and qualifications

- pp11-12 Hypothesis4.6 is not implied by Schottky vanishing on J4. Restricting a Fourier series can produce cancellations among nonzero coefficients; it does not make the coefficients themselves vanish. The known hyperelliptic theorem forces a failure of the hypothesis and/or subsequent coefficient extraction, not a new exception to Codogni.
- Hypothesis4.6 concerns rank4 root configurations with zero-sum pairs/triples; the proof then discards all rank4 contributions without proving all relevant rank4 configurations or analytic corrections are covered.
- p14 Theorem4.11 identifies a multi-valuation coefficient with star-frame contribution without proving no other configurations share that valuation. Having a nonzero ambient coefficient is insufficient after substitution into constrained period matrices.
- p12 Component-Partition Theorem is correct unconditionally as a statement about rank5 six-root zero-sum cycles inD4^6. It does not control lower-rank contributions to the full degeneration coefficient.
- p10 rank<=3 annihilation follows from the Siegel operator and is sound. p9 degeneration genus calculation6-2+1=5 is sound. The crucial issue is loss of terms/cancellation after specialization, not these elementary facts.
- p3 equation1.4 has a normalization conflict with the companion compression manuscript's explicit definition J8=Theta(E8^2)-Theta(D16+): for h30 direct-sum multiplication gives F30=Theta(E8)J8 exactly, while15*m30/4=3360. A different J8 normalization must be explicitly declared to reconcile the formulas. In compression p9 the explicit definition and factor occur together and are inconsistent.
- h30 nonvanishing on generic J5 and vanishing on H5 are established. Primary Grushevsky-Salvati Manni source: https://arxiv.org/abs/0809.1391 , full https://arxiv.org/pdf/0809.1391 . Their result identifies the trigonal zero locus. The manuscript's exact multiplicity2 assertion was not independently verified here; do not promote it merely from the correct set-theoretic statement. Hyperelliptic curves have a degree3 pencil WITH a base point; distinguish trigonal closure/Brill-Noether locus from actual base-point-free trigonal curves.

## 3. Algebraic compression and exact ADE moment threshold

The polynomial moments pp3-4 are taken as DEFINITIONS in our computation. The paper does not provide a full derivation from an independently defined frame generating function; consequently the tested theorem is a theorem about these specified polynomial invariants of formal ADE multisets. It is not yet a theorem about arbitrary lattice theta coefficients without that additional identification.

Independent `audit_ade.cpp` enumerates every canonical nondecreasing multiset of86 types through rank43 and sorts exact integer4-signatures. It is a separate enumeration algorithm rather than downloaded authors' code. Results:
- ranks1..42:5,565,832 configurations, no repeated signature within a rank;
- rank43:1,556,185 configurations, exactly one duplicate, signature(318,-4110,268824,-2433804);
- cumulative7,122,017, agreeing with paper pp6-7.
The third stored coordinate is3c3. `controls.py` verifies the claimed pair with c3=89608 and fifth-moment separation158592. This independently certifies the finite polynomial rigidity claim and its minimality.

R1=2A1+2A4+A5+2A7+A8+D6;
R2=2A2+A3+2A6+A9+3D5.
Both ranks43. c5(R1)=93582144/5; c5(R2)=94375104/5.
The Coxeter-normalized defect is exactly-10h(h-2)(h-4), verified symbolically. The Euler-Boole finite-difference proof pp4-5 is sound: degree<=M+1 polynomial invariants annihilate under the(M+2)nd difference and produce equal-rank examples at(M+2)2^M. This yields rank96 atM4 for the restricted construction, not a relation proving existence of a rank96 extremal lattice.

Imported sectors pp7-9 remain imported, not consequences of the ADE computation. The Leech commutant needs the omitted actual character decomposition/rationality certificates. A rational commutant only forces rational scalars for a rational equivariant operator; a real equivariant Hessian needs rationality of that operator separately. The screening-ratio arithmetic reduces correctly to730/3597 if the input eigenvalues are correct, but no mechanism links that reduction to repunit rarity. The claimed weighted15-design uses two radii and therefore needs explicit normalization to a common sphere and harmonic verification; the paper only imports it. These should not inherit the 'certified' label from our ADE audit.

## 4. Useful research integration, without numerology

1. Treat generating functions as maps with explicitly stated domain and kernel. The exact ADE threshold is a rigorous finite model of how four aggregated moments can appear faithful over a large search range and then fail at a minimal collision. Transfer the audit methodology to tournament invariant candidates: enumerate isomorphism classes, record full coefficient signatures, exhibit first collision, and test the next refinement.
2. Distinguish ambient invariants from restricted evaluations. Niemeier forms give a particularly strong cautionary example: even all genera restricted to the hyperelliptic locus miss five pairs. Generic Jacobian genus5 succeeds by a different theorem. Any tournament specialization or lattice-to-modular projection requires an injectivity proof on the intended class.
3. Formal modular feasibility does not imply geometric realization, just as matching several tournament/deck statistics does not reconstruct an object. The weight48 formal theta is uniquely fixed and arithmetically consistent through the checked range; constructing a lattice requires exact Gram certification plus exhaustive minimum proof.
4. The number24 is structurally relevant through rank24 Niemeier classification, weight12 theta series, and the dimension24m bound. Rank96 in the Euler-Boole construction is a different mechanism; equality of numbers is not an identified mathematical bridge.
5. For1729=7*13*19 and repunits, a legitimate connection must run through an explicit multiplicative-order/cyclotomic identity or a proved invariant correspondence. The appearances of7 in this lattice construction, and7*11*13=1001 in an imported fraction reduction, do not supply such a theorem. Keep repunit rarity quantitative and conditioned on modulus/order, not inferred from these coincidences.

## Reproduction

Compile and run: g++ -O3 audit_ade.cpp -o audit_ade; ./audit_ade .
Run python controls.py (SymPy), root_counts.py (standard library), check_matrix.py (NumPy/SymPy; local downloaded old matrix required).
`ade_audit_results.txt`, `controls_results.txt`, `root_counts_results.txt`, `matrix_results.txt` preserve outputs.
The short-vector/LLL heuristic scripts are exploratory only and carry no exclusion certificate. No original PDFs were uploaded to Git, and no repository operations were performed.

## Addendum: explicit algebraic mechanism for the rank43 collision

Let d_h be the difference, R2 minus R1, of total component rank carried by Coxeter number h. The nonzero masses are d2=-2,d3=4,d4=3,d5=-8,d6=-5,d7=12,d8=1,d9=-8,d10=3. Their generating polynomial factors exactly as

P(z)=sum_h d_h z^h=z^2(z-1)^5(z+1)(3z^2+4z+2).

Consequently sum_h d_h h^j=0 for every j=0,1,2,3,4. The per-rank A-series cumulants through order4 are polynomials of degree at most4 in h. Their contributions therefore cancel by this fifth-order finite-difference identity. For orders1-3 the D per-rank values coincide with A. At order4 the remaining corrections cancel independently: R2's three D5's give15 units of rank at h8; R1's D6 gives6 at h10; and15*delta4(8)=6*delta4(10)=-28800. Thus the full4-moment collision has a short exact explanation, not only a machine witness.

The factorization proves this particular collision and illuminates its finite-difference mechanism. It does NOT prove global minimal rank43 or uniqueness: those are the exhaustive audit's separate results. No assertion that this11-component-support relation is a linear-matroid circuit is made; there are more columns than moment constraints. Exact symbolic output is in collision_mechanism_results.txt.

### Completed cheap hostile checks on old96 matrix

Every signed unit-coordinate vector with support2,3,or4 was checked exactly (up to global negation), together with all basis vectors. None has norm below20. A separate bounded rounded-eigenvector heuristic and numerical Gram-LLL control likewise found no norm<=8 witness. These checks do not disprove or establish min10: they leave the claim uncertified. They must never be described as exhaustive lattice enumeration, since arbitrary integer coordinates and larger support remain untested.

The fifth-order collision split is explicit: A-proxy contribution72576 plus D correction86016 gives158592. The normalized fifth-order defect is(28/5)h(3h^3-24h^2+67h-58).
