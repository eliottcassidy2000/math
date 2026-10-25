# Session57 wave1: recovered square-sum / magic-square interfaces

Scope: repository inheritance at commit 73d9abcd5ddc516a4142079bd0ccd7c4022d10fa. Status labels distinguish recovered results, elementary new consequences, untested possibilities and explicit failed transfers. This packet makes no integer magic-square existence claim, literature priority claim, or search-record claim. SOURCE_PINS.json records exact paths, hashes and permanent source links. Source labels below refer to those snapshots.

External anchor: Onno Cain, *Gaussian Integers, Rings, Finite Fields, and the Magic Square of Squares*, arXiv:1908.03236v2 (2019), https://arxiv.org/abs/1908.03236v2. Its abstract identifies a quartic/factorization reformulation and a Gaussian search specialization, and distinguishes finite fields from quotient rings. Precise internal results below are recovered from repository proofs, not attributed wholesale to Cain.

## Immediate recovery ledger

- S07 already proves the exact target: D_e={d>0: e²−d and e²+d are integer squares}; find x>y>0 with x,y,x−y,x+y in D_e and x≠2y. The grid is e² plus the deviations (x,−x−y,y; −x+y,0,x−y; −y,x+y,−x). All nine values, not signed roots, must be distinct.
- S01 already classifies every positive AP triangle in the square-sum graph: 2g²(C²−4AB,C²,C²+4AB), with A²+B²=C² primitive and C²>4AB. This is NOT an AP of square vertex labels.
- S02 counts unfiltered Gaussian/Pythagorean petals by (product_(p=1 mod4)(2v_p(H)+1)−1)/2. Its numerical-height obstruction concerns many petals sharing one center; it is not a bound for arbitrary square-sum paths.
- S04 recovers exact endpoint-marked two-path-cover birth gluing, the realized symmetric-difference hitting-set characterization, and a cotree parity decoder with degree/connectivity gates. S05 records first endpoint-plus-ordered-birth-color ambiguity at order36.
- S03 identifies the three Q_12 islands as 3,4,5 vertices with founders 2,4,1 respectively; closure holes {1,4,6} are a different object. Q_n stays connected from14 onward, so there is no recurring three-island phase.
- S08–S11 already retain named equality partitions for CRT, full multiplicities for counting, and shared original choices for overlapping joins. Existence antichains, partition histograms and correlated assignments solve different retention problems.

## Ranked connection board

Rank is by cheap discriminating test and utility, not probability of solving the integer problem.

### 1. Exact short-gap-to-AP-petal correspondence — PROVED here from RECOVERED formulas

For fixed e>0 and d in D_e, map

    d -> T_e(d)=(2e²−4d, 2e², 2e²+4d).

The three pair sums are 4(e²−d),4e²,4(e²+d), hence squares. This triangle is positive exactly when d<e²/2. Conversely a positive AP square-sum triangle centered at 2e² has difference4d for an integer d in D_e, by S01: e=gC and d=2g²AB. The map is a bijection between short gaps D_e∩(0,e²/2) and positive AP square-sum triangles with that center. Distinct d gives distinct petals.

Preserves: exact roots, common denominator/center root, gap identity. Loses if only the petal stock is retained: all long gaps d>e²/2. The boundary d=e²/2 cannot occur for positive integer e because e²/2 cannot be an integer square. All short gaps form only a subproblem: a diamond there would give a magic square, but no proof says every diamond can be moved there. A common-center petal Hamiltonian path does not impose additive diamond closure.

Cheap test: compare independently enumerated D_e to AP triangles via their roots; retain positivity explicitly. e=5 has D_5={24}, but no positive petal. e=41 has D_41={720}, giving (482,3362,6242). The common divisor-capacity product is explained by this map, not a numerical coincidence.

### 2. Opposite-pair reflection versus square-shell reflection — PROVED obstruction / OPEN repaired transfer

A magic grid with square center e² has opposite-pair sum2e². This is not an integer square for e>0. Therefore every opposite pair is a NONEDGE in the square-sum graph on the nine entry values. A bridge that directly declares the four opposite pairs to be square-sum edges is REFUTED universally.

Both carriers have reflections t->s−t, but the permitted shell differs: s=2e² for magic opposites, s=k² for Q_n edges. Retain the actual shell and its square class. The two-step Q walk translation is l²−k²; forgetting which shell acts is the lost coordinate. A repaired graph may use typed edges or lift to the petal map in route1, but must prove new incidence statements rather than inherit Hamiltonicity.

Cheap test: the nonsquare shell proof eliminates this proposed direct route before search. It also prevents identifying AP-petal center2e² with magic center e² by ordinary square rescaling.

### 3. Diamond hyperedges, not arbitrary Hamiltonian paths — RECOVERED iff / OPEN algorithmic bridge

Make a four-uniform hypergraph on D_e with edge {y,x−y,x,x+y} whenever x>y>0 and x≠2y. The integer problem at e is precisely nonemptiness of this hypergraph. Preserve which gap occupies each named role, or canonically normalize x>y. A graph of additive compatibility can accelerate discovery, but unlabelled edges, high degree or a Hamiltonian path need not force a common four-role witness.

Cheap test: S06's e=1105 has13 gaps and no diamond, so capacity alone is refuted. Compare a proposed graph criterion against this stock before extending it. Search membership in x−y first, then x+y: S06's benchmark reports this already beats its additional Parker collision sieve. Do not re-advertise that sieve as a speedup.

### 4. Named-port path gluing and named-cell arithmetic gluing — RECOVERED theorem / OPEN common implementation

Path birth at v is H_(G+v)=sum_a x_va E_a(G)+sum_{a<b}x_va x_vb F_ab(G), where F_ab counts two path components with a,b in DIFFERENT components. CRT combines local arrays by equality-partition meet on the SAME nine cells. Both operations require owners, but the owners differ: component endpoints versus cell identities.

Exact shared language: relations of assignments on named variables, natural join, and once-only weights (S11). Preserve square-root/carry variables as well when moving from residues to integers. Lost by marginalizing first: component ownership, actual CRT choices, or cross-block equations. A join tree is sufficient only with running intersection and actual local constraints; an arbitrary collection of projections can omit a global parity constraint.

Cheap test: route candidate interfaces through S04's disjoint edges1–24 and4–5 at order25, individually used4 and3 times but never jointly; and S11's [2,5,3] overlap hostile, where true mass112 falsely becomes784 and fabricates128 distinct rows. A proposed implementation should reject both.

### 5. Sensor banks and arithmetic collision sieves — RECOVERED / OPEN restricted compression

For fixed endpoints, differences of realized Hamiltonian paths form a hypergraph of symmetric-difference edge sets; a coordinate bank decodes iff it hits every hyperedge (S04). For CRT, each retained modulus/choice must separate all36 cell pairs globally; selected collision masks intersect to zero (S08). Both are finite relational hitting problems, but one hits path differences and one eliminates surviving equalities.

Preserve actual edge/cell labels and the family of REALIZED candidates. Lost by using only block sizes, degree sequences or arbitrary cycle-space elements: which collisions occur together. Cheap test: {18,13} endpoint ambiguity at Q_26 is resolved by bit3–6; same-sized CRT partitions A,A fail while A,B succeed at35. OPEN: whether a useful candidate-conditioned sensor bank cheaply selects root-carry checks. Such a bank must be certified on the chosen candidate family, not promoted to a universal arithmetic theorem.

### 6. Gaussian norm stock versus additive closure — RECOVERED mechanism / OPEN obstruction

S07 gives e=k(m²+n²), d=4mn(m²−n²)k² and |D_e|=(product(2alpha+1)−1)/2. This is a norm-representation enumeration. The desired x±y memberships are additive relations among the resulting norm differences, not consequences of multiplicative prime allocation.

Retain the actual signed Gaussian angle/exponent data until evaluating the gaps; count-only projections lose the additive geometry and, for route1, the positivity angle. Cheap test: e=1105, then equal-capacity centers with different short-gap counts. OPEN: prove a restriction on additive diamonds in specific Gaussian prime-support families. No arbitrary norm-factorization-to-Hamiltonicity transfer is currently justified.

### 7. Affine-plane partitions versus projective polarity — RECOVERED exact limited map

The two three-block partitions used at modulus35 are complementary coordinate partitions of nine labelled cells, expressible as affine F3 coordinate fibres (S06,S08). Their meet is discrete because each block intersection has one cell. This is a legitimate affine incidence interpretation.

S13's F64/F4 trace polarity instead has21 points and21 lines and requires a typed sheet switch. It is not this nine-cell geometry. Forgetting sheet type previously broke polarity; forgetting shared cell names breaks CRT. Their common lesson is a loss ledger, not an isomorphism. Cheap test: incidence cardinalities9 versus21 and the actual binary operation. Four affine parallel classes on nine points require four hierarchy trees in S12; one arbitrary LCA tree cannot preserve all crossing partitions. OPEN only after an explicit operation-preserving map is supplied.

### 8. Linear forests, path covers and degree fibres — RECOVERED / REFUTED count conflations

S14's tournament linear-forest expansion is a signed arc-subset identity for Hamiltonian path counts; reversing odd runs cancels their contributions. S04's square-sum two-path covers are unsigned undirected spanning objects with exposed endpoint owners. S12's score-fibre count is the number of unrooted spanning forests of K_n minus a locked Hamiltonian path, while det(I+L) counts rooted/component-weighted forests. These are different species.

Preserve orientation, signs, component count, roots and locked-edge ownership. Cheap tests: triangle rooted-forest count16 versus unrooted7 (S12); fixed-score tournaments with the same flux but11 versus13 Hamiltonian paths. Degree fibres are not isomorphism classes and degrees do not certify additive diamonds. OPEN: carry an edge-marked path-cover generating function into an arithmetic constraint factor graph, retaining all square-root and no-repeat constraints. This would be a generic exact solver interface, not a solved structural equivalence.

## Recommended next wave

Prioritize route1's lossless coordinate map and route3's labelled diamond incidence; test whether proposed path intersection statistics retain x−y and x+y jointly. Keep route4's shared-choice joins as the implementation discipline. Retire direct opposite-edge and same-center identifications immediately. Historical 3–4–5 islands, polarity and forest vocabulary should enter only when an explicit typed map survives their cheap hostile tests.
