# Session 55: Parker fields, faithful lifts, and residual score constraints

25 October 2026. Final research packet. The integer magic-square-of-squares existence and least-height questions remain OPEN. No computational record or external priority claim is made.

Convention: square residues include zero; arrays are labelled and all nine values must be distinct unless stated otherwise.

## Results and reading routes

- [Finite-field theorem and classification](finite-fields/THEOREM.md), [independent audit](audit-finite-fields/AUDIT.md): every characteristic-two field is Parker. The complete odd-order exception set is {3,5,7,9,11,13,17,19,23,25,27,31,43,47,67,243}. The proof combines a six-factor character estimate giving existence for odd q>=14400 with independent exact checks of all 1730 smaller odd prime powers. This reproduces the classification in a located 2026 Lai paper, not a first-classification claim. The character estimate depends on the cited Weil theorem.
- [Integer-lifting search](integer-search/RESULTS.md), [independent audit](audit-integer-search/REPORT.md): prime 59 is the first for the WEAK canonical no-wrap quadratic-residue predicate. Its grid has ordinary sum 48 but only two actual integer-square entries. Faithful no-wrap lifting is equivalent to the open integer problem. A complete gap-based search through center root 100000 is a reproducibility demonstration, far below historical bounds.
- [Exact six-coordinate carry equation](carry-lattice/PROOF.md), [independent audit](audit-lifting-scores/AUDIT.md): the integral linear magic lattice splits; modular root carries must cancel the residue-grid syndrome. Linear modular equations alone do not preserve actual integer squareness or height.
- [Residual tournament score constraints](residual-scores/THEOREM.md): removing a Hamilton cycle gives the user's total n(n-3)/2 and exactly n/2 high residual scores for even n. Removing the fixed Hamilton path instead leaves one additional edge, so the same two-value prescription has n/2+1 highs. All labelled high sets are realizable in the stated all-size ranges, with the sharp order 4 cycle exception; internal balanced ports preserve both endpoints. These are score-support statements, not full wiring or equal realization counts.
- [Recovered connections](recovered/ATLAS.md): common Gaussian norm representations give an exact four-offset closure test, not automatic magic completion; center root 1105 has 13 offsets but no required diamond. Square-sum all-order existence is a CITED Gerbicz result with full construction unreplayed here, not a conjecture inferred from a short window. Extension fields and quotient rings have different reduction maps.

## Quotient-ring classification and exact counting

- [Complete quotient-ring classification](crt-wiring/CLASSIFICATION.md), [independent audit](audit-odd-rings/AUDIT.md): exactly 136 moduli n>=2 are Parker, with largest 3216. Exactly 28 odd moduli n>=2 are Parker, largest 201. Including the trivial modulus 1 changes these counts to 137 and 29. The reduction covers every positive modulus, leaving only 203 explicitly certified candidates; it is not an extrapolation from searching small n. This establishes the bounds stated as Cain's Conjectures 9.3 and 9.4, with the audited prime-field theorem as an explicit dependency. External priority is not claimed.
- [CRT compatibility theorem](crt-wiring/THEOREM.md): keep the actual repeated-cell relation on nine fixed labels. Component masks intersect under CRT; global distinctness means empty intersection. Two Parker components can combine into a non-Parker product, first at 35=5*7. Minimal collision-mask antichains preserve feasibility, not multiplicity.
- [Prime-power lift towers](padic-lifts/PROOF.md), [independent audit](audit-padic-lifts/AUDIT.md): Z/2^kZ is non-Parker iff k>=7, and Z/3^kZ iff k>=3. Explicit fixed positive ordinary magic grids have square residues through arbitrarily deep towers without having nine actual integer-square entries. Unit-root solutions at an odd prime have exactly p^3 children per level.
- [Exact square-count generating functions](prime-power-counts/PROOF.md), [independent audit and ordinary-modulus extension](modular-count-audit/PROOF.md): if N_p(k) counts labelled distinct square-residue magic arrays modulo p^k, the series are

  N_2(z)=1536z^7/[(1-z^2)(1-2z)(1-4z)(1-8z)],

  N_3(z)=216z^3(1+27z)/[(1-z^2)(1-3z)(1-9z)(1-27z)].

  The normalized limits N_2(k)/2^(3k) and N_3(k)/3^(3k) are 1/504 and 27/728. They count grids of residues, not integer squares, root choices, scaling classes, or tournament classes.
- The ordinary distinct-magic-array count modulo any n has a proved degree-three quasipolynomial of minimal period 60, obtained from the exact Smith invariants of the collision arrangement. This ordinary count omits the square-entry restriction; see the [formula and boundary cases](modular-count-audit/PROOF.md).

## Final extensions

Read the [session synthesis](SYNTHESIS.md) for the mechanisms, lost coordinates, and remaining integer frontier.

- [Dyadic re-entry](dyadic-reentry/THEOREM.md): the complete classification has a short seven-layer profile. Exact local identities explain non-Parker 50 -> Parker 100 -> non-Parker 200. Parker status is not monotone under taking multiples.
- [Weighted CRT counting](weighted-crt/THEOREM.md), [independent audit](audit-weighted-density/AUDIT.md): equality-partition histograms convolve under meet; the zeta transform makes CRT multiplication coordinatewise. Counts at 35 and 55 are 24 and 40, or 3 and 5 dihedral position-orbits. Dropping existence-dominated states loses multiplicity.
- [Field/ring density bridge](density-bridge/PROOF.md): odd-field labelled counts are q^3/512+O(q^(5/2)), proved uniformly by character sums rather than independence. Comparing exact dyadic/ternary ring limits yields 64/63 and 1728/91, with their actual dimension/content mechanism stated.
- [Ordinary magic generating function](ordinary-count-gf/PROOF.md): an explicit reduced rational function has cyclotomic pole orders 4,3,2,2,2,2 at orders 1,2,3,4,5,6. Its coefficient quasipolynomial has minimal period 60 despite no primitive 60th-root pole. Square restrictions are NOT included in this ordinary count.
- [Collision sieve](collision-sieve/RESULTS.md), [independent audit](audit-collision-sieve/AUDIT.md): the exact modular filter removes 347087 of 347561 candidate pairs, but is slower than the direct difference-gap lookup on this benchmark. No speedup is claimed.
- [Exact SVG diagrams](diagrams/README.md) show the fixed-label CRT repair, field/ring contrasts, and dyadic re-entry, with deterministic source-pinned replays.

## Scope boundary

Finite-field classification, integer existence, native metagraph wiring, and square-sum Hamiltonicity are different predicates. The exact bridges retain root coordinates, carries, labelled endpoints, or the specified additive closure pattern. Matching small-size transitions alone is not a proof transport.

## Verification

Each packet contains its executable checker and exact outputs. The sibling audit packets give independent arithmetic/cut/root enumeration paths. Publication excludes private full third-party paper text, caches, and startup transcripts. The final session status is recorded separately by the verified publication receipt and explicit end message.
