# Independent audit: one-apex Hamilton spectra and bounded-run feedback grammars

Status: ACCEPT after the corrected README authorship. Mathematical statements in `proofs.md` and `run_density.md` pass proof-level audit. All reported finite controls pass independent replay. The global odd Hamilton-path spectrum conjecture remains OPEN.

Audit date: 2026-10-20. Source directory: `session16-hamilton-spectrum`. This report covers the separately frozen initial packet and bounded-run addendum, not a subsequently proposed intercyclic reduced-word grammar. No producer mathematical source, repository file, or remote was changed by the auditor.

## 1. Correction, scope, and packaging

The original README misattributed the 2026 special-tournament paper to Nino Samashvili. Independent primary-source checking identifies Maria Abi Aad through the author's SSRN record. The producer corrected the README and appended its hash to `SHA256SUMS`; that corrected README is approved. The original snapshot is retained only as superseded audit evidence, not as a current bibliography.

The mathematical source hashes were unchanged by this correction. `RUN_SHA256SUMS` does not itself seal `probe_affine.py`, although `verify_runs.py` imports its matrix builder. The full audit manifest therefore includes this dependency, both probes, all approved proof sources, and normal/optimized outputs. Publish the dependency with the verifier. This is a packaging observation, not a mathematical defect.

No claim is made that a finite computation proves an all-order theorem, that a connected metagraph supplies unit steps, that scalar products imply SCC decomposability, or that a finite-state numerical image is automatically complete.

## 2. Initial packet: proof-level verification

### Exact one-apex formula

After deleting the apex from a Hamilton path, the remaining vertices form two increasing blocks. Their last/first boundary conditions give the stated ordinary cuts and inversion contributions. The inversion with endpoints i<j has exactly 2^(j-i-1) assignments of its intervening vertices. The identity N01-N10=b_m-b_1 supplies the additional adjacent-inversion contribution, giving weights w(1)=1 and w(d)=2^(d-2) for d>=2. The empty-chain and one-letter cases are covered.

The four-coordinate update is consistent with those geometric weights; it retains the preceding bit separately because the adjacent weight differs from the otherwise geometric rule. All counted objects are ordinary labelled vertex-word paths, not isomorphism or automorphism orbits.

### Strong normalization and extremal boundaries

A leading zero is a source and a trailing one is a sink, so their deletion gives an actual count-preserving bijection. A nonempty reduced core begins with one and ends with zero and has the explicitly displayed spanning directed cycle. This proves the all-order reduction, rather than merely making it plausible from small words.

For a fixed number r of ones, every decreasing interior exchange strictly lowers H by the exact displayed difference. Bubble sorting reaches the unique minimum word 1 0^a 1^b 0; reverse sorting reaches the unique maximum word 1^r 0^s. Uniqueness is correctly restricted to fixed-content pointed words.

The product formula (2^a+1)(2^b+1) requires a,b>=1. At a=0 or b=0 the separate value 2^(m-1)+1 is necessary. At ordinary order five there are exactly three minimizing reduced words, 1000, 1010, and 1110, all with H=9; 1100 has H=11. The balanced interior minimum is strictly below both boundary values only from m>=5. Both global fixed-order extremal formulas and the order-three/four exceptions are correct.

The interior-swap lower bound contains both fixed endpoint contributions. For m>=5 at least one is positive, giving a jump of at least four. Thus the natural adjacent-exchange graph has no H-jump two even though each fixed-content fibre is connected. This does not exclude other moves or other metagraphs.

### The all-order H=13 omission

The minimum for strong one-apex cores at order six is 15 and increases thereafter. Smaller reduced cores have only H in {3,5,9,11}; transitive tournaments give one. Therefore no one-apex tournament of any order has H=13. The explicit unrestricted five-vertex tournament with rows (22,12,24,17,2) is strong and has exactly thirteen paths by literal permutation counting, so this is a construction-class obstruction, not a global missing value. Products and source/sink padding cannot first create the prime thirteen.

### The arbitrary fixed-feedback gap automaton

Fixing the feedback permutation partitions a path into k+1 increasing ordinary blocks. For each gap, emptiness and the current last-vertex validity are sufficient: the left endpoint is checked exactly once, on the first assignment, and the right endpoint is checked at termination using the most recent assignment. Empty internal gaps require the actual feedback arc. Summing transition multiplicities is indispensable and is done correctly.

Every path has one feedback permutation and one assignment of ordinary vertices to its gaps; conversely every accepted assignment yields one path. This bijection proves exactness uniformly in chain length. Taking the direct sum over k! permutations gives the stated upper bound k! 3^(k+1), including the k=0 degenerate case. The last gap actually needs only two states, but using three gives a valid upper bound. The code's two-feedback realization uses 36 states, not the looser 54-state bound.

### Poisoning-DAG hostile

The nine rows specify a valid tournament with exactly the six listed directed triangles. Two disjoint triangles exist, but three cannot, and every vertex belongs to a triangle. For the chosen matching, the remaining vertices are {2,7,8}; the poisoning arcs are exactly 7->8 and 8->7. Thus the historical DAG lemma is false, specifically because its asserted inequality w_k!=w_2 fails for a directed two-cycle.

The final safe-deletion conclusion is not refuted. In fact vertices 1,2,4,5,6 are all safe deletions in this witness. Its Hamilton-path count is 75 by literal enumeration of 9! permutations. The correct scope is a failed proof step with an explicit hostile, not a refutation of the final dichotomy or the permanent omissions.

## 3. Bounded-run addendum: proof-level verification

### Commuting gap projectors and exact spectrum

For a fixed letter, assigning to one gap is an idempotent partial deterministic map: either an initially forbidden empty gap is killed, or the gap is reset to the fixed terminal status; a second reset changes nothing. Operators on distinct gap coordinates commute even in the killed-state cases. The same identities hold blockwise over feedback permutations.

Commuting rational idempotents have a common rational eigenspace decomposition. Hence the letter matrix is diagonalizable and its spectrum is contained in {0,...,k+1}. The addendum correctly says contained, not that every value occurs. The square-free annihilator and absence of polynomial run-length factors follow. The zero-eigenvalue contribution at r=0 requires the stated convention 0^0=1.

### Positive surjection expansion

Expanding a power of the sum of projectors and grouping assignments by the set of used gaps gives exactly Surj(r,|J|) copies of each projector product. Every contextual coefficient is a nonnegative integer, including after summing over feedback permutations. This positivity, not an arbitrary exponential fit, rules out cancellation. It is the essential mechanism behind the effective compression bound.

### Compression and density

If some coefficient for t>=2 is positive and r>=k+1, fixing the first t choices surjectively and using either of two gap choices thereafter gives H>=2^(r-t)>=2^(r-k-1). Thus a run longer than k+1+floor(log2 H) must have all higher-support coefficients zero; replacing it by one copy preserves the value exactly. The prefix and suffix are frozen during each application and may be updated for the next run. This makes successive compression valid without any global stability assumption on coefficients.

The resulting finite encoding bound counts feedback tournaments, run-letter sequences, and run lengths. It proves O_(k,R)((log X)^R) distinct values up to X for fixed k,R, including the empty-chain case. It does not assume strong connectivity survives vertex deletion, and it deliberately does not cover unbounded run complexity.

The multiplicative-closure corollary is also correct: a prime in a product must already occur as one factor. The provided elementary prime-count lower bound is valid, because each prime-power floor difference in the central binomial valuation is 0 or 1, hence p^v_p(binomial(2n,n))<=2n. Together with the central binomial lower bound this beats every fixed polylogarithm. Infinitely many odd primes therefore remain absent from such a bounded grammar even after arbitrary products.

### Intercyclic reduction and its boundary

The cited Bang-Jensen--Bessy theorem really does give feedback number at most two for an intercyclic tournament, with a sharp example. Since every tournament cycle contains a directed triangle, the no-two-disjoint-cycles and no-two-disjoint-triangles conditions agree. The source provides no run bound. The frozen addendum correctly refrains from inferring one and makes no claim of a complete intercyclic numerical spectrum.

## 4. Reproduction and independent controls

Producer replays, each normal and optimized, byte-match the frozen output:

- `python3 source_snapshot/verify.py`: 511 one-apex controls, 4,097 exact adjacent exchanges, 256 product-family checks, and 1,272 k=0,2,3 general automaton controls.
- `python3 source_snapshot/verify_runs.py`: all twelve k=2 gap projectors, 36 commutators, four annihilators, homogeneous powers through length nine, 1,764 contexts, and 933 long silent-run checks.

The clean-room `independent_checks.py` imports no producer or repository mathematics implementation and uses Python arbitrary-precision integers throughout. Its normal and optimized JSON outputs byte-match. It checks:

- 1,335 literal permutation counts against literal gap assignments, with feedback vertices labelled before the chain rather than after it;
- fixed-content uniqueness and global extremizers through chain length ten, retaining all order-five minimizing words;
- the explicit strong H=13 witness and all six triangles, poisoning arcs, safe deletions, and literal H=75 for the nine-vertex hostile;
- for k=0..4, 22,395 idempotence identities, 108,615 commutation identities, 4,665 exact annihilator columns, and 984 positive Stirling row expansions;
- 2,510 prime-power valuation controls for the elementary prime-count lower bound.

These are independent hostile and boundary controls, not substitutes for the arguments in Sections 2--3. No large tournament census was run. The two negative affine-probe conclusions remain restricted to their explicitly specified finite context/block universes.

## 5. Independently checked primary sources

1. Bang-Jensen and Bessy, *Small cycle transversals in tournaments*, author-hosted preprint dated January 16, 2013. Theorem 1.8 on PDF page 2 supplies the sharp feedback-two bound; the proof begins on page 3. Source: https://www.lirmm.fr/~bessy/publis/SmallTransversTournament.pdf . Direct PDF retrieval succeeded.

2. Arthur H. Busch, *A Note on the Number of Hamiltonian Paths in Strong Tournaments*, Electronic Journal of Combinatorics 13 (2006), N3. Corollary 1 on PDF page 4 gives the exact unrestricted strong minimum, which is 75 at order nine. Source: https://www.combinatorics.org/ojs/index.php/eljc/article/download/v13i1n3/pdf . Direct PDF retrieval succeeded.

3. Maria Abi Aad, *About the number of directed Hamiltonian paths in special tournaments*. The author-deposited preprint record, posted November 25, 2025, identifies the author: https://papers.ssrn.com/sol3/papers.cfm?abstract_id=5804403 . The publisher search extraction identifies Discrete Applied Mathematics 394, December 15, 2026, pages 220--228, DOI https://doi.org/10.1016/j.dam.2026.07.017 . Its described scope is special reversed-arc families, not all-odd spectrum completeness. Direct publisher and SSRN page opening returned 403/redirect failures; primary-source indexed text was available. Malformed mathematical extraction is not used as a proof input, and no priority claim is approved.

4. Amine El Sahili and Maria Abi Aad, *About the number of directed paths in tournaments*, Discrete Applied Mathematics 277 (2020), pages 115--126, DOI https://doi.org/10.1016/j.dam.2019.09.005 . Primary publisher indexed text confirms the earlier special-family work. No theorem here depends on reconstructing malformed extracted formulas.

## 6. Remaining boundary and audit handoff

The accepted claims give exact representations and sharp obstructions for restricted construction grammars. They do not prove spectrum completeness at k=2, do not resolve either infinite prime lane in THM-4094, and do not settle the global conjecture. A later intercyclic reduced-word classification needs a separate audit and manifest entry.

The bounded startup was read-only. The shared repository remained clean. The optional processor status read was unavailable because no `.machine-id` exists; no registration was invented. This coordination-only limitation does not affect any mathematical or replay result. Root/publisher owns all repository publication and integration of the historical-correction lineage.
