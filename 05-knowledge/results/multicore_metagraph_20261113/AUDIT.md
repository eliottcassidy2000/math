# Scope and reproduction

The mathematical proofs are in GROWTH, INCIDENCE, and FAMILIES. No theorem
number is reserved, no Lean verification is claimed, and the conjectural
extension of the six inflation exclusions is explicitly OPEN.

## Exact universes and independent paths

|Producer|Universe|Independent control or hostile|
|---|---|---|
|rooted_growth|Every end-rooted class through8 and every attachment to it|Full endpoint DP; all path permutations through5; every tiling through6; ordinary orbit weights|
|cycle_index|Odd cycle-type partitions through30|Rational Burnside integrality; independent agreement with generated rooted and ordinary classes through8|
|growth_audit|Every rooted deletion card through child7; every presentation lift profile through parent6|Different incidence formula G2; first ordinary-node nonlumpability witness at4|
|component_growth|Every end state through7, every attachment, all component-word masses through8|Independent reachability SCC decomposition and reciprocal-series weights|
|blue_incidence|All tilings and complement matches through7|Literal reversal/converse test; loop half-degrees; black multiplicity parity; ordinal and strong cyclic boundary hostiles|
|family_search|233 complete rooted cores; all fixed binary contexts of length≤12; both middle bits|Complete inherited ranker; no random sampling; coefficient and interface keys retained separately|
|family_audit|All28,776 retained recipes; all265 edges at three clocks|Exact matrix annihilators and positivity; literal graph maps;16 independent full-count/readiness checks|
|next_target|All3,524,577 unique-start strong parents through18|Per-target published order bound; graph-count checks of the resulting source witnesses|
|clock_probe|Every vertex of both ready881 strong sources|Independent literal TT2 counts; sole undershoot TT3 check; all-order source completeness proved separately|
|divisor_entry|All unique-start cores through15,37 specified target pairs, every compatible empty-prefix right word|Target divisibility determines the entire word;147 candidates, not a sampled word list|
|inflation_descent|All16 strong ready sources for the first six405-ray values; every potentially successful module-size vector|Different all-clone DP; positive inflation budget;40 exact candidates; singleton-SCC cost criterion over every ordinary class through7|

The core census total in the next-target row is the sum of the complete
inherited Fibonacci rank sizes at orders3,...,18. `next_target.json.gz`
retains those sizes and the order bound separately. None of the family
entries is counted as an unmarked isomorphism class merely because its
coefficient pair is different.

## Arithmetic and implementation safeguards

Most calculations use Python integers. The accelerated clone counter uses
signed64-bit integers only on old orders≤19, and candidate Hamilton counts
on orders≤20. An ordered two-bin cover of an n-vertex tournament is bounded
by(n+1)n!, so at n≤19 it is at most20!<2^63. All summands are nonnegative.
The candidate Hamiltonian table is bounded by20! as well. Larger inputs
are rejected by the experiment's explicit order guard.

The first SCC helper accidentally used a mutable list inside that same
list's sorting key. Python exposes the list as empty during its sort, so
components were left in label order. The per-state transition test rejected
it immediately. The repaired key reads an immutable snapshot; all2961
profiles then pass. This was a local implementation error, not a refutation
of the suffix-merger theorem. Its mechanism is recorded in MISTAKES.

The readiness type is E_y(2)−H. Its equality with deletion minus start count
is used only with the required indegree-one guard. The selected-cut defect
H−W(O) remains a different observable. Intermediate macro surgery cards
are not assumed ready. Colour/partner tests retain the same-vertex alignment
and never substitute the graph of single-arc reversals.

## Reproduction and source integrity

From the repository root:

    python -B 05-knowledge/results/multicore_metagraph_20261113/reproduce.py

The audited environment uses Python 3.13.14, NumPy 2.3.2 and Numba 0.64.0.
NumPy and Numba are required by the accelerated inflation checker;
Matplotlib is needed only to regenerate the separate figure. The receipt
checks do not depend on Matplotlib or on regeneration of that illustration.

The11 producer pairs are rerun normally and with Python optimization;
`reproduction.json` checks byte-identical deterministic receipts. The
optimization-only option compares a fresh optimized run with existing normal
receipts. This is a robustness check, not a substitute for the independent
mathematical engines listed above.

`SOURCE_PINS.json` records imported repository code and the exact inherited
proof files used here. `MANIFEST.json` hashes the completed packet, excluding
itself. The figure was rendered and visually inspected; it depicts the exact
order-four multigraph and the generated growth-state counts.

Primary literature was checked directly: Busch, *A Note on the Number of
Hamiltonian Paths in Strong Tournaments*, EJC13(2006), N3, Theorem1 and
Corollary1; and Chan–Pak, *Computational Complexity of Counting Coincidences*,
Conjecture6.2. A000568's published b-file supplies the ordinary census control.
The inherited complete ready grammar, orbit dictionary, and incidence
pullback retain their original provenance and are not claimed as new here.
