# Independent almost-regular controls, Session30

Status: FINITE-EXACT for integer H, exhaustive small universes, and forest coefficients; floating-point VERIFIED for spectral/propagator evaluations. No asymptotic theorem is claimed. No parameters fitted. All scripts are independent implementations, not imports from producer scripts.

## Definitions and source boundary

Audited sources: `01-canon/theorems/THM-4757-linear-forest-loop-expansion-of-hamiltonian-paths.md`, `05-knowledge/hypotheses/HYP-9365-two-loop-law-for-hamiltonian-paths-in-regular-tournaments.md`, and `05-knowledge/results/sierpinski_audit2_20261022.md`.

S=A−Aᵀ, M=S/n, δ=M1; H counts directed Hamiltonian vertex orderings, κ=2^(n−1)H/n!, ρ=exp[(logdet(I−M²)−logdet(I+M²))/2], y=n(1−κ/ρ). Baseline y_reg is the regular two-loop expression from HYP9365. Tested near-regular extension is y_path=y_reg+nδᵀ(I+M²)^−1δ. This source extension is a known near miss: no second-order nonregular claim is assumed. For nonregular tournaments, ρ is a defined spectral quantity, not the complete excess-zero cluster sum.

Source guards: HYP9365's uniform remainder is open; universal regular limit e is false. HYP9367 later engine additions were not covered by audit2.

## Exact exhaustive controls

All labelled tournaments were enumerated at n=2,4,6 and retained iff every signed score is ±1. Counts: 2,24,2640. Full permutation orbits give respectively 1,1,5 isomorphism classes. Independent brute-force permutations validate H for every representative.

At n6, class H-values are 41,43,43,45,45; labelled orbit sizes 720,720,720,240,240, so automorphism sizes1,1,1,3,3. Two classes with H45 have distinct spectral y:1.00187583474646 and0.83694619443673. Thus equal H does not mean equal ρ or spectral residual. The complete bit strings are in small_classes.json.

## Frozen larger controls

Deterministic rotational deletion, fixed-seed circulant deletion, and triangle-reversal-chain deletion families at every even n8..20; Paley deletions at n6,10,18. Each is almost regular. The chain is not claimed to sample uniformly. See PREREGISTRATION.md, written before computation.

|family|n|exact H|y|y_path|n(y−y_path)|
|---|---:|---:|---:|---:|---:|
|rotation deletion|8|629|1.636465702|2.467733773|−6.650144565|
|rotation deletion|12|513669|1.951910681|2.595673243|−7.725150739|
|rotation deletion|16|1482122869|2.142916985|2.672170039|−8.468048858|
|rotation deletion|20|11166605274181|2.272576960|2.722311052|−8.994681846|
|chain deletion|12|525179|1.357496506|1.758080408|−4.807006817|
|chain deletion|20|11186410301151|1.432005841|1.690294446|−5.165772084|
|Paley deletion|6|45|0.836946194|1.569031608|−4.392512480|
|Paley deletion|10|15745|1.000227624|1.524856805|−5.246291811|
|Paley deletion|18|117266659317|1.212698899|1.502307143|−5.212948377|

The observed residual is of the scale expected from a missing second-order term. It is family-sensitive: subtracting a single empirically fitted constant/n would not remove both flat and rotational discrepancies. This is a diagnostic observation, not proof of a limit or bound.

Hostile transitive controls have H=1 and signed-score norm squared n(n²−1)/3, rather than n. Their predictor failure grows rapidly, e.g. n12 gives y≈11.99998 versus y_path≈83.41072. The almost-regular leaf-cost premise is absent. These refute unqualified extension to arbitrary tournaments, not an almost-regular asymptotic statement.

## Independent forest cross-check

forest_exact.py computes signed Hamiltonian sums on every odd subset by integer endpoint DP, then constructs forests by the unique component containing the least remaining vertex. No cluster formulas are used. The exact identity κ=Σ N_L/(n)_L is verified on every small isomorphism representative and all frozen larger/hostile/regular controls through n12. forest_controls.jsonl contains every N_L and the exact rational coefficients of log Φ.

Elementary degree-two identity (direct proof):

N₂=Σ_(distinct i,j,k) S_ij S_jk = 1ᵀS²1−tr S² = n(n−1)−||S1||².

For almost-regular n this becomes n(n−2), independent of isomorphism type. Reversing an odd-length component proves N_L=0 for odd L. These facts are checked exactly. This is an independent control, not a novelty claim against the repository.

## Reproduction

Compile `g++ -O3 -o count count.cpp`; run `OPENBLAS_NUM_THREADS=1 python experiment.py`; `python classify_small.py`; `python forest_exact.py`.

count.cpp uses unsigned128 integer endpoint subset DP, safely beyond n20!; Python permutations supply independent controls. No floating-point computation enters H or N_L. Spectral outputs use NumPy double precision, not interval certificates. Rounded values do not justify exact equality claims. The optional n22 Paley case in preregistration was not run; bounded resource scope is n≤20 for H, n≤12 for full forest controls.

Files: PREREGISTRATION.md, REPORT.md, experiment.py, count.cpp, classify_small.py, forest_exact.py, data.jsonl, summary.json, small_classes.json, forest_controls.jsonl, run.out, forest.out. Compiled count binary and __pycache__ excluded from publication.

## Exact isospectral obstruction

Already at n6, almost-regular bit strings 000110001000000 and 000110011000000 have H41 and H45, respectively, yet both have characteristic polynomial λ⁶+15λ⁴+15λ²+1. Thus no function of the skew spectrum alone can recover H, even within almost-regular tournaments of a fixed order. Here bits record S_ij=+1 for bit1, i<j in lexicographic order. Full permutation enumeration confirms both H-values. Exact Newton identities and independent SymPy determinants confirm the polynomial.

The shared ρ is exactly66637/27755. But the path-tree quantity differs:42372/27755 versus36/35. This displays information missing from eigenvalues: how the all-ones vector couples to the skew eigenspaces. Additional exact witnesses at n8,10,12 are in isospectral_witnesses.json. Their H values use exact endpoint DP, not all-permutation validation at those larger orders.

## Lollipop extension: finite held-out success

After the first frozen outcomes, the theory lane supplied L=−2diag(H−G)ᵀG1 and y_lollipop=y_path−nL. All original cases are explicitly retrospective. A new family was preregistered before its outcomes: a non-vertex-transitive regular7 base with H171, lexicographic product with C3, then deletion of one vertex in each base block (n20). Seven deletions comprise two orbit types, not seven independent isomorphism classes.

For the exceptional block, exacty1.792896927, pathprediction2.124692123, lollipopprediction1.872538829. Other blocks: exacty1.871523823, path2.215094209, lollipop1.965194921. Error shrinks from0.332/0.344 to0.080/0.094 with no fit. This supports the structural correction but does not validate a second-order formula or uniform error bound. The formula's proof and core audit are the theory lane's responsibility; this lane independently tests its numerical consequences.
