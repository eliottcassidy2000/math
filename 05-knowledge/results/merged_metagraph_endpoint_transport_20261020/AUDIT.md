# Independent audit: endpoint transport and complement leakage

Status: PASS for the general proofs in their stated scopes, with independent exact controls at n=3,4,5,6. No sufficiency or unique-arrangement statement is endorsed.

## Proof audit

The endpoint-current Laurent polynomial and folded score-profile identities have the correct Hamiltonian-path/automorphism factor, converse factor, endpoint orientation, and loop normalization. The total triangle-current is constant on a merged fibre exactly when both extremal SCCs are regular. This concerns the total current, not merely its blue restriction.

The regular-to-almost-regular shell constraint and its reverse high-to-low endpoint formula are correct. The almost-regular shell is strong for n>=5. The entire triangle-current-rigid node set is independent for n>=4, giving the stated colour-by-colour half-mass bound. The n=3 exception is essential.

The complement leakage covariance, exact rank q-c-b, pseudoinverse colour decomposition, and equality conditions are valid. These are generic compression identities unless coupled to the endpoint inventories. The explicit strong-four witness proves triangle leakage at least 12/5 for every n>=5; the separate n=4 value is 14/5.

## Independent exact checks

`independent_transport.py` uses a separate all-permutations classifier in `independent_classifier.py`, independently constructs the descending-path cube from an ascending implementation, and independently assigns first-occurrence node labels. It imports no primary-implementation module or atlas. It enumerates all Hamiltonian paths and automorphisms of one representative per merged node.

Across all 49 merged nodes at n=3,4,5,6, every colour-resolved triangle-current row and every folded-score-profile row agrees. Independent reachability-based SCC decomposition verifies the regular-extremal-component criterion. Exact rational matrices verify the colour law of total variance.

The n=5 hostile switch preserves the claimed global current and score-pair data, colours, masses, and loop data, but changes the two prescribed blue row polynomials. Its support is not preserved: two existing contacts disappear and two other existing contacts become double. No previously forbidden contact is introduced. This caveat is correctly explicit in the paper.

The hostile is a necessary-constraint separation, not a proof that the remaining constraints characterize realizability.

## Reproduction

Keep both audit Python files together; Python 3 and SymPy suffice.

    python3 independent_transport.py

The generated `independent_transport.json` records all nodewise HP/Aut and colour-current inventories. The program runs without primary census files; if the matching local primary census JSON files are present, it additionally checks direct equality against them. No broad Walsh-distance audit output is needed for this transport receipt.

Packaging note: the independently audited classifier functions were extracted unchanged into `independent_classifier.py`; unrelated audit entry points and unused NumPy imports were omitted. The portable programs were rerun together. Original source fingerprints are retained in `DEPENDENCIES.json`.
