# Session215: fixed cyclic-cap composition

**Scope:** rigorous specialization/application of inherited tournament path-cover composition. No Chan–Pak completeness, coverage improvement, or generic new tomography claim.

## Results

- Fixed nonempty cyclic caps A,C read only F_B(1..min(|A|,|C|)+1).
- An ordinal word B=D_1=>...=>D_r has exact finite-mode output `sum_j beta_j product_i Q_Di(j)`. Singleton blocks recover the companion multi-exit integer modes; arbitrary tournament seeds replace these by Q_D(j).
- An explicit five-vertex seed has H=11,Q(2)=100. It yields strong tournaments of order 5r+2 and H=100^r+11^r for every r>=1.
- Symmetric caps C3[TT_m,T,TT_m] give triangular profile probes with diagonal (m!)². This is an **alternative actual-context family**, not the first physical tomography theorem: Session176 already proves actual probes and contextual equivalence; the later proof-principles packet already proves parity-compressed probe recovery. THM-3166 supplies earlier negative-evaluation recovery.
- Fixed caps are essential to the finite-mode statement. Growing all three blocks needs the inherited THM-3202 moving jet.

## Files and replay

`PROOF.md`: definitions, proofs, exact inheritance, hostile controls, limitations.
`probe.py`: subset path and set-partition profile helpers.
`verify.py`: deterministic verifier with independent permutation and forward Held–Karp controls.
`receipt.json`, `receipt.optimized.json`: matching outputs.
`SOURCES.json`: pinned inherited-source Git blob IDs and content SHA256.
`SHA256SUMS`: final artifact hashes (regenerated after audit corrections).

From this directory:

    python3 verify.py > replay.json
    python3 -O verify.py > replay.optimized.json
    cmp replay.json receipt.json
    cmp replay.json replay.optimized.json

The verifier checks 160 physical cyclic compositions, 76 independent profile comparisons, 76 triangular recoveries, and 79 compressed recoveries. All 75 labelled tournaments of orders 1..4 are included, with explicit larger cases through 13 output vertices. It also checks the independent multi-exit bridge beta=(87,-189,72,36). Checks remain active under -O.
