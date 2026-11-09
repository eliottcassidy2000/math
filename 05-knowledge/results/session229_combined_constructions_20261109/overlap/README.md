# Session229: exact overlap mechanisms and barriers

**PROVED within stated block templates; Chan–Pak completeness remains OPEN.**

- Every nonempty block kills nonadjacent equal-(H,Q,alpha) defect interactions. All scalar words have an exact nearest-neighbor quadratic formula, with efficient extrema and exact set-valued target compilation.
- Equal-U alpha words on repeated modules are classified by bounded integer carries. Coprime H,Q force each alpha digit; noncoprime carries have |c|<=gcd(H,Q)-1.
- Literal interacting H19 masks and a distinct H9/Q84 compensating carry provide positive physical witnesses, not only obstructions.
- Growing common-divisor cost is inherited from Session228 and specialized with a direct proof.
- The all-mask, all-length strong-flanked H19 template has only O(X^(3/4)) attainable scalar values <=X. The proof includes arbitrary adjacent interactions and changing sidecars, and extends to finite block banks under an explicit sufficient entropy condition.

See `PROOF.md` for complete statements, scope, mechanism, provenance, and next questions. No interval of unrestricted coverage is claimed. `verify.py` is standalone standard-library code. Run:

    python3 verify.py > replay.json
    cmp replay.json CHECK.json
    python3 -O verify.py > replay.optimized.json
    cmp replay.optimized.json CHECK.optimized.json

Both receipts are deterministic and identical. `probe.py` and `scout.py` are exploratory provenance, not dependencies or canonical checks; the publication packet need only include README, PROOF, verify, and the two receipts. They should not be cited for an exhaustive larger search.
