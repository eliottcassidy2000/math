# Session133: exact core floor and finite kernel budget

**PROVED (audited algebra):** [THEORY.md](THEORY.md) refines the inherited
singleton-SCC bound to an exact maximum core size under `H<=M`. It also
identifies the exact minimum H at each feasible core order, classifies the
minimizers for N>=6, and counts their labelled/unlabelled classes.
Corollary C adds supplied-tournament threshold recognition in O(n^2) graph
operations plus O(B^(3 log_5 2) log^3(B+2)) bit operations. This is not
prescribed-value existence and not polynomial in log B.

```
K(M)=max_{b=0,1,2;5^b<=M} [4b+3 floor(log_3(M/5^b))].
```

K(7)=4; K(21)=7. The extremizer for M=21 has **H=15**, not21.
The all-order exclusions of7 and21 and finite spectrum decidability are
inherited, explicitly cited; unrestricted completeness remains **OPEN**.

Files:

* `THEORY.md`: proofs, primary-source import, lineage, scope and hostiles.
* `check_core.py`: dependency-free integer verifier (Python3).
* `normal.json`, `optimized.json`: identical deterministic receipts.
* `SOURCE_PINS.json`: exact repository sources at pinned main commit.
* `SHA256SUMS`: payload integrity.

Reproduce:

```
python3 check_core.py --out normal.json
python3 -O check_core.py --out optimized.json
cmp normal.json optimized.json
sha256sum -c SHA256SUMS
```

Finite checks: 33,868 labelled tournaments through order6 (including empty),
2,131,019 induced-subset H-equality tests, literal permutation agreement
through order5, independent component-partition DP through core order400,
resource-budget comparison through M=100000, and constructed extremizers.
The all-order proofs do not extrapolate those tests.

Transitive-module hostile: blowing up one vertex of C3 to a transitive
module of size t gives H=2^t+1, but contraction gives3. Only adjacent
vertices of a transitive module are directed twins relative to every other
vertex; arbitrary members agree only outside the module. An initial test
incorrectly demanded arbitrary-pair twin equality; that test failed at t=3,
was corrected, and both complete runs were repeated. No theorem depended
on that failed test.
