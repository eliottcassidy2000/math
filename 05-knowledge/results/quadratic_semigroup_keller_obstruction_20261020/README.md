# All quadratic-semigroup Keller sources

20 October 2026. **PROVED; independently reviewed PASS.**

Let P,Q belong to C[x,y], with P of positive y-degree and nonzero constant leading coefficient, and J(P,Q) in C*. Then the degree semigroup of C(x)[P,Q] cannot be <2,m> for any odd m>=3. Equivalently, a proper primitive generic source algebra cannot contain a quadratic.

The [full proof](QUADRATIC_SEMIGROUP_KELLER_OBSTRUCTION.md) combines constant-leading target reduction, cited characteristic-degree bounds, an elementary filtration induction, polynomial descent of the quadratic center, and the cited one-line injectivity theorem. It does not assume m prime, saturation, unique embedding, or constant special-fiber behavior.

The key safeguards are retained: primitive m and the later-entry bound eliminate two-even terminal degrees; the one-line theorem applies only after an actual polynomial source shear. The two-even primitive hostile (t^4,t^6+t) shows why parity of initial degrees alone is insufficient.

See [sources](PRIMARY_SOURCES.md), [dependency boundary](DEPENDENCIES.md), and [review](review/AUDIT.md). The all-m conclusion is mathematical proof, not a finite census. The 111 exact symbolic controls pass normally and under Python -O with identical outputs.

From this packet directory, reproduce with:

    python3 verify_quadratic_semigroup.py > replay.json
    python3 -O verify_quadratic_semigroup.py > replay-optimized.json
    cmp replay.json results.json
    cmp replay-optimized.json results-optimized.json

The planar Jacobian conjecture remains open. No historical novelty claim is made.
