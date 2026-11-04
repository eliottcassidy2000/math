# Independent audit of the session179 crossfield packet

Verdict: **PASS after two wording clarifications**, now applied by the author. No mathematical counterexample found within the stated scopes. Collatz and odd-perfect-number conjectures remain open.

Audited: `../crossfield/REPORT.md`, `verify.py`, `valuation_dual.py`, and their JSON outputs. No tested source module was imported by the independent checker. The independent script uses only standard-library exact integers and fractions, including exact integer square roots for primality checks.

## Proof audit

- **Valuation bridge A1:** equality of all prime valuations is equivalent to the displayed equality of positive integer products. Multiplicativity of sigma applies to actual distinct-prime factorization blocks, while repeated same-prime blocks are explicitly a relaxation. This distinction is correctly preserved.
- **Separation A2:** the fixed dual has nonnegative pairing with all30 columns and negative pairing with e2, so it excludes nonnegative real combinations, hence nonnegative integer combinations. The prime and exponent ranges remain finite. The 23 listed coordinates and five zero-pairing columns were independently confirmed.
- **Max-product B1:** cycle inequalities are necessary by multiplication. For sufficiency, deleting a cycle of normalized product<=1 cannot decrease the path product; only finitely many simple paths need be considered. The resulting positive rational h obeys every edge inequality. The finite recurrence correctly detects stabilization within |V| updates.
- **Finite-core entry:** the weighted affine recurrence has strict descent outside B/(1-rho). For positive integer arithmetic states there are finitely many states (vertex,n) below any fixed starting V, so indefinite strict decrease outside the core is impossible. This proves entry or language exit, not universal language coverage or Collatz convergence. The author explicitly restricted the arithmetic domain to positive integers following audit feedback; unrestricted integers would make the core infinite.
- **Spectral versus worst path:** K^L has rank one with nonzero eigenvalue ((a0+a1)/2)^L; K has sole nonzero eigenvalue (a0+a1)/2. The author clarified this antecedent following audit feedback. The half-moment example has average <1 but an expanding all-odd loop, correctly separating quantifiers.
- **Parity-language boundary:** every finite parity word has arbitrarily large positive representatives. Thus a finite-prefix language covering all positive starts cannot exclude any word. This does not rule out height-dependent observers or liveness conditions, as the report correctly states.

## Independent finite checks

`independent.py` separately:

1. Reconstructs the width-five graph by four-bit integer shift arithmetic and verifies all supplied h inequalities and additive bounds.
2. Enumerates all44 directed simple cycles (one representative at minimum vertex) and verifies each exact slope product<1.
3. Enumerates sources modulo32 and computes actual five-step shortcut orbits, independently recovering all26 words and sharp intercept19/8 with slope27/32.
4. Recomputes every sigma(p^a) by brute summation over all divisors of p^a, independently verifies prime coordinate primality, rebuilds each rational ratio, and checks all30 dual scores.

Normal and optimized independent outputs byte-match (`independent.json`). Author scripts replay successfully; graph outputs match stored `results.json` and optimized output, valuation output matches stored `valuation_results.json` and optimized output.

Reproduction from any directory:

    python /workspace/shared/session179/algorithms_audit/independent.py

The checker uses a relative sibling lookup for the crossfield packet. These computations support the finite identities, not an unbounded prime/exponent claim or universal Collatz statement.
