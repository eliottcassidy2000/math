# CatalanEllipticAudit: finite kernel certificates

**VERIFIED scope:** twelve closed arithmetic statements in Lean4.30.0,
using only Std and ordinary `decide`. Every audited theorem is axiom-free.
This package does not formalize an elliptic-curve classification, the
bounded-period census, Catalan's theorem, or Collatz convergence.

From this directory run:

```text
python verify.py
```

The verifier requires the declared toolchain, inventories every public
theorem and Lean source, checks the absence of proof escapes and custom
axioms, runs `lake clean` then `lake build`, and imports the public root
from `AxiomAudit.lean`. Each subprocess must return zero. It records
the complete axiom audit and SHA256 hashes in `verification.json`, with
the command transcript in `verification.log`. A failed command cannot
produce a PASS record. There are no external package dependencies.

The constants `aa,bb,cc` are copied literally from the elliptic lane's
[exact experiment](../../experiments/catalan_elliptic_20260921_elliptic.py).
For positive a,b,c, define

```text
numerator=a(a+b)(a+c)+b(a+b)(b+c)+c(a+c)(b+c),
denominator=(a+b)(a+c)(b+c).
```

The rational fruit sum is `numerator/denominator`. The certificates prove
that the pasted triple fails `numerator=4*denominator`, while the repaired
triple `(aa,bb/10,cc/10)` satisfies it. Separate theorems prove that both
divisions by10 are exact natural-number divisions and that every pair sum
in both triples is positive. Thus no floor-division or zero-denominator
interpretation is hidden in the repaired witness. The universal statement
relating a rational fraction equation to its cleared form is elementary
written mathematics; these Lean declarations certify the concrete cleared
equations and denominator facts.

Further certificates verify `(10,49)` on `y^2=2x^3+4x^2+1`, the clock gap
`2^11-3^7=-139`, its ordered carry2363, the cancellation `2363=17*139`, and
the corresponding fixed-point equation. The seven cleared transitions of
the known negative Collatz cycle and the oddness of all successors are
also checked. These are witnesses, not a completeness assertion about
integer points, integer cycles, or higher periods.

Public entry point: `import CatalanEllipticAudit`. The proof-bearing module
is `CatalanEllipticAudit/Certificates.lean`; the root imports it directly.
