# ArithmeticSeams

Small standalone Lean 4.30.0 package using `Std` only. The public root
`ArithmeticSeams.lean` imports every proved statement, and `AxiomAudit.lean`
imports that root before inspecting every theorem.

The operations are explicit proposed models of shifted arithmetic:

```text
shift a = a+1
augAdd a b = a+b+1
augMul a b = a*b+a+b
augPow a 0 = 0
augPow a (n+1) = augMul a (augPow a n)
```

The package proves the two shift identities, commutativity, associativity,
both distributive laws, multiplicative identity0, diagonal formulas, and
`augPow a n + 1 = (a+1)^n`. It also proves the exponent-addition law and
that augmented addition has **no additive identity on Nat**. The shift's
image is precisely the positive naturals. Thus this is not a claimed
semiring on all Nat with additive zero0: an additive identity would need
the additional label-1, outside this package's domain.

The 20 theorem statements include an exact finite control
`augMul 2 3 = 11`, `augPow 1 6 = 63`, `augPow 1 16 = 65535`.
Their proofs are kernel checked. No Collatz, prime, quadratic-cycle,
Hamiltonian-path, or sphere theorem is asserted by this package.

From this directory run:

```text
python verify.py
python -O verify.py
```

The verifier rejects unexpected Lean source files, external Lake packages,
missing root imports, omitted theorem audits, and proof escapes. It clears
the prior status before rebuilding, fails on every nonzero process exit,
and records source hashes plus the actual axiom report in
`verification.json` and `verification.log`. Normal and optimized runs pass.
The reported dependencies are at most Lean's `propext` and `Quot.sound`;
the finite control has no axioms. No custom axioms or unproved declarations
are used. The source mathematics is in
`05-knowledge/results/arithmetic_seams_20260921_operations.md` at the repo root.
