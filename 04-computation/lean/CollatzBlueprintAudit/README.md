# Collatz blueprint: a checked premise audit and descent interface

**Status: PROVED in Lean 4.30.0 in the scopes below. Collatz convergence remains OPEN.**

This independent package uses Lean's bundled `Std` library, with no external
packages or mathlib download. Its explicit public root,
[`CollatzBlueprintAudit.lean`](CollatzBlueprintAudit.lean), imports
[`CollatzBlueprintAudit/Basic.lean`](CollatzBlueprintAudit/Basic.lean).

The inherited useful mechanism is finite-word Collatz iteration; the missing
coordinate is a stopping time for each individual starting number. The prior
squarefree-prefix results permit arbitrarily long growth and do not discharge
this requirement. The corrected statement here is a strong-induction reduction,
not a deduction from squarefree density or a monodromy analogy.

## Checked statements

- The blueprint premise that **every integer** has residue `1`, `3`, or `5`
  modulo `6` is impossible: substitute zero. A structure containing that field
  is consequently uninhabited. This checks a well-typed isolated version of the
  problematic field; the pasted template itself also contains syntax/type errors.
- Restricting the premise to odd integers yields a proved statement, including
  negative odd integers with Lean's Euclidean remainder convention.
- For a map `f : Nat → Nat` preserving positive inputs, every positive start
  eventually reaches `1` **if and only if** every start `n > 1` has some iterate
  at a **positive time** strictly below `n`. The time may depend on `n` and
  intermediate values need not decrease. Positive preservation keeps the lower
  endpoint in the domain of strong induction.
- The standard unaccelerated map
  `collatz n = if n % 2 = 0 then n / 2 else 3*n+1` preserves positivity, so the
  equivalence specializes to it. The theorem
  `collatz_reachesOne_of_strictDescent` retains eventual strict descent as an
  explicit hypothesis. Nothing here proves that hypothesis.
- If an actual orbit segment satisfies `2^K*y = 3^L*n+B`, then `y<n` is
  equivalent to the full inequality `3^L*n+B < 2^K*n`. The theorem
  `strictDescent_of_affineCertificates` accepts such identities and inequalities
  at a positive time for each start. Neither arbitrary word realization nor
  disappearance of the carry term `B` is assumed.
- Exact controls: `7` reaches `1` in `16` steps; zero remains zero and never
  reaches `1`. The constant-one map witnesses satisfiability of the general
  criterion, while the constant-zero map descends from every `n > 1` but does
  not always reach `1`, showing why positivity cannot be discarded.
- The odd step `19 → 29` increases the value. It takes two standard-map steps,
  and `(3*19+1) % 4 = 2` verifies that exactly one factor of two is removed.
- For the standard signed map on `Int`, `-1`, `-5`, and `-17` return to themselves
  after `2`, `5`, and `18` steps respectively. These are return equalities, not a
  classification of all signed cycles and not a minimal-period theorem.

`EventuallyReachesOne` means hitting `1` at some finite time, allowing time
zero for the start `1`. It does not require `1` to be fixed: the standard map
continues along `1 → 4 → 2 → 1`.

## Reproduction and trust boundary

From this package directory, run:

```text
python verify.py
```

The verifier runs `lake clean`, `lake build`, and
`lake env lean AxiomAudit.lean`, and saves the actual output to
[`verification.log`](verification.log). The axiom audit imports only the public
root, so it also checks root-import reach. All public theorems in `Basic.lean`
must appear in the audit.

[`verification.json`](verification.json) records the toolchain, source hashes,
theorem count, and actual dependencies. The only permitted foundational axioms
are Lean's standard `propext` and `Quot.sound`; several exact computations and
iteration lemmas use none. There are no custom axioms, admitted proofs, or
native evaluation trust shortcuts. `decide` checks the small finite examples
through kernel reduction. The verifier's own checks remain active under Python
optimization.

The finite controls and the proved equivalence do not prove Collatz, exclude
unknown positive cycles, classify every negative cycle, or establish any claim
about odd perfect numbers.
