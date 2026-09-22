# Collatz blueprint: a checked premise audit and descent interface

**Status: PROVED in Lean 4.30.0 in the scopes below. Collatz convergence remains OPEN.**

This independent package uses Lean's bundled `Std` library, with no external
packages or mathlib download. Its explicit public root,
[`CollatzBlueprintAudit.lean`](CollatzBlueprintAudit.lean), imports
[`CollatzBlueprintAudit/Basic.lean`](CollatzBlueprintAudit/Basic.lean) and
[`CollatzBlueprintAudit/ClockAudit.lean`](CollatzBlueprintAudit/ClockAudit.lean).

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

## Actual step counts and the follow-up certificate

The follow-up attachment defines the **standard unaccelerated** map while
calling it Syracuse. Its iteration time therefore counts both odd and even
steps. The new module computes `evenSteps n t`, `oddSteps n t`, and
`standardCarry n t` from the actual standard-map orbit. It proves

```text
evenSteps n t + oddSteps n t = t,
2^(evenSteps n t) * iterate collatz t n
  = 3^(oddSteps n t) * n + standardCarry n t.
```

At an even step the carry is unchanged. At an odd step, with prior halving
count K, it becomes `3*B+2^K`. This is a derived invariant, not an assumed
certificate identity. A root-imported theorem makes convergence equivalent
to the positive-time margin for these actual counters at every `n>1`.
The global margin remains unproved.

The attachment's separate local certificate leaves K and B freely chosen.
It is satisfiable at `n=2` with total time `t=1`, `K=3`, and `B=2`, even
though the actual counters are `K=1`, `oddSteps=0`, `B=0`. With the true K,
using the incorrect exponent `3^t` would demand `2=6+B` and is impossible.
The package checks both facts, so it does not mistake unguarded existential
metadata for a true orbit tally.

One certificate at a single start is only local. The positive-preserving map
`localTrap n = if n=2 then 1 else n` has the same local certificate at two,
but three is fixed and prevents global descent. This refutes the generic
local-to-global inference; it does not refute the Collatz-specific equivalence
from the attachment, whose forward direction would require proving the open
conjecture. The attachment's separate all-natural-number convergence claim
is directly false at zero, as formally checked here.

Exact paired controls expose another lost coordinate:

| start | standard time | endpoint | even steps | odd steps | carry |
|---|---:|---:|---:|---:|---:|
| 23 | 10 | 5 | 7 | 3 | 19 |
| 95 | 6 | 323 | 3 | 3 | 19 |

Both starts are `5 mod 9` and have the same carry and odd-step count. The
halving totals differ, and one segment descends while the other grows.

## Reproduction and trust boundary

From this package directory, run:

```text
python verify.py
```

The verifier runs `lake clean`, `lake build`, and
`lake env lean AxiomAudit.lean`, and saves the actual output to
[`verification.log`](verification.log). The axiom audit imports only the public
root, so it also checks root-import reach. All public theorems in the audited
library modules must appear in the audit.

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
