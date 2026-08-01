# Comparator challenges

A **statement/proof separation** layer for this repo's Lean results, imported
from the pattern used by
[`github.com/openai/ten-proofs`](https://github.com/openai/ten-proofs)
(fetched and audited 2026-08-01, klein-S428).

## The problem this solves

A `sorry`-free, kernel-pure Lean proof certifies that *the stated proposition
follows from the axioms*. It certifies nothing about whether the stated
proposition **is the theorem**. Every remaining risk in a green formalization
concentrates in the statement, and specifically in the definitions the statement
is phrased with.

This repo has already been bitten by exactly that:

* `THM-3018` stated `FC(n)` **with a homogeneity hypothesis** that the actual
  factorial conjecture does not have — a scope error invisible to the kernel.
* `GMC2Main`'s own header records that the residual risk is *upstream* of the
  wiring, not in it.
* `MISTAKE-190` is the same shape: a conditional result recorded as
  unconditional.

## The pattern

Each challenge is a pair.

| file | role |
|---|---|
| `<Name>Challenge.lean` | **imports only Mathlib.** Defines everything it needs from scratch, states the theorem, and closes it with `sorry`. Short enough to read in full. |
| `<Name>Solution.lean` | imports the repo and discharges the challenge theorem, *same name, same type*. May be arbitrarily long. |
| `<Name>.json` | names the theorem(s) that must match and pins `permitted_axioms`. |

An auditor reads **only the challenge file**. If they accept it as a faithful
statement, and the solution builds with no `sorry` and no axiom outside the
permitted list, the result is certified — without reading the proof or any repo
definition.

## GMC(2)

`GMC2Challenge.lean` / `GMC2Solution.lean` / `GMC2.json`.

The risk in `GMC2Main.GMC2.gmc2` is the functional `E`. It is defined in
`GMC2Reduction.lean` as a sum over a polynomial's support against a weight `wt`,
and checking that this is the Gaussian expectation means checking repo code.

So the challenge **does not mention `GMC2.E`**. It quantifies over every
functional satisfying

```text
  additivity      E (P + Q) = E P + E Q
  homogeneity     E (C c * P) = c * E P
  Wick moments    E (X 0 ^ a * X 1 ^ b) = a! * [a = b]
```

Auditing the statement therefore reduces to accepting one textbook identity —
the moment law of a standard complex Gaussian, with `X 0 = Z` and `X 1 = conj Z`.
Since such a functional is unique, this is the same theorem, only independently
checkable.

The bridge lemma `GMC2Challenge.wick_eq_E` is where the work happens: it
*derives* the Wick weight `wt` from the moment law rather than assuming it, and
`wickFunctional_E` shows the hypothesis is non-vacuous by exhibiting `GMC2.E` as
an instance.

## Building

```sh
lake build ComparatorChallenges.GMC2Challenge   # the statement (has the sorry)
lake build ComparatorChallenges.GMC2Solution    # the discharge
```

`enable_nanoda` is `false` here because this repo does not vendor an external
kernel; the axiom check is the `#print axioms` at the foot of the solution.
Turning on a second, independent kernel would be a strict improvement and is the
obvious next step if anyone wants it.

## Adding a challenge

Good candidates are results whose *statement* carries risk: anything with a
hypothesis that could silently narrow scope, anything whose definitions are
repo-local, and anything a reflection describes as "kernel-pure but the risk is
upstream". A challenge is cheap — it is a definition and a `sorry` — and it
converts an unfalsifiable claim of correctness into a checkable one.
