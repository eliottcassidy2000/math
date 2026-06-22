---
id: HYP-2874
status: SUPPORTED proof-order / closure lemma
source: codex-2026-06-22-S97
tags: [tournaments, h-spectrum, productive-closure, boundary-functions, strong-components, ocf, forbidden-H]
related:
  - THM-079
  - THM-115
  - THM-200
  - THM-343
  - HYP-2872
  - HYP-+2873
---

# HYP-2874: the productive H-closure is an admissible boundary closure, not a wild quotient closure

The useful closure for the permanent `H` gaps is target-specific and
`H`-preserving:

1. Condense a tournament into strong components.  The OCF/strong-component
   ledger gives

   ```text
   H(T) = product_i H(C_i)
   ```

   up to unit factors from acyclic/singleton components.

2. If a target value remains as one strong component, use the strong-core
   lower-bound boundary.  For `H=21`, THM-115 records the Moon pancyclic
   odd-cycle bound

   ```text
   alpha_1 >= (n-2) + sum_{odd L, 5<=L<=n} ceil(n/L),
   H >= 1 + 2 alpha_1.
   ```

   This gives `H>=25` for every strong core with `n>=9`.

This is the productive closure: every candidate route for a fixed low value
is either a factor route through smaller strong values, or a single strong-core
route closed by the finite base plus Moon boundary.

## Exact low-value ledger

The S97 script checks the factor routes for `{7,21}`.

For `H=7`, the only non-unit route is a single strong component with `H=7`.
THM-200/THM-343 close this route.

For `H=21`, there are two non-unit routes:

- `(3,7)`: closed because the factor `7` is forbidden by THM-200/THM-343.
- `(21)`: closed by the exhaustive base `n<=8` in THM-079 Part G and the
  `n>=9` Moon boundary in THM-115.

Thus `{7,21}` are closed by the admissible product/core boundary, without
using even-graph minor moves or any non-preserving quotient.

## Boundary-function reading

The analogy from boundary theory is not cosmetic.  Classical boundary-value
theorems distinguish admissible approaches, where the limit object is stable,
from arbitrary curvilinear approaches, where almost any boundary behavior can
be forced by choosing enough freedom in the approach.

In this `H` setting:

- strong-component condensation is the admissible approach;
- the OCF packet value is the boundary value;
- Moon/Busch-style strong-core lower bounds are the regular boundary theorem;
- even-graph smoothing and GF(2) contraction are wild curvilinear approaches:
  they preserve parity/cycle-space data but not `H`.

This explains why HYP-2872's even-graph minor quotient is useful as an address
but cannot carry the proof.  The preservation theorem is missing, and the audit
shows it is false for those moves.

The incoming S33 minor-order phrasing should be read through this lens:
`contract` means strong-component condensation with an `H` factor ledger, not
arbitrary even-graph contraction.  Any suppression/removal move must prove what
it preserves before it is allowed into the productive closure.

## Guardrail: not the false `7Z` ideal

Incoming S33 corrected the tempting but false ideal picture.  Strong
tournaments at `n=7` already realize values divisible by `7`, including

```text
35, 49, 133, 147, 175.
```

So the permanent obstruction is not an infinite ideal `7Z`.  The right
Kuratowski/Wagner analogy is a finite low-boundary forbidden set, with
strong-core re-entry later in the spectrum.

## Artifacts

- `04-computation/h_productive_boundary_closure_codex_s97.py`
- `05-knowledge/results/h_productive_boundary_closure_codex_s97.out`
- `07-reflections/h-spectrum-admissible-boundary-closure-codex-s97.md`
