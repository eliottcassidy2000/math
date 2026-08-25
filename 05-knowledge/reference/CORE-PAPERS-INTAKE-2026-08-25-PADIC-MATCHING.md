# August 25, 2026 source intake: p-adic zeta and matching logic

**Status:** source pin and import guardrail. External theorem headlines are
`PREPRINT CLAIM / UNDER AUDIT`; exact in-repo consumers retain their own
narrower statuses.

## Long -- *Hybrid arithmetic holonomy and twenty-two individual p-adic zeta values*

- **Primary / freshness:** [GitHub research draft at commit `b46a177`](https://github.com/octonion/p-adic-zeta-irrationality/tree/b46a1770901551961710e155d775aae7c5ea39e7),
  dated 2026-08-24. **EXTERNAL PREPRINT CLAIM / UNDER SPECIALIST AUDIT.** No
  separately indexed paper was located, the repository has only three
  same-day commits, and its README explicitly requests independent review.
- **Claimed result:** individual irrationality of `zeta_2(s)` for odd
  `3<=s<=29`, `zeta_3(s)` for odd `3<=s<=11`, and `zeta_5(3)`, `zeta_5(5)`,
  `zeta_7(3)`, using a hybrid Hasse/de Rham-torsor arithmetic cost, modular
  Jensen energy, and a Bost slope inequality.
- **Verified import:** the standard-library rational interval certificate
  replays all 22 stated terminal margins exactly; its smallest lower endpoint
  is positive in the `(p,s)=(5,5)` cell. This is **FINITE-EXACT** for the
  implemented formulas only. THM-4089 proves an independent elementary
  optimizer, four next-weight global no-go results, and the separate
  `p=13,s=3` all-radius no-go; THM-4091 proves the coordinate-change
  denominator boundary.
- **Does not prove in this repository:** any of the 22 irrationality claims
  until the global BGG/descent, Hasse-source, product-digit Cartier,
  fixed-bundle CDT, continuation-radius, and Bost-filtration interfaces receive
  specialist audit. The manifest binds only the script/output, not the TeX or
  proof. Python `assert` guards disappear under `-O`, and line-ending conversion
  changes the advertised raw hashes on a default Windows checkout.

## Chen--Rosu -- *Completeness and incompleteness of basic matching logic*

- **Primary / freshness:** [arXiv:2608.13306v1](https://arxiv.org/abs/2608.13306v1),
  submitted 2026-08-13. **CITED / UNREFEREED ARXIV V1.** This is the paper at
  the user-supplied arXiv link; it is unrelated bibliographically and
  mathematically to Long's p-adic-zeta draft.
- **Imported role:** proves global completeness for arbitrary theories over
  one-sorted finitary signatures in basic, definedness-free, fixpoint-free
  matching logic; gives a many-sorted incompleteness witness; and proves that
  validity ceases to be recursively enumerable after a least-fixpoint
  extension. Its positive mechanism is semantic localization to a
  backward-closed core followed by a two-sheet countermodel construction.
- **Repo consumer:** THM-4090 independently sharpens the many-sorted witness to
  exactly two sorts and one unary symbol. One sort is minimal only in the exact
  basic language/system scope above; nominals, set variables, definedness, and
  fixpoints are outside that statement.
- **Does not prove:** any p-adic localization, arithmetic double-cover,
  irrationality, interval-certificate, or effective proof-search statement
  about the finite certificate in Long's draft. The only legitimate bridge is
  a representation-loss audit that retains sort/weight/lattice sidecars.

## Consumers and reproduction

- THM-4089: exact optimizer, four next-weight failures, `p=11` chamber wall,
  and formal `p=13,s=3` one-power margin no-go.
- THM-4090: exact two-sort global-completeness obstruction.
- THM-4091: exact LCM/coefficientwise-depth coordinate-change boundary.
- THM-4093: determinant/Bowen--Lanford diagonal gauge and exact p-adic
  tournament-zeta tangent boundary; unrelated to Kubota--Leopoldt zeta.
- Full audit reflection:
  [`padic-holonomy-matching-logic-and-coordinate-depth-frontiers-codex-20260825.md`](../../07-reflections/padic-holonomy-matching-logic-and-coordinate-depth-frontiers-codex-20260825.md).
