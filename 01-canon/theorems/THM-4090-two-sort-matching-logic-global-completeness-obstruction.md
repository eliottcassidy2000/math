---
id: THM-4090
title: "Two-sort matching-logic global-completeness obstruction"
status: >
  PROVED FINITE SORT-FLOW SHARPENING + CITED SOUNDNESS INPUT + INDEPENDENTLY
  RULE-BY-RULE AUDITED, relative to the many-sorted definedness-free
  fixpoint-free Figure-2 matching-logic semantics and calculus of
  arXiv:2608.13306v1. Two nonempty sorts and one unary symbol
  f:b->a already admit a satisfiable theory Gamma and a closed b-sorted
  consequence phi with Gamma models phi but Gamma does not derive phi. This
  sharpens the preprint's three-sort Proposition 30 to a proved two-sort
  failure. Together with the preprint's claimed one-sort theorem, this would
  make its sort-count boundary sharp. It refutes global completeness of that
  calculus, not axiomatizability or all calculi.
source: codex-padic-zeta-tournament-20260825
depends_on: []
related:
  - THM-4088-order-tournament-arithmetic-type-blindness-and-lrc-margin-density
external_source: https://arxiv.org/abs/2608.13306v1
---

# THM-4090 -- global-completeness failure already occurs with two sorts

**PROVED finite sort-flow sharpening + CITED soundness input + INDEPENDENTLY
RULE-BY-RULE AUDITED.**

The source
[“Completeness and incompleteness of basic matching logic,” arXiv:2608.13306v1](https://arxiv.org/abs/2608.13306v1)
claims global completeness for the one-sorted, definedness-free,
fixpoint-free system and gives a three-sort counterexample in Proposition 30.
The elementary sharpening below needs only two sorts and one unary symbol.
See the [local source audit](../../05-knowledge/reference/arxiv-2608-13306-matching-logic-source-audit.md)
for the theorem ledger and imported inputs.

The status of the preprint's own new results remains **PREPRINT CLAIM / UNDER
AUDIT**. This theorem uses its explicitly displayed semantics and Figure-2
rule typing as definitions, plus the earlier soundness theorem cited there,
and gives a separate finite sort-flow argument.

## 1. The two-sort theory

Let the carriers `M_a,M_b` be nonempty, let the signature consist of

```text
f:b -> a,                                                   (1)
```

and use the many-sorted matching-logic semantics of the source. Element
variables denote singleton subsets of their carrier, symbols extend
pointwise to subsets and propagate the empty set, and a pattern is satisfied
when its denotation is its entire result-sort carrier.

Put

```text
Gamma={forall x:b forall y:b. f(x and y)},                 (2)
phi=forall x:b. x.                                         (3)
```

The outer sort of the hypothesis in `(2)` is `a`; quantification does not
change the outer sort. The conclusion `(3)` has outer sort `b`. Both are
closed.

### Theorem 1.1

```text
Gamma is satisfiable,       Gamma models phi,       Gamma does not prove phi. (4)
```

Consequently the many-sorted Figure-2 calculus is not globally complete
already for two sorts.

## 2. Satisfiability and semantic consequence

Take both carriers to be singletons and interpret `f` totally. Then `(2)` is
total, so `Gamma` is satisfiable.

Conversely, suppose a model satisfies `(2)`. If `r,s in M_b` were distinct,
choose a valuation with `x=r,y=s`. Since element variables denote
singletons,

```text
rho(x and y)={r} intersect {s}=emptyset.                   (5)
```

Pointwise symbol extension gives `f(emptyset)=emptyset`, which cannot equal
the nonempty carrier `M_a`. This contradicts totality of `(2)`. Hence every
model of `Gamma` has singleton `M_b`.

For a singleton `M_b={w}`,

```text
rho(forall x:b. x)=intersection_(z in M_b){z}={w}=M_b.    (6)
```

Thus `Gamma models phi`. On a two-element `b`-carrier the intersection in
`(6)` is empty, so `phi` is not valid.

The source's longer conclusion

```text
forall x:b forall y:b. (x iff y)                           (7)
```

works as well. Formula `(3)` is the minimal one-variable version.

## 3. The sort-flow non-derivability proof

Say a sort feeds another if they are equal or a positive-arity symbol has an
argument of the first sort and result of the second, and take transitive
closure. In signature `(1)`,

```text
a => a,             b => b,             b => a,
a does not feed b.                                           (8)
```

Every inference rule of the source's many-sorted Figure-2 system has premises
whose sorts feed the sort of its conclusion:

| rule family | sort behavior |
|---|---|
| tautologies, quantifier axioms, propagation axioms, existence, singleton-variable axiom | premise-free at their result sort |
| modus ponens | all patterns have one outer sort |
| existential and derived universal generalization | outer-sort preserving |
| framing | moves from a symbol argument sort to its result sort, exactly one feed edge |

Induction on a finite derivation from `Gamma` now proves:

```text
every b-sorted derived line is derivable from the axioms alone. (9)
```

Indeed, the only theory-dependent initial line has sort `a`. It cannot be a
premise in the ancestry of a `b`-sorted conclusion by `(8)`, while all
premise-free lines remain available.

If `Gamma proved phi`, equation `(9)` would give `prove phi`. Soundness of
the displayed calculus would make `phi` valid, contradicting the two-element
model above. This proves the last assertion in `(4)`.

## 4. Conditional sharpness and exact scope

Independently of the preprint's positive theorem, this example proves failure
already at two sorts, so the one-sort hypothesis cannot be weakened merely to
“finitely many sorts.” If the preprint's claimed one-sort completeness theorem
is correct, the two statements together make the boundary sharp in **number
of sorts**. This is not a claim that every two-sort signature fails.

The same proof survives:

- arbitrary extra nullary symbols, because they create no feed edge or
  application context with a hole;
- extra positive-arity symbols only when they create no path `a=>b`.

It does not cover definedness, fixpoints, free set variables, or nominals.
Nominals are singleton-interpreted model-fixed constants, not ordinary
nullary symbols. It also does not prove that the semantic consequence relation
is unaxiomatizable. As the preprint stresses, many-sorted model consequence
translates to first-order consequence under effective hypotheses; a different
complete calculus may exist.

## 5. Relation to the repo frontiers

The result is a formalization guardrail, not a mathematical bridge to the
repo's headline problems.

- Encoding LRC(14), planar Keller pairs, or a p-adic irrationality statement
  as `Gamma models phi` would still require a faithful semantics-preserving
  encoding and a proof of the semantic consequence. Completeness alone
  constructs neither a lonely time nor an arithmetic contradiction.
- Natural recursive encodings may require fixpoints, where the preprint's
  Theorem 19 instead claims non-recursive-enumerability of validity.
- The argument-reachability relation in the paper may have loops, missing
  pairs, and both directions. It is not a tournament. The sort-flow relation
  in `(8)` is a preorder whose strict quotient is a DAG, also not a tournament.
- The exact shared mechanism is controlled forgetting: localization or sort
  projection preserves reachable consequences and loses an unreachable
  global carrier. That analogy transfers only after the source, target, map,
  and missing sidecar are specified.

Thus THM-4090 sharpens one external logical frontier while leaving LRC(14),
JC(2), the 22-value p-adic-zeta claim, and every named tournament inequality
open.
