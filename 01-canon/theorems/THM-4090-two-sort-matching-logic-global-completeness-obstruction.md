---
id: THM-4090
title: "Two-sort matching-logic global-completeness obstruction"
status: >
  PROVED + CITED MINIMALITY INPUT + FINITE-EXACT + INDEPENDENTLY AUDITED.
  The standard Hilbert system for basic many-sorted, definedness-free,
  fixpoint-free matching logic is already globally incomplete for a
  satisfiable theory over two sorts and one unary symbol. Two sorts are
  minimal existentially relative to the cited preprint-v1 one-sort
  completeness theorem. This is a calculus obstruction, not
  nonaxiomatizability.
source: >
  codex-padic-zeta-tournament-20260825; Chen--Rosu,
  arXiv:2608.13306v1, Corollary 15 and many-sorted system of Figure 2
depends_on: []
related:
  - arXiv:2608.13306v1
script: 04-computation/two_sort_matching_logic_obstruction_thm4090.py
output: 05-knowledge/results/two_sort_matching_logic_obstruction_thm4090.out
independent_audit: .scratch/matching_two_sort_referee_20260825/REPORT.md
script_sha256: ff64ea8dadc2fd2d4718ce303e93a87c22e95ad6db2a86cacc3a907a48bca758
output_sha256: 3ac734da6800af7278ad9e13c85e848773548d8869d5680b694f09cd9759f6b3
hash_basis: raw LF bytes
---

# THM-4090 -- two sorts already obstruct global completeness

**PROVED + CITED MINIMALITY INPUT + FINITE-EXACT + INDEPENDENTLY AUDITED.**

## 1. Exact setting and statement

Use the standard Hilbert system displayed in Chen--Rosu,
[*Completeness and incompleteness of basic matching
logic*](https://arxiv.org/abs/2608.13306v1), with nonempty carriers and
set-valued symbols.  Work in its basic many-sorted fragment: no definedness,
fixpoints, nominals, or set variables.  Take two sorts `b,a`, one symbol

```text
f : b -> a,
```

and the closed patterns

```text
Gamma = { forall x:b forall y:b. f(x and y) },
phi   =   forall x:b. x.                                (1)
```

Then `Gamma` is satisfiable and

```text
Gamma |= phi,              but              Gamma !|- phi. (2)
```

Consequently the displayed standard Hilbert system is not globally complete
in general with two sorts.  Two is the least possible number of sorts for
such an example relative to the paper's **CITED PREPRINT-v1** Corollary 15,
which claims global completeness in the same basic one-sorted fragment.  The
paper's new positive results have not been independently formalized here; see
the [bounded source audit](../../05-knowledge/reference/arxiv-2608-13306-matching-logic-source-audit.md).

The minimality statement is existential.  It does not say that every
two-sort signature or theory is incomplete.

## 2. Semantic mechanism

Element variables denote singletons, conjunction denotes intersection, and
a positive-arity symbol propagates an empty argument to the empty set.  If
the `b`-carrier contains distinct `r,s`, the valuation `x=r,y=s` makes

```text
[[x and y]] = empty,        [[f(x and y)]] = empty.      (3)
```

Since the `a`-carrier is nonempty, the universal intersection in `Gamma`
cannot then be total.  Every `Gamma`-model therefore has `|M_b|=1`.

Conversely, if `M_b={r}`, then `Gamma` is total exactly when
`f_M(r)=M_a`.  This gives a model for every nonempty `M_a`, so the entailment
is not vacuous.  Finally

```text
[[forall x:b. x]] = intersection_(r in M_b) {r},        (4)
```

which is total exactly when `|M_b|=1`.  Equations (3)--(4) prove
`Gamma |= phi`.

## 3. Sort-flow proof of nonderivability

Orient a feed edge from each input sort of a positive-arity symbol to its
output sort.  Here the only nontrivial feed is

```text
b -> a;                         there is no path a -> b. (5)
```

Induct on a putative finite derivation and strengthen the assertion to:

> Every derived line of sort `b` has an empty-theory derivation.

This is a rule-by-rule invariant of the displayed calculus.

- The sole hypothesis has sort `a`, so it cannot be such a line.
- Propositional, quantifier, propagation, Existence, and
  singleton-variable axiom instances already have empty-theory derivations.
- Modus ponens and existential generalization preserve the pattern sort, so
  their premises can be replaced by the empty-theory derivations supplied by
  induction.  The derived universal-generalization rule is handled the same
  way.
- Framing is the only primitive inference that may change sort.  It follows
  an input-to-output feed edge; for `f:b->a` it can conclude at `a`, never at
  `b` from an `a` premise.

Thus `Gamma |- phi` would imply `|- phi`.  But in a model with a two-element
`b`-carrier, (4) is empty, so `phi` is not valid.  Soundness of the displayed
system gives the contradiction and proves the second half of (2).

This is the two-sort analogue of the paper's sort-flow induction around
Corollary 35; it is not literally an instance of that corollary's fixed
three-sort signature.

## 4. Equality, scope, and tournament firewall

The one/two-sort boundary is sharp only for this exact language and calculus.
A genuinely global or cross-sort rule can break the induction, and the proof
does not rule out a different complete recursively enumerable calculus.  In
fact the paper notes that fixpoint-free global consequence remains reducible
to first-order consequence.  Adding a symbol with an input feed from `a`
toward `b` also destroys the specific separation (5), without by itself
proving completeness.

The feed carrier is an intrinsic directed graph, not generally a tournament:
it may omit pairs, contain both orientations, and have loops.  Completing it
to a tournament invents reachability and can invalidate the proof.  The
load-bearing invariant is the absent feed path `a -> b`, not a chosen total
orientation.

## 5. Verification record

The primary exact companion enumerates every set-valued interpretation
`f:M_b -> P(M_a)` with

```text
1 <= |M_b| <= 4,                  1 <= |M_a| <= 3.      (6)
```

It checks `5,050` models, finds exactly three `Gamma`-models (one for each
tested size of `M_a`, always with `|M_b|=1`), checks (4), and independently
closes the feed graph.  This finite audit is a hostile control for the
all-cardinality symbolic proof, not its replacement.

An independent referee rebuilt the denotations with literal finite sets,
audited every primitive rule and the one-sort minimality import, and returned
`PASS`.  Ordinary and optimized runs are required to agree byte-for-byte.

## 6. What this does not prove

THM-4090 proves no nonaxiomatizability result, no result about the paper's
least-fixpoint fragment, and no arithmetic theorem.  Treating rational,
irrational, and transcendental numbers as three sorts supplies only types;
transfer among them requires actual total operations and true axioms.  Bare
sort or tournament labels cannot certify arithmetic type.
