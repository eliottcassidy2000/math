# Independent referee report: two-sort matching-logic obstruction

## Verdict

**PASS**, with two non-blocking wording corrections before promotion.  The
semantic entailment, satisfiability witness, nonderivability induction, and
minimal positive sort count are correct for the exact standard many-sorted
system defined in `ml-onesort-localization.tex`.

The candidate is

```text
sorts: b, a
symbol: f : b -> a
Gamma = {forall x:b forall y:b. f(x and y)}
phi   = forall x:b. x
```

Both displayed patterns are closed.  `Gamma` has sort `a`; `phi` has sort `b`.

## Exact semantic audit

The paper's clauses (source lines 1580--1612) give nonempty carriers, element
variables as singletons, conjunction as intersection, pointwise symbol
extension with an empty result on an empty argument, and universal
quantification as intersection.

For a valuation with `x=r`, `y=s`,

```text
[[x and y]] = {r} if r=s, and empty otherwise.
```

Consequently, if `M_b` contains distinct `r,s`, the `(r,s)` factor in the
universal intersection defining `Gamma` is empty at sort `a`.  Since `M_a` is
nonempty, `Gamma` is not total.  Thus every `Gamma`-model has `|M_b|=1`.

Conversely, if `M_b={r}`, then

```text
[[Gamma]] = f_M(r),
```

so `Gamma` is total exactly when `f_M(r)=M_a`.  This supplies a model for every
nonempty `M_a`.

Finally,

```text
[[forall x:b. x]] = intersection_{r in M_b} {r},
```

which equals `M_b` exactly when `|M_b|=1` (carriers are nonempty).  Hence
`Gamma |= phi`, and satisfiability is non-vacuous.

## Complete primitive-rule audit

The exact many-sorted system is Figure 1 read with sorts (source lines
1623--1633).  Induct on a finite derivation, strengthening the claim to:

> Every line of sort `b` has an empty-theory derivation.

The cases are exhaustive:

| Item | Audit for a sort-`b` conclusion |
|---|---|
| Hypothesis | Impossible: the sole hypothesis has sort `a`. |
| (1) propositional tautology | Already an empty-theory axiom instance. |
| (2) modus ponens | Both premises and the conclusion have the same sort; replay MP on the empty-theory derivations supplied by induction. |
| (3) existential-quantifier axiom | Already an empty-theory axiom instance. Variable substitution is same-sorted and is not a premise-bearing cross-sort rule. |
| (4) existential generalization | Premise and conclusion have the common pattern sort; replay the rule after induction. The quantified variable may have either sort, but quantification never changes the pattern sort (lines 1592--1593). |
| (5)--(7) propagation | These are empty-theory axiom instances. |
| (8) framing | Its premise has an argument sort and its conclusion the symbol's result sort. The only symbol is `f:b->a`, so framing concludes only at `a`, never at `b`. |
| (9) Existence | `exists x:s. x`, not `forall x:s. x`; it is an empty-theory axiom. |
| (10) singleton-variable | An empty-theory axiom for every well-sorted context instance. |
| derived universal generalization | Sort-preserving because quantification does not change the body sort; replay it after induction. |

Thus `Gamma |- phi` would imply `|- phi`.  Many-sorted soundness is stated at
source lines 1631--1633.  But on `M_b={r,s}`, with any nonempty `M_a` and any
interpretation of `f`, `[[phi]]=empty`, so `phi` is not valid.  Therefore
`Gamma !|- phi`.

This also resolves the possible `forall x:b. x` trap.  Universal
generalization can indeed be applied to hypotheses (the paper explicitly notes
`{x} |- forall x.x` at lines 357--360), but there is no hypothesis or derived
line `x:b` here.  The sole hypothesis has sort `a`, and no primitive inference
has an `a` premise and a `b` conclusion.  `forall x:b.x` is not the Existence
axiom, and soundness plus its two-point countermodel independently rules out an
empty-theory derivation.

## Minimality and scope

The paper's one-sort global-completeness corollary (source lines 1059--1066)
rules out a one-sort counterexample in the same basic, definedness-free,
fixpoint-free system.  This two-sort witness therefore proves that two is the
least positive cardinality of a sort set for which **some signature and
theory** make this standard system globally incomplete.

This is existential, not uniform over all two-sort signatures or theories.  It
does not cover extensions by nominals, set variables, fixpoints, definedness,
or genuinely global/cross-sort rules.

## Requested corrections

1. `REPORT.md` line 251 should not call the two-sort theorem literally an
   "instance" of the paper's Corollary 35.  That corollary is stated over the
   particular three-sort signature of Proposition 30 (source lines
   1813--1824).  Say instead: **"The same sort-flow induction proves the
   following two-sort variant of Corollary 35"**, or state the general
   sort-separation lemma explicitly.
2. In the minimality sentence, retain the full language/system scope:
   **basic, definedness-free, fixpoint-free matching logic, with ordinary
   set-valued symbols and no nominals or set variables, using the displayed
   standard Hilbert system**.  Without this qualifier, the paper's one-sort
   nominal counterexample with the same conclusion `forall x.x` would make the
   sentence false.

The sentence about "necessitation" at report lines 248--249 is unnecessary:
necessitation is not a separate primitive rule in the displayed calculus.
The paper's derived necessitation proof (source lines 756--770) also uses a
tautology, modus ponens, `Propagation_bottom`, and propositional composition;
framing is its only sort-changing step.  If retained, say that derived
necessitation can follow a well-typed symbol/feed path but cannot reverse one.

## Computation and repository checks

The supplied bit-mask probe is exact for its stated finite universe.  It
enumerates all `5050` maps `M_b -> P(M_a)` for `1<=|M_b|<=4` and
`1<=|M_a|<=3`, finds exactly the three singleton-`b` models, and has no
entailment failure.  Its recorded source/output SHA-256 values match, and
normal/optimized outputs are byte-identical.

`independent_probe.py` repeats the semantic audit with literal frozenset
denotations, checks the closed forms on every interpretation, and separately
computes the feed closure `{a->a, b->b, b->a}`.

Repository search found no matching `MISTAKE-*` entry.  The only exact
THM-4090 canon file is still `RESERVED / UNPROVED EMPTY STUB`; it contains no
live theorem statement or proved dependency.  After fetching, `origin/main`
added only THM-4092 and did not alter this conclusion.
