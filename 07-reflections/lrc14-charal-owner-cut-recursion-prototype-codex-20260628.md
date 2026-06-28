# LRC14 Charal Owner-Cut Recursion Prototype Reflection

Date: 2026-06-28

This pass started after pulling concurrent HYP-3410, which had already executed
the broad Bring/Schwarz-Christoffel/BDH/Menger/charal synthesis requested by
the user.  The productive move was to narrow that synthesis into an executable
finite lemma API rather than create another atlas of analogies.

After rebasing over the concurrent S258/HYP-3415 critical-path map, the right
positioning is sharper: this is not the main LRC14 completion route.  The main
route is q-witness for non-covering sets, LRC<=13 induction on `Q`, and the
single uniform decorrelation floor `|SPEC| < product`.  The owner-cut recursion
is useful as finite exception routing around that floor: it says which endpoint
owner labels, accessory debts, height factors, or named residuals must be
kept visible before any averaged or scalar inequality is trusted.

The later rebase over HYP-3416 and HYP-3417 sharpened the role again.  HYP-3416
is the recursive quotient-ladder abstraction; HYP-3417 is the signed
owner-current certificate synthesis; this HYP-3419 artifact is best read as
the charal decision-tree and purity-test harness that supplies finite cut
trees to HYP-3417 or reports the first quotient failure back to HYP-3416.

S259/HYP-3418 then changes how to read the final branch of the `10->20` tree.
The label `2:g2` is not just a convenient third separator; it is the visible
even-cover coordinate in the frontier cut.  If the covering floor is 2-adic,
then future owner-cut experiments should track whether every hard frontier
requires an even-cover label before the unit-petal residual disappears.

HYP-3410's strongest datum is the owner-cut readout:

```text
height leak:               one-label cut 5:g1
persistent owner leak:     one-label cut 1:g1
10->20 frontier leak:      size-3 cut, no common core label
```

That kills a tempting but false simplification.  The next theorem should not
say "there is always one owner label that separates the petal exit."  The
latest known frontier already asks for a bounded multi-label cut.  The exact
prototype decision tree for the `10->20` fiber is:

```text
test 13:g1
  yes -> positive-Haar-open
  no  -> test 11:g1
    yes -> positive-Haar-open
    no  -> test 2:g2
      yes -> positive-Haar-open
      no  -> unit-petal-named
```

This is a much more proof-shaped object than the famous-function prompts.  It
is finite, checkable, and compatible with a Menger/Farkas dual statement: every
mixed charal fiber should have a bounded endpoint-owner separator, a dual
owner-current certificate, or a named residual exit.

The charal-purity table also gives a warning.  On the small represented sample,
turn words and richer charal words often split the rows already.  That does not
make them proof predicates.  The full charal signature is close to a row label,
while a small owner-cut sidecar is a theorem-shaped compression.  The proof
should therefore prefer bounded cut sidecars over full row retention.

The outside ideas now have sharper roles:

- Bring: branch alphabet only.
- Schwarz-Christoffel: turn word plus accessory owner debt.
- BDH: finite label-variance priority order.
- Menger: actual bounded-cut theorem shape.
- Krasner: local owner/contact stability gate.
- Sophie Germain: height/flex factor channel after a live height debt is named.
- HLW: no-scalar-shadow guardrail.
- Ramanujan-Soldner: zero-level hygiene.
- Meissel-Mertens: tail entropy only after finite exceptions are named.

The new proof pull is concrete: extend the HYP-2963/HYP-3406 bank beyond
`(72,20)` and compute the minimal owner-cut size of the first
`residue+owner_support` failure, or of the first failure of any weaker charal
quotient.  If the minimal cut stays bounded by a small constant, try to prove a
finite owner-cut theorem.  If it grows, record the first growth pattern as
Schwarz-Christoffel accessory debt, tropical/off-grid debt, state-lift debt, or
a new finite residual.

Read in the HYP-3415 critical-path frame, the most useful follow-up is not a
parallel census.  It is a stress test for the floor theorem: identify the first
finite packet where a compressed charal/owner quotient would hide a theorem
exit, and either give a bounded cut-code certificate or mark the exact debt
that the final `R' > 0` proof must discharge.

Read in the HYP-3418 2-adic frame, the same follow-up should distinguish
binding-owner labels from even-cover labels.  A bounded owner tree that always
needs one even-cover coordinate would point toward the 2-adic descent route;
a tree that avoids even labels would be evidence that the obstruction has
already moved to a finite accessory sidecar.

Assumption challenged: famous named structures are not the vertices.  The
vertices that move the proof are proof modules: bounded owner-cut theorem,
decision-tree API, finite BDH label ranking, owner-stability gates, accessory
reconstruction, height-factor channels, and scalar-firewall hygiene.
