# LRC14 Two-Adic Descent Loss Ledger

This session made the HYP-3418 correction concrete.  The promising object is
not "odd runners are safe at `t=1/2`."  That is precisely the trap: at
`t=1/2` every even speed is dead.  The promising object is a descent packet
that remembers the halved even child, the odd blockers at the displaced
witness, and the owner-current labels that explain why a quotient is allowed
to forget anything else.

The exact AP-collar scout is not a covering-floor proof.  It uses HYP-3410's
and HYP-3417's non-covering mixed fibers because those rows have rich labelled
owner data and exact row names.  That guardrail matters: all `23` reconstructed
rows still have the elementary `q=14` witness.  The contribution is not direct
coverage; it is a test of what a future covering packet must carry.

The output is sharper than expected:

```text
half_witness_failure_rows=23/23
rows_with_even_binders=22/23
rows_where_even_child_carries_M_at_u=2t*=22/23
```

So the raw half witness is always wrong on this substrate, but the halved even
child almost always retains the actual bottleneck.  The lone exception,
`two drop(5,10)->add(15,20)`, lands in the new class
`odd_binder_after_even_shift`: after shifting off the even death point, the
actual binders are odd.  That should become a named discharge case rather than
a failure of the 2-adic route.

The owner-current bridge also became more literal.  HYP-3417's hardest current
is `{2:g2,11:g1,13:g1}`.  The `2:g2` label is the even-cover hinge; the two
`g1` labels are binding-side owners.  This says the proof interface should not
be just "v2 layer" or "owner cut."  It should be:

```text
halved even child + odd-blocker ledger + bounded labelled owner current.
```

This also absorbs the q-Pochhammer/modular-function prompt in a conservative
way.  The useful analogy is not a new special function calculation; it is the
finite negative-tail rule.  A q-expansion may forget everything except a finite
principal part only after the cusp data is declared.  Likewise, a halving
quotient may forget raw row order only after its finite loss ledger is
declared.

After fetching mainline, that exact lift has already started in a stronger
form.  HYP-3422 audits genuine covering packets by splitting

```text
S = O union 2E
```

and moving from the dead `t=1/2` witness to two odd branches under `u=2t`.
HYP-3425 then rewrites the target as a Helly-style claim that `E_safe` is not
contained in the two-color odd bad core.  A later fetch added HYP-3426, which
uses the mirror involution to reduce the target to one branch, and HYP-3427,
which replaces bare positive survivor mass by exact wall-signature
certificates.  That mirror/wall chain is now the direct proof route.

So the role of this HYP-3428 ledger is now narrower and more useful: it is the
controlled-forgetting audit for the relocation proof.  If a covering packet
fails the HYP-3426 mirror target or the HYP-3427 wall-word target, classify the
failure as one of the named loss classes here, an owner-current/even-hinge
debt, an off-grid debt, a sheet debt, or a state-lift debt before adding any
new analytic analogy.
