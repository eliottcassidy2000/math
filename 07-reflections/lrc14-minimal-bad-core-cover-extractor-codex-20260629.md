# LRC14 Minimal Bad-Core Cover Extractor Reflection

HYP-3436 turns HYP-3435's "emit a minimal cover if a component is covered" hook
into an exact object.  The useful carrier is not a row scalar and not a runner
tournament.  It is the local ledger

```text
bad-core component
+ minimal B0 odd-owner cover
+ minimal B1 odd-owner cover
+ even endpoint walls
+ survivor gaps.
```

The audit on `135` primitive covering rows found survivors in every row.  That
is still evidence, not a proof.  The improvement is that the obstruction side is
now finite and labelled: `11670` bad-core components, mostly `(1,1)` cover
pairs, with rare hard components requiring total cover size `6`, and endpoint
support size never above `3`.

The tight canonical row `{1..11,13,84}` is especially informative.  It has
`22` fully bad even-safe components and only `4` mixed components, yet those
mixed components leave survivor mass `1/105`.  This matches HYP-3431's
corridor-fence picture: proving the floor is not about making the bad core
small; it is about proving the local bad-core covers cannot glue across every
even-safe component without leaving a survivor corridor or paying named
sidecar debt.

## Proof Pull

The next theorem should be a local-to-global gluing obstruction:

```text
For primitive covering rows, the emitted minimal B0/B1 cover ledgers cannot
cover all E_safe components compatibly unless the row routes to an existing
sidecar/debt.
```

Candidate exits are HYP-3431 corridor-fence structure, HYP-3429 endpoint
spines, HYP-3428 two-adic loss debt, HYP-3427 wall words, HYP-3426 mirror
reduction, HYP-3417/HYP-3420 owner-current, HYP-3423 topology legality, and
HYP-3421/HYP-3129 signed-SPEC/Rprime.

After rebasing over HYP-3437, the handoff is clearer: HYP-3436 supplies the
bad-core atoms and owner-cover signatures, while HYP-3437 should build the
incidence/Menger graph that explains negative-slack overlap tax and multi-owner
no-gluing.

## Assumption Challenge

Alternate vertices considered: runners, gaps, fixed sections, section
boundaries, wall crossings, residues, endpoint owners, branch bad intervals,
even-safe components, bad-core components, survivor gaps, and proof
obligations.

Chosen quotient: exact bad-core components plus minimal odd-owner cover words.
It preserves the two-adic branch-overlap predicate.  It destroys global row
geometry unless paired with endpoint gates and survivor-gap placement.  Raw
bad-core measure, harmonic tails, and topology summaries remain guardrails or
ordering sidecars, not certificates.
