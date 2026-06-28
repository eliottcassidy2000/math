# LRC14 Gamma Harmonic-Sieve Remainder

The useful Euler-Mascheroni move was not to import a constant into the proof.
HYP-3430 already makes that a scalar firewall: harmonic intercepts calibrate
tails but do not determine endpoint-spine class.  HYP-3431 then proves the
canonical corridor-fence case, and HYP-3432 shows reciprocal endpoint budgets
are only wall-ranking sidecars.  HYP-3433 then makes the canonical endpoint
tail finite part legal only after labels return.  The next question is what a harmonic
compression is forgetting in the general one-branch target from HYP-3426.  The
exact answer is now visible:

```text
branch0_mass = naive_slack + overlap_tax.
```

The naive slack is what remains after scalarizing the odd bad intervals by a
sum over owners.  It fails on `59/150` audited rows.  Every failure is repaired
by a positive overlap tax, and the tightest repair is the canonical
`{1..11,13,84}` row with tax/deficit ratio `1.090875`.

This is a good finite lemma target because it names the associativity failure
of the compression.  The set identity is perfectly associative, but once the
union is replaced by a scalar sum, intersections have been erased.  The erased
information is not philosophical; it is the exact rational overlap-tax term.

Euler-Mascheroni calibrates the denominator-prefix side of the story:
`H_N-log N-gamma` and the odd-prefix analogue stay in a small bounded band on
the audited rows.  That suggests harmonic endpoint budgets may be useful for
normalizing candidate inequalities, but the proof still has to be exact and
local.

The next route I would try is graph-theoretic.  Make a bipartite incidence
graph between `E_safe` components and odd bad intervals, weight incidences by
exact intersection mass, and view overlap tax as shared-incidence conductance.
HYP-3429's endpoint-spine rank `<=2` then looks like a bounded Menger cut or
low-rank conductance certificate.  A Green-current or algebraic-connectivity
bound could be relevant, but only if it retains endpoint labels and negative
leakage sidecars.

The clean theorem shape is:

```text
either naive_slack >= 0,
or a rank-2 endpoint spine certifies overlap_tax > -naive_slack.
```

That turns the vague "harmonic sieve" idea into a finite-ruler lemma with an
explicit correction term.
