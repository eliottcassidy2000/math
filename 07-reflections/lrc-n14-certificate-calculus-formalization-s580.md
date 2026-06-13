# LRC n=14 certificate calculus formalization (S580)

The most useful thing formalization did today was stop the proof program from
being a cloud of nearly-right metaphors.

After the upstream S572 work, the global spine is no longer vague:

```text
no multiple of n -> t=1/n witness
multiple of n    -> Cprime: positive measure
Cprime           -> LRC(n)
```

THM-398 then proves the dominance/long-interval part of Cprime, leaving only
the all-short small-multiple residual.  In parallel, HYP-2095/HYP-2100/HYP-2101
say that the n=14 fixed-boundary route is usually cheap before it is delicate:
unblocked small pairs fire in the transversal census, the unit-lift audit, and
the apex-lift sheaf site.

The rebase brought in HYP-2104, which names the same Cprime split as the Vitali
handoff: `n|v` is the transition from construction (`t=1/n`) to measure
(Vitali covering).  That makes the all-short residual more precise: it is
short-interval arc alignment against the `1/(nw)` lattice.

Those two stories look different only before formalization.

## The common residual

The all-short Cprime/Vitali residual says:

```text
one AP of thin arcs from v=nw is aligned enough to cover every component of G(S')
```

The failed-gluing sheaf residual says:

```text
every local cheap-pair section is blocked by shields or endpoint anchors
```

Both are cover statements.  Both should expose endpoint owners.  Both should
feed THM-397/HYP-2095 style private pivots.  So the sharper extension is:

```text
endpoint-cover circuit positivity:
  failed cover/gluing -> owner circuit -> private pivot -> peel -> positive measure.
```

That is a better proof target than "finish the 64 rows."  The 64 fixed fibres
still matter, but only after they carry owner labels and section restrictions.

## Why the script is small

`04-computation/lrc_n14_certificate_calculus_s580.py` does not search speed
sets.  It records the current gates as formal objects:

```text
SpeedRow -> CertificateSection or NamedResidual
```

with section kinds:

```text
witness_1_over_n
cheap_pair
positive_measure
ledger_failure
residual
```

This is the user's cascade language in proof form: each gate contributes a
conditional clearance factor.  The level-3 product is the repeated multiplication
closure of these local clearances.  A zero factor is no longer "mysterious bad
row"; it must be named as a residual and pushed through restrictions.

## Tournament Analysis choice

I deliberately did not use runners as vertices.  I also did not use the 64
classes as vertices.  For this session the vertices are proof obligations:

```text
G0 no-multiple n-clock
G1 Cprime reduction
G2 dominance/long-interval dodge
G3 all-short Cprime residual
G4 paired-or-anchored cheap-pair gate
G5 fixed-boundary owner functor
G6 unit-lift cheap sieve
G7 apex-lift certificate sheaf
G8 endpoint-cover circuit positivity
G9 tie-wall limit functor
G10 certificate-calculus closure
```

The pair observable is:

```text
(open_rank, residual_weight, dependency_count, -evidence_strength, tie_order)
```

The switch orients toward the harder unresolved obligation.  The fingerprint is
transitive: no directed 3-cycles, singleton SCCs, and one Hamiltonian path.  That
is expected.  This is a proof ledger, not yet the internal monodromy graph.  The
interesting cycles should appear only after S579's certificate transports are
attached inside fixed-owner fibres.

## Concrete next tests

1. Build an all-short `v=14` Vitali-alignment table.  For every component of
   `G(S')`, record which `14`-arc covers it and which other-speed endpoints are
   adjacent.  Check whether every full cover produces a THM-397 owner circuit.
2. Build the promised fixed-fibre sheaf table for AP, `V*`, transversal flips,
   and minimal `gcd 3/gcd 9` lifts.  Columns should be cheap pair, reduced sum,
   mod-7 lane owner, mod-27 shell, shield blockers, endpoint anchors, D/U/N
   private pivots, exact `M`, and safe-measure flag.
3. Treat S579 union-only `{ledger,ledger}->cheap` rows as two-step composites.
   Those are not defects in the sheaf picture; they are evidence that local
   ledger failures can glue into a cheap global section.
4. Keep the analytic branch alive.  The component-discrepancy route might prove
   Cprime without the fixed-fibre machinery, but it needs a real bound on the
   number and denominator structure of components of `G(S')`.

## Assumption challenged

The tempting assumption is that the proof vertex set is predetermined: runners,
arcs, or self-converse classes.  It is not.  This quotient preserves exactly the
current predicate:

```text
every row exits by witness, cheap section, positive measure, or named residual.
```

It destroys exact speed values and round-class identity.  That loss is useful
because the frontier is now about how certificates transport, not about how
large the ambient class space is.

## Handoff

HYP-2111 should be read as a formalization and extension program, not a theorem.
The best next proof attempt is endpoint-cover-circuit positivity, because it is
where HYP-2103/HYP-2104's all-short Vitali residual and HYP-2101's sheaf residual
meet.
