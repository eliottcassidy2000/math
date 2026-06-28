# LRC14 AP-collar finite lemma certificate

This session pushed HYP-3401 from a quotient-mixing scout into a finite lemma
certificate.  The useful shift is small but important: the AP collar is no
longer just evidence that a quotient leaks.  It is now a finite theorem target
with exact witnesses.

The executable certificate checks `924` AP one-swap collar rows through
replacement speed `84`.  AP and Goddyn-Wong `12->24` are the only
boundary-tight rows.  Every other row has a rational strict-open interval, and
the uniform strict-open mass floor is exactly `1/1260`, uniquely at `12->36`.
The certificate digest is
`c40c24d7746f05a708a9b625afeedcfae5d6fff8e8e39ba892b4325cd5b1e148`.

The information-theory read is direct.  A quotient is a compression function;
the target is `exit_status`.  The HYP-3311 nonunit-height packet still has a
mixed boundary/strict fiber of size `31`, so `exit_status` is not a function of
that packet.  The sharp pair is AP versus strict-open `13->27`.  They share the
six unit-contact status, C3 skeleton, quadratic `Q(sqrt(-7))` character,
covering layer, and nonunit height data.  Their first visible difference is a
unit-height lift:

```text
(13,0) -> (13,1)
```

This is the concrete failure-of-compression analogue of the older law-defect
stories.  The destroyed property is not commutativity, associativity,
invertibility, or a field shadow.  It is a unit-height coordinate.  The repair
matrix says `unit_height_flex`, `full_height_flex`, and `height_completed_packet`
kill the mixed fiber, while unit contacts, covering layer, and nonunit height
do not.

The next proof route should use HYP-3405 as a base lemma, not as a terminal
invariant.  Full height retention is too close to row identity for the global
proof.  The right global replacement should be a chamber theorem: every
height/flex move either lands on AP/GW boundary, opens strict safe mass, enters
`Phi14d` equality, discharges through Toeplitz/Green/root motion, lifts to
state debt, or names a residual.

Assumption challenge: the tournament vertices were not runners, residues,
gaps, replacement speeds, or fixed circle sections.  Those are coordinates
inside the certificate.  The preserved predicate is finite-lemma proof
legality, so the vertices were proof carriers: strict-open witness
certificate, height-completed oracle, unit-height obstruction vector, sidecar
repair matrix, boundary atom classifier, and raw mixed-fiber scout.

The immediate next finite step is to turn the generated certificate into a
short formal lemma statement that a Lean-side or paper-side proof can consume:
two named boundary rows, a rational interval witness for every non-boundary
row, and a one-line obstruction vector showing why nonunit-height compression
cannot be used without a unit-height sidecar.
