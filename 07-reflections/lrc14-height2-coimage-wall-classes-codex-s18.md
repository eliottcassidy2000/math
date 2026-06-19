# LRC14 Height-2 Coimage Wall Classes - Codex S18

The coimage lens paid off again, but with a warning label attached.

The first temptation was to ask for more exact row clearing: enumerate
height-2 walls, compute `meas(S7)`, hope the ledger closes.  That is still a
good finite task, but it is not the most informative first question.  The
better question is: which analytic coimage classes do low-height walls already
explain?

That question has a surprisingly sharp answer.  For the `k=8` and `k=9`
ambient layers, height `<=2` one-large support-six walls touch every nonzero
projective coimage class.  In other words, the entire nonzero signed-mass
coimage table is already visible from very small additive relations.  The
tail, if it survives after finite wall accounting, cannot hide behind a new
residue address in those dimensions.

The `k=10` layer is different, and more interesting.  Height `<=2` walls hit
about `84%` of signed coimage mass but miss `31` nonzero classes.  The largest
misses are almost all repeated-residue packets:

```text
(1,1,1,1,a,a)
(1,1,1,1,a,b)
```

plus a small zero-cusp halo.  That feels like the true shape of the remaining
analytic theorem.  Four equal residues form the persistent cusp; the last two
coordinates carry the actual oscillation.  A generic lattice-volume theorem is
still too blunt for that shape.  The proof should specialize to repeated-root
cotangent or Dedekind sums.

This also reframes the large-absolute/tiny-signed clue.  The absolute mass is
not merely too large; it is large in exactly the classes where the coimage has
forgotten most of the row geometry.  The signed object only remembers the
projective residue packet.  Once we see that, the "tail" is not a dark ocean.
It is a repeated-residue packet with high cancellation ratios.

The caveat matters.  A coimage class can be wall-addressed and still occur in a
high-height tail.  So S18 does not delete those classes from analysis.  It says
which finite wall ledger must be accounted for first and which classes deserve
the next analytic theorem.

KPS S9's THM-539 update makes the warning more structural.  In the max-min
spectrum, the useful escape is not "large numbers are bad" but a precise
binding-pair resonance `{2a-1, a(k-1)}` plus an lcm-killer condition on the
coarse clocks.  The LRC14 coimage tail feels like the same phenomenon after
quotienting: repeated packets `(a,a,a,a,b,c)` are where many support clocks
have lost individuality, leaving a small signed residue problem behind.  That
suggests the next proof attempt should be less like a bigger bounded-height
search and more like a resonant tail theorem for repeated coimage roots.

Tournament-wise, the useful vertices are not runners.  They are not even raw
supports.  They are proof quotients:

```text
height-2 wall-addressed classes
> height-1 wall-addressed classes
> repeated-residue tail packet
> coimage atlas
> signed reciprocal theorem
> raw supports
> raw runners.
```

That quotient preserves the support-six tail predicate and discards the time
geometry.  For the current proof obligation, that is the right loss.
