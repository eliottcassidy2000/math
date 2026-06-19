# LRC14 Block-Frequency Transfer - Codex S24

The useful reframe from this pass is simple enough to say out loud:

```text
do not count the six support variables first.
split the exact relation into a core channel and a pair channel first.
```

The old absolute envelope saw a six-dimensional reciprocal cusp and panicked.
The block-frequency transfer sees a finite list of exact additive channels
`s`, each carrying a `6 x 6` residue matrix.  That keeps HYP-2632's signed
`chi_7` / affine / `Q` packet alive long enough for cancellation to happen.

The strongest numerical clue is the zero-lane row.  At `H=24`, the raw atom
absolute envelope is about `1420` times the signed mass.  If the atoms are first
summed into exact `s` channels and residue matrices, the entrywise block
envelope drops to `18.9`.  That is still too big to be a theorem by itself, but
it is the right kind of collapse: large absolute mass, small signed mass, and a
concrete side channel explaining where the mass went.

I added a same-residue spread core `(1,8,15,22)` to make sure the scout was not
only exploiting four identical exact speeds.  The collapse weakens, as it
should, because the core now spreads across many more `s` channels.  But it is
still there: raw/signed sits around `185-302`, while block `L1/signed` stays
around `14-17`.  That makes the reframe feel structural rather than cosmetic.

This also explains why the KPS S11 distribution scout feels adjacent.  There,
consec/AP wins only after keeping the full empty-sector distribution `p_t`;
here, the two-large tail only simplifies after keeping the exact frequency
channel and residue matrix.  Same grammar:

```text
retain the coimage/distribution/cochannel
then apply the signed estimate
then scalarize.
```

The later KPS `p_6` maximality scout strengthens that bridge rather than
changing it: AP/consec extremality is a correlation endpoint of the whole
distribution.  The reciprocal-tail analogue should therefore look for the
endpoint inside the channelized transfer operator, not in the raw atom envelope.

The next proof object should be a pair-line lemma.  For fixed `A,B`, each
`A*x+B*y=-s` channel is a one-dimensional arithmetic progression.  Away from
`s=0`, `1/(xy)` splits into two arithmetic harmonic sums; at `s=0`, it becomes
a quadratic line.  That is the place for the cotangent/Dedekind bound the user
has been pointing at.  The finite HYP-2632 packet supplies the phase.  The pair
line supplies the analytic handle.  The core table supplies the remaining
successive-minima accounting.

I do not think this finishes LRC(14) yet.  I do think it removes a false
difficulty: the residual is not a formless six-variable cusp.  It is a transfer
operator with small finite residue rank, and the proof should be written in
that language.
