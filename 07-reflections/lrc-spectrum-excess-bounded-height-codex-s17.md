# LRC Spectrum Excess and Bounded Height - Codex S17

The useful scalar is not "below the doubled-top mediant."  It is the excess

```text
e = p(k+1)-q, where M(S)=p/q.
```

That one integer separates the two phenomena that had been blending together.
Rows with small absolute gap are rows with large `q/e`; the really dangerous
ones are unit-excess rows with denominator proportional to a large height.

The AP-defect ladder is behaving like a bounded-height constant improver, not
like an asymptotic threat.  The known branches

```text
2/(2k+1), 3/(3k+2), 4/(4k+3)
```

all have `e=1`, so they improve the reciprocal depth constant cleanly.  But in
the S17 probe, trying to push past `r=4` immediately raises the excess:
`r=5,6` give `e=3,5`, `r=7` falls into a `1/k` row, and larger `r` loses
normalized depth.  That looks like a local residue ladder, not a runaway
sequence.

The pulled KPS denominator lemma is the reframing I needed: if `M=p/q`, then
`q <= 2 max(S)`, so

```text
M - 1/(k+1) >= 1/(2 max(S)(k+1)).
```

So any true `o(1/k^2)` dip has to make `max(S)/k` grow.  This converts the
lower-bound question from "search all families" into "find small-excess,
height-escaping extremizers."  The bounded AP-defect ladder can keep sharpening
constants forever and still not change the order unless `r` grows with `k` and
keeps excess under control.  S17 saw no sign of that.

There is also a nice correction to the user's first heuristic.  The doubled-top
family is universal as an upper packet, but the true second point is already
sporadic in small exact boxes: KPS S9 finds `k=6` at `5/33` and `k=7` at
`3/23`, below the doubled-top mediant.  That does not weaken the main lower
bound story; it makes the excess ledger more necessary.

The next proof-shaped move is modest but real: prove the upper cover for the
`r=3` residue classes.  The witness half now has explicit formulas:

```text
k == 7  mod 30: (3k - 1)/(5(3k+2))
k == 13 mod 30: (6k + 7)/(5(3k+2))
k == 19 mod 30: (6k + 1)/(5(3k+2))
k == 25 mod 30: (3k + 5)/(5(3k+2))
```

Those are the numerator-addresses for the lower bound `M>=3/(3k+2)`.  A modular
cover proving no crossing clears more than `3/(3k+2)` would turn a numerical
branch into a theorem.

Tournament-wise, this is another case where runner vertices are a distraction.
The right vertices are proof obligations: excess-one classifier, height-escape
filter, residue witness, upper-cover certificate.  That quotient preserves the
predicate we care about, the gap scale above the AP floor, while discarding the
microgeometry that does not affect the asymptotic question.
