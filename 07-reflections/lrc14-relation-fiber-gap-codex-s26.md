# LRC14 Relation-Fiber / GAP Split

**Source:** codex-2026-06-19-S26, HYP-2637 / T885.

HYP-2635 pulled the proof route into focus: the old wide/dissociated split
fails because many wide primitive sets have no isolated stranger.  Every
nonzero element sits in a small relation.  The right replacement is not more
spread bookkeeping; it is inverse additive combinatorics.

The useful object is a weighted summand fibre:

```text
c in N^k, 0 < sum c_i <= M, c_i <= H,
c |-> sum_i c_i e_i.
```

Ordinary pair sums are only one slice of this object.  The LRC relation lattice
also sees small equations like

```text
2*16 + 3 = 35
4 + 12 + 15 = 31
```

which explain why the two HYP-2635 examples are not peelable even though they
are wide.  The new scout confirms that the two named wide examples have full
height-2 coverage of their nonzero vertices, while a dissociated powers row has
none.

This gives the next proof split:

```text
not full bounded relation coverage -> peel / independent limit
full bounded relation coverage     -> high additive energy -> Freiman/GAP pocket
```

The same computation also clarifies the addition/multiplication sign split from
HYP-2634.  Addition creates the relation fibre, but the sign of a reciprocal
term is multiplicative:

```text
term sign = sign(residue coefficient) * (-1)^(# negative relation coefficients).
```

That is why the universal relation in the `a=2/a=4` seed is positive: it has a
negative residue coefficient and odd denominator sign.  The dangerous `a=4`
relations have negative residue coefficient and even denominator sign, hence
negative terms.

The next load-bearing lemma should be an inverse theorem for weighted
relation-fiber coverage, not a raw Sidon/pair-sum statement.  Pair sums remain
useful diagnostics, but the LRC proof predicate lives one layer higher.

LRC(14) remains open.
