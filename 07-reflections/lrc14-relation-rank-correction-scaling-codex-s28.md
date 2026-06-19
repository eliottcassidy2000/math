# LRC14 Relation-Rank Correction Scaling

**Source:** codex-2026-06-19-S28, HYP-2640 / T888.

The useful disappointment from the rank atlas is that raw relation rank is too
eager.  It saturates on AP, near-AP, GAP, and the wide third-pocket rows.  That
would be bad news if the intended theorem were "many relations imply a large
correction."  But it is good news for the proof split: rank is not the meter;
it is the gate.

The picture now feels like this:

```text
no small relation rank
  -> orbit is close to independent
  -> large baseline margin;

saturated small relation rank
  -> inverse-combinatorics applies
  -> now inspect the signed visible quotient.
```

The reason this matters is that the third pocket is exactly saturated but not
dangerous.  At `k=8`, AP and third-pocket A both have exact height-2 rank `6`
on the nonzero coordinates.  Their `L_y` corrections are not close:
`0.308965` for AP and `0.013547` for third-pocket A.  The difference is carried
by coherence: AP has `12` fold motifs and `1786` exact height-2 relations;
third-pocket A has `3` folds and `326` exact relations.  The third pocket has
the rank of a relation-rich object but the signed visible mass of a loose one.

That distinction is the coimage lesson in a cleaner language.  The integer
relation lattice maps to a finite mod-27 coimage, and then to the seven-sector
signed functional.  Rank only sees the size of a domain or image.  LRC sees the
trace after the signed functional.  A large kernel can have huge absolute mass
and still small signed correction; this is exactly the same pattern as the
two-large reciprocal tail, where the raw absolute envelope is enormous until
we keep additive channels and characters visible.

So the next proof should not try to squeeze everything through one scalar.  It
should be a two-stage statement:

1. low relation rank or uncovered weighted-fibre vertex gives an independent
   peel;
2. saturated rank gives inverse structure, and non-AP inverse structure must
   lose visible fold/coimage coherence.

The AP is then not merely "the densest relation lattice."  It is the row where
the relation lattice, observer folds, and mod-27 signed coimage all point in
the same direction.  Non-AP rows can keep one or two of those resources, but
the atlas suggests they cannot keep all three.  That is the shape of the
remaining lemma.
