---
id: THM-372
name: tournament-analysis-switchboard-well-defined
status: PROVED
date: 2026-06-01
session: codex-2026-06-01-S500
depends_on: []
---

# THM-372: Tournament Analysis switchboards are well-defined tournaments

## Statement

Let `V` be a finite labelled set.  Fix a strict total order `<_P` on `V`,
called the tie path.  For every unordered pair `{i,j}`, suppose we are given a
switch value

```text
s_ij in {-1,0,+1}
```

where `+1` means "orient by the displayed order `i -> j`", `-1` means
"orient by the displayed order `j -> i`", and `0` means "tie".

Resolving ties by the fixed path order

```text
i -> j  iff  i <_P j
```

produces a unique tournament on `V`.

If the switch values depend on a parameter `t` and each pair has only finitely
many switch-change times on an interval `I`, then the resulting tournament
`T(t)` is constant on the connected components of `I` minus those change times.

## Proof

For each unordered pair `{i,j}`, exactly one of three cases holds.

If `s_ij=+1`, orient the pair one way.  If `s_ij=-1`, orient it the other way.
If `s_ij=0`, the strict total order `<_P` chooses exactly one of `i<_P j` or
`j<_P i`, hence exactly one orientation.

Thus every unordered pair receives exactly one directed arc, which is precisely
the definition of a tournament.

For the parameter statement, outside the finite set of change times every
pairwise switch value is locally constant.  Since the tie path is fixed, every
arc is locally constant.  Therefore the whole tournament is constant on each
connected component of the complement of the change set.

## Significance

This is the formal base of the repo's Tournament Analysis pipeline:

```text
pairwise data -> switch functional -> tie path -> tournament movie.
```

It separates the mathematical object from the choice of gauge.  Rank-style
gauges, edge-local gauges, basketball pass gauges, LRC pressure gauges, and
metric threshold gauges all become instances of the same finite construction
once their pairwise switch values and tie path are declared.

## Related

- HYP-1931
- HYP-1932
- HYP-1940
- HYP-1941
- `01-canon/definitions.md`
- `04-computation/tournament_analysis_framework_s471.py`
- `04-computation/tournament_analysis_switchboard_s454.py`
