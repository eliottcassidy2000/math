# LRC14 True-Wide Cap-Floor Gate - Codex S60 Reflection

The useful move this session was not a new inequality, but a change of
currency.  HYP-2693 says true-wide rows should satisfy the raw Bonferroni4 cap
gate `U4<=cap_k`.  THM-535 says every cap has a universal subadditive floor
`(k-6)/7`.  The exact audit shows that true-wide rows appear to spend only
that floor for `k>=9`; the exact cap dividend is only needed at the exceptional
`k=8` edge.

That is a real structural simplification.  A proof of

```text
second_largest(E)>14, |E|=k>=9 => p0+p5+5p6 <= (k-6)/7
```

would close the HYP-2693 gate without knowing the cap minimizer.  This matters
because the cap minimizers are not uniform: `k=8` is anomalous, while `k=9..12`
follow the `{1}+top cluster` family.  The audit's behavior matches that split
exactly.

The exact scan result is sharp enough to be useful:

```text
k=8:  3 true-wide floor failures in B18, 0 cap failures.
k=9:  27020 true-wide rows in B18, 0 floor failures.
k=10: 3432 true-wide rows in B16, 0 floor failures.
k=11: 3003 true-wide rows in B16, 0 floor failures.
k=12: 2002 true-wide rows in B16, 0 floor failures.
```

The worst `k=8` floor debt is `107/8820`, while the dividend is `563/5880`.
The `k=9` leader has floor slack `29/980`; the `k=10` leader has floor slack
`29/588`.  Those two matching `29` numerators are probably signal: both leaders
are near-boundary true-wide rows with tail `1/14`, and both look like a
finite-scale overlap of a dilated AP core with the barely-far pair `{15,16}`.

The proof route should now avoid three mistakes:

1. Do not call this a proof of LRC14.  It is a sharper proof obligation.
2. Do not mix THM-534's optimal moment dual with THM-556's raw Bonferroni4
   `U4`.  Here `U4=p0+p5+5p6`.
3. Do not force boundary/AP rows into this gate.  AP9 and its doubled boundary
   copy fail raw level 4; HYP-2691's transfer-DP template branch is the right
   place for them.

The next concrete lemma should be a true-wide high-tail suppression statement.
HYP-2683's coarse state-mass profile suggests that high true-wide risk forces
low missed-state entropy/support.  HYP-2684 supplies the nonresonant BV/Fourier
decorrelation route.  HYP-2676/HYP-2677 supply finite routers for same-sign
packet pockets.  The cap-floor gate is where those tools can meet: prove the
floor on the decorrelated/high-state complement and reserve finite ledgers for
the low-state resonant rows.

Tournament Analysis in the script used proof obligations as vertices rather
than runners.  The resulting tournament was transitive, with AP/boundary floor
failures ranked first, the three `k=8` true-wide floor debts next, and the
`k=9` true-wide leader only after all cap-dividend rows.  That quotient matches
the mathematical split: boundary first, `k=8` dividend second, `k>=9` floor
gate third.
