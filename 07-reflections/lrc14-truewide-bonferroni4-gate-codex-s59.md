# LRC14 True-Wide Bonferroni4 Gate - Codex S59 Reflection

The incoming sector-sieve work showed the right warning: level-4 Bonferroni
fails on AP-like rows but succeeds on the named true-wide leader.  The useful
reframe is that level 4 is not just a numerical sieve.  For six inner sectors,

```text
1 - S1 + S2 - S3 + S4 = p0 + p5 + 5p6.
```

So level 4 cancels the entire middle profile `p1,p2,p3,p4`.  The only price
above the desired `p0` is the high missed-sector tail.  This gives a clean
proof-shaped split:

```text
true-wide branch:  prove p0+p5+5p6 <= cap,
boundary/AP branch: route through finite low-state templates.
```

The exact scan supports this split.  AP9 and its doubled boundary copy have
the same missed-sector profile and fail the level-4 gate, with tail
`53/392`.  The k=9 true-wide leader has tail `1/14` and still leaves exact
slack `3338/35035`; the k=10 true-wide leader also has tail `1/14` with slack
`629/7644`.

This connects back to HYP-2691.  The AP-prefix transfer DP already showed that
the largest local residual pressure is a finite AP6 append template, not a
generic true-wide obstruction.  HYP-2693 says the final-row sieve has the same
shape: AP-like dilation can preserve large high-tail mass, while true-wide
rows appear to suppress it enough that four Bonferroni levels suffice.

Incoming HYP-2692 adds a compatible leading-residual view: the apex-divisor
arithmetic organizes the resonance classes, but the lever is the leading
summed residual.  This high-tail gate is the final-row sector-sieve form of
the same split after the low-state and leading-residual pieces are routed.

The main proof obligation is now more exact:

```text
second_largest(E)>14, span(E)>14
    => p0(E)+p5(E)+5p6(E)<=cap_|E|.
```

Any proof of this should keep sector-state information long enough to detect
the finite resonant exceptions.  Cardinality distribution is enough for the
upper bound, but not enough for classification.
