# LRC14 AP84 Tail Component Phase

This session tightened the HYP-3451 graph target by taking its "AP-with-84m
tails are the danger family" sentence literally.  The useful surprise is that
there are two different clocks.

The proof clock starts at `m=5`.  From that point onward the best
component-cover escape is the same interval HYP-3433 found:

```text
[(14ceil(48m/7)+1)/(588m), (14ceil(48m/7)+13)/(588m)]
```

with labels `E:84m/E:84m` and length `1/(49m)`.  This is the good clock,
because it retains endpoint labels and a concrete safe component.

The danger-score clock is residue-noisy.  Extending beyond HYP-3451's first
twelve AP-tail rows moves the checked raw dead-fraction maximum to `m=35`,
where the escape-count correction is smallest.  That does not make `m=35` the
right theorem base.  It means raw dead fraction is an illegal quotient unless
it carries the mod-35 boundary count and endpoint labels with it.

The new candidate AP-tail proof shape is:

```text
1. prove m=1..4 as finite mixed E/B1 transient cases;
2. prove the rank-one E:84m/E:84m interval for all m>=5;
3. prove the mod-35 Beatty escape count from the low corridor floor counts;
4. use connected dead projection and pair-rank <=2 as graph sidecars, not as
   replacements for the interval certificate.
```

This is smaller than the full HYP-3451 saturation theorem but more surgical.
It should be the canonical base family for any Menger, Green-current, or
conductance proof of the component-cover obstruction.
