# LRC14 AP84 Mod-35 Floor Count

HYP-3456 is the cleanest AP-tail improvement in this run.  HYP-3454 had a
checked mod-`35` escape clock, but it was still expressed as an empirical
correction vector.  The new script derives that vector from the HYP-3431 low
corridors and the moving `42m` grid.

The key observation is that, once the low core is reduced to

```text
[8/49,6/35] and [29/35,41/49],
```

the high even half-speed contributes only safe gaps

```text
G_k(m)=[(14k+1)/(588m),(14k+13)/(588m)].
```

Counting gaps that intersect `[8/49,6/35]` gives

```text
N(m)=floor((504m-6)/70)-floor((96m-13)/14),
```

and the mirror corridor doubles it.  This recovers the HYP-3452/HYP-3454
escape count through `m=70` with no component-audit failures and gives the
closed shift `escapes(m+35)-escapes(m)=24`.

HYP-3457 now closes the finite `m=1..4` mixed cases as explicit four-window
packets.  The remaining AP-tail proof is mostly carrier bookkeeping: prove or
import the HYP-3431 fixed-corridor identity as the low branch-union carrier,
use HYP-3454 for the endpoint interval, use HYP-3456 for the escape count, and
splice those named sidecars into HYP-3439.  This is still local to the AP84
tail, but it removes the sampled ingredients from the HYP-3439 descent.
