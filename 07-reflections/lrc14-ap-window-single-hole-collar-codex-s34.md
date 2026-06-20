# LRC14 AP-Window Single-Hole Collar

This session closed the first exact proof obligation left by HYP-2651.  The
right object was not the hole position by itself and not a scalar invariant of
the core.  The useful proof carrier is the addressed wall gap

```text
R(v,a) -> L(w,b).
```

That is the bounded opposite-sign pattern in its most literal form: one danger
tooth ends with a right wall, the next begins with a left wall, and the safe
component is the signed determinant left between them.

For the drop-6 core the surviving components have determinant numerators

```text
[3, 5, 5, 3],
```

owned by the wall chain

```text
13 -> 12 -> 11 -> 12 -> 13.
```

This is more informative than saying the safe measure is `7/858`.  It says the
collar is a small signed boundary-determinant phenomenon in the high-speed
mouths opened by deleting `6`.

## What Changed

THM-541 now proves:

```text
min_e meas(G_{[1,13]\{e}}) = 7/858,
```

uniquely at `e=6`, with next value `426/35035` at `e=12`.

The exact finite certificate is:

```text
04-computation/lrc14_ap_window_single_hole_certificate_codex_s34.py
05-knowledge/results/lrc14_ap_window_single_hole_certificate_codex_s34.out
```

The theorem is local.  It does not close OPEN-Q-108, but it removes the first
AP-window proof obligation from the queue.

## Incoming Work Integrated

HYP-2650 and HYP-2652 both said that scalar summaries become faithful only
after retaining an address.  THM-541 is a small exact instance of that rule:
the measure table follows only after the boundary-owner ledger is kept.

The THM-538 dispute also matters indirectly.  It is another reminder that
dropping inactive or zero coordinates can fake a cancellation.  Here the
certificate avoids that mistake by carrying the actual fixed-observer walls and
their wrapped periodic addresses.

## Next Target

The next proof obligation is the near-collar mouth-retention theorem:

```text
If meas(G_C) < 426/35035 for a primitive positive 12-core C,
then C retains the drop-6 mouth template.
```

HYP-2653 corrected the naive exact-row version: the AP-tail row
`(1,2,3,4,5,7,8,9,11,12,13,20)` is already below `426/35035`, but it keeps the
old drop-6 mouths undamaged and adds `1/980` new mouth mass.  The certificate
therefore suggests a sharper route than raw enumeration: classify possible
low-measure R-to-L wall chains and prove that old-mouth damage, or the absence
of the central missing clock, costs at least the AP one-hole second value.
