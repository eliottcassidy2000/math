# LRC14 Few-Apex Lift Packets

This pass attacks the covering bucket left by the labelled-packet theorem.
The important correction is that the `|14Z cap S|>=7` branch should not be
treated as the active live family: THM-571 already closes the apex-majority
case, modulo the accepted below-frontier input.  The live covering side is:

```text
S = R union 14Q,  1 <= |Q| <= 6.
```

I added:

```text
04-computation/lrc14_few_apex_lift_packet_probe_codex_s152.py
05-knowledge/results/lrc14_few_apex_lift_packet_probe_codex_s152.out
05-knowledge/hypotheses/HYP-2968-lrc14-few-apex-lift-packet-bridge.md
```

The exact packet is the lift `u=14t`.  For each `u`, the original time circle
has fourteen lifts:

```text
t = (u+k)/14,  k=0,...,13.
```

For `14m in 14Q`, the condition is just `||m u||`.  For a residual speed
`r`, lift `k` has exact danger intervals:

```text
||r(u+k)/14|| < 1/14
  <=> u in ((14n-1)/r-k, (14n+1)/r-k) cap [0,1].
```

So this is a labelled Borel/Baire/Haar packet: Q comb, residual deletion
intervals, fourteen lift labels, and exact rational mass.  A positive lift
component is already a strict LRC14 witness.

Stored run:

```text
audited qdiv>14 rows          8190
zero strict lift mass         0
no positive lift              0
smallest lift-safe mass       563/105105
```

Tightest exact-M fallback:

```text
drop(12)->14*6       M=7/89
drop(6)->14*1        M=2/23
drop(6)->14*2        M=2/23
drop(12)->14*12      M=14/173
```

All are above `1/14`.  The minimum open-lift count also rises as the apex
count grows: in this bank it starts at `2` for `k14=1` and reaches `8` for
`k14=5,6`.

The theorem target is now:

```text
primitive qdiv>14, 1<=|14Z cap S|<=6
  -> positive regular-open lift packet
     or K33 / HYP-2908 / THM-572 state lift.
```

This is the F5/F6 boundary-moment bridge in packet form.  It keeps the labels
that raw runner sets and raw apex residue tournaments destroy.
