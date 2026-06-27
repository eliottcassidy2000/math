# LRC14 NORK Pinch-Template Atlas

HYP-2966 attacks the remaining HYP-2956 F6 bucket as NORK:

```text
No Open Residual Kernel
= qdiv >= 14,
  no strict safe open interval,
  and not AP/Goddyn-Wong.
```

New computation:

```text
04-computation/lrc14_nork_pinch_template_audit_codex_20260624.py
05-knowledge/results/lrc14_nork_pinch_template_audit_codex_20260624.out
```

Default bank:

```text
AP
one-swap add<=420
two-swap add<=60
three-swap add<=34
four-swap add<=24
```

Readout:

```text
generated rows          705940
exact qdiv>=14 rows     141351
non-AP/GW F6/NORK rows  0

F1 AP/GW boundary       2
F2 unit-petal positive  28762
F3 K33 positive         340
F4 q14-front positive   78651
F5 covering positive    33596
```

So in this enlarged AP-neighborhood atlas, AP and GW are still the only
zero-open threshold-support rows.  Everything else creates a positive open
front.

The useful move is not the count.  The useful move is the **pinch template**:
for each positive row, keep the shortest strict safe interval, exact endpoints,
endpoint owners, width, slack, q-class, C27-normalized owner labels, and atom
keys.  Scalar mass is too blunt; the endpoint owners say why the mass cannot
collapse to a hidden F6 kernel.

Recurring pinch shapes:

```text
13L -> 12R    q14/off-apex
14L -> 13R    covering
11L -> 16R    q14/four-swap
7L  -> 20R    unit-petal
5L  -> 36R    K33
```

This suggests the next theorem:

```text
NORK pinch theorem:
Every primitive non-AP/GW AP-source-core packet creates a positive
endpoint-owner pinch template, unless it constructs HYP-2908/THM-572
or a genuinely new F7 Johnson-harmonic sector.
```

Tournament Analysis was run on proof carriers and pinch templates, not runners.
The retention tournament is transitive: preserve F6/NORK status and boundary
skeletons first, then pinch templates, fixed-margin packets, C27/K33 labels,
qdiv gates, scalar mass, and only last the raw runner set.

Creative synthesis note: this is the same moral split as the fixed-margin
swap-chain paper's count sector versus Johnson-harmonic sectors.  The LRC
scalar sector is `qdiv/M/Haar`; the non-scalar sector is owner-labelled:
boundary endpoints, C27/unital address, K33/state-lift flag, and
boundary-moment route.  If the proof closes, it closes in the labelled sector
first and scalarizes only after discharge.
