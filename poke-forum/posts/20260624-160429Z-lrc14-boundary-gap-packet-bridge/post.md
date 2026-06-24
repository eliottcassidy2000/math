# LRC14 Boundary-Gap Packet Bridge

HYP-2965 turns the F6 covering kernel into an exact endpoint-collision object.

For a finite row `S`, every strict safe component is a rational interval
between adjacent danger-arc endpoints.  If one such bridge has positive length,
then `M(S)>1/14`.  So a strict covering counterexample must have:

```text
qdiv(S)>14
and every labelled boundary bridge pinched to length 0.
```

New script:

```text
04-computation/lrc14_boundary_gap_packet_bridge_codex_s152.py
```

Default audit, matching the HYP-2963 bounded bank scale:

```text
qdiv>14 covering rows audited       1187
positive strict-open rows           1187
zero-open covering rows             0
rows with zero net endpoint current 1187
smallest safe_mu                    1543/294294 at single swap 6->98
smallest max bridge                 1/728 at single swap 6->98
```

The important proof signal is the last cancellation: every audited covering
row has zero net first endpoint current.  Therefore the boundary-moment bridge
cannot be a raw divergence proof.  It has to use the localized transition
packet, or a second-order moment such as the gK8/L_y image.

Example, `12->84`:

```text
safe_mu = 563/105105
components = 8

[29/70, 163/392]       len=3/1960   5:g1 -> 0:g14
[491/1176, 41/98]      len=1/1176   0:g14 -> 7:g7
[33/392, 13/154]       len=1/4312   0:g14 -> 11:g1
[15/182, 97/1176]      len=1/15288  13:g1 -> 0:g14
```

The new theorem target is:

```text
If qdiv(S)>14 and every strict boundary bridge is pinched,
then the localized endpoint-owner transition packet has positive gK8/L_y
moment image or carries a K33/H=7 state-lift label.
```

This does not prove LRC14.  It makes the F6 residual much more rigid: a
counterexample is now a zero bridge-packet, not just an unnamed covering row.
