# LRC14 C27 Unital Block Lift

The q=3 unital does lift the C27 transfer picture, but only branch-locally.

The clean obstruction is the shared H12 hole pair.  In the raw residue-pair
model,

```text
GW H12->D3   = {3, 12, 15, 24}
K33 H12->D9  = {9, 12, 15, 18}
```

These two candidate unital blocks share `{12,15}`.  A q=3 unital is
`2-(28,4,1)`, so two blocks cannot share a pair.  Therefore one unital chart
cannot contain both the tight GW transfer and the K33 near-miss transfer.

That sounds negative, but it is actually useful.  It says the q=3 unital is
not a universal C27 atlas; it is a branch chart.  The tight chart
`GW + 10-petal + 13-petal` embeds.  The K33 chart
`near-miss + 10-petal + 13-petal` embeds.  The S138 two-hole rows are exactly
two-block splices:

```text
drop(10,12)->add(20,24) = petal H10->D7 + GW H12->D3
drop(10,12)->add(20,36) = petal H10->D7 + K33 H12->D9
```

So the unital should be used as a local pair-unique grammar.  If a future proof
wants both `12` branches in one object, it has to split the H12 pair by an
extra branch coordinate or use multiple charts.  Either way, the preserved unit
has changed and must be named.

