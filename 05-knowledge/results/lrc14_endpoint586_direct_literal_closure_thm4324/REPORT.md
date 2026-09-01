# Endpoint-586 unchanged-carrier literal closure

**FINITE-EXACT / PROVED RELATIVE TO THM-4318. LRC(14) REMAINS OPEN.**

The unchanged 3,925-mask THM-4318 carrier was replayed on all twelve rows of
the post-endpoint-587 endpoint-586 frontier.  The exact scan checks
`12*C(30,9)=171,685,800` labelled pair/body cases.  Eleven rows close
outright.  Only `(50,586)` has carrier misses: 4,090 bodies with body FNV64
`14ce094f4ab4ba94`.  Aggregate exposed bodies are 525,048, hit incidences are
12,125,735, global failure FNV64 is `ffb884b2b17e6ef4`, and pair-ledger FNV64
is `46a2c17caecc55df`.

The exact direct-wall audit evaluates the rigorous lower bound

```text
L(B)=sum_{F & B = 0, popcount(F)<=9} w(F) <= M(B).
```

For `(50,586)`, the wall grid is 26,723,298,545,143,200, with 8,381 open
cells, 5,802 pair-safe cells, 2,424 full failure classes, and 2,371 low-rank
classes.  Every missed body has strict truncated surplus.  The minimum is

```text
63*L(B)-4D = 136,976,666,519,138,544
at B = 011c6405.
```

The full mass is also positive.  Truncated and full masses agree on 3,602
bodies and differ strictly on 488.  O2/O3 carrier and primary outputs are
byte-identical.  The aggregated-class primary and separate raw-cell
implementation agree byte-for-byte on all 4,090 mass rows, detail SHA256
`6dbfc6178c7da5c844abc95e980f2696cf02c3845a83bdfff95f881c12c9c4fd`.

Thus endpoint 586 closes without carrier surgery.  Consuming it gives typed
union 2,229 (FNV64 `035ebf12f02ecc62`) and residual 20,418 (FNV64
`89b73e31224821c4`).  The next endpoint is 585 with twenty-three rows and
FNV64 `8f1b7c8db8fd5e87`.

This is finite-exact for the frozen pool, carrier, typed universe, and layer.
It proves no physical owner/entry map, terminating descent, or LRC(14).
