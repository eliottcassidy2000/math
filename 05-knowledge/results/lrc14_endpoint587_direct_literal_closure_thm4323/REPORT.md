# Endpoint-587 unchanged-carrier literal closure

**FINITE-EXACT / PROVED RELATIVE TO THM-4318. LRC(14) REMAINS OPEN.**

The unchanged 3,925-mask THM-4318 carrier was replayed on all ten rows of the
post-endpoint-588 endpoint-587 frontier.  The exact scan checks
`10*C(30,9)=143,071,500` labelled pair/body cases.  Nine rows close outright.
Only `(50,587)` has carrier misses: 41 bodies, body FNV64
`c76719ced1d5c52b`.  Aggregate exposed bodies are 53,771, hit incidences are
1,587,855, global failure FNV64 is `bae65dc3d3bd34d0`, and the pair-ledger
FNV64 is `0062cc50be726e54`.

The exact direct-wall audit uses the same rigorous truncation as THM-4320 and
the endpoint-588 scratch result:

```text
L(B)=sum_{F & B = 0, popcount(F)<=9} w(F) <= M(B).
```

For `(50,587)`, the wall grid is 53,537,802,887,368,800, with 8,385 open
cells, 5,792 pair-safe cells, 2,420 full failure classes, and 2,365 low-rank
classes.  All 41 missed bodies have strict truncated surplus.  The minimum is

```text
63*L(B)-4D = 283,424,219,270,897,292
at B = 11186405.
```

The complete mass is also strictly positive.  Truncated and full masses agree
on 31 bodies; full mass strictly dominates on ten.  O2/O3 carrier and primary
outputs are byte-identical.  The aggregated-class primary implementation and
independent raw-cell implementation agree byte-for-byte on all masses, detail
SHA256 `5a89b23b60fce1ec4be6d03bdd92c4a3cf3a861afe3d30226f8bfe831154c895`.

Thus the endpoint-587 layer closes without carrier surgery.  Consuming it
gives typed union 2,217 (FNV64 `e6592cbef9b616d8`) and residual 20,430
(FNV64 `4710f750dfcf91ea`).  The next endpoint is 586 with twelve rows and
FNV64 `a1b617faa2e7f63f`.

This is a finite result for the frozen pool, carrier, typed universe, and
endpoint layer.  It proves no physical owner/entry map, terminating descent,
or LRC(14).
