# Endpoint-588 unchanged-carrier literal closure

**FINITE-EXACT / PROVED RELATIVE TO THM-4318. LRC(14) REMAINS OPEN.**

This scratch packet audits the sixty-six rows of the post-THM-4320 typed
endpoint-588 frontier.  It retains the 3,925-mask THM-4318 carrier unchanged.
The exact all-body replay checks

```text
66 * C(30,9) = 944,271,900
```

labelled pair/body cases.  Fifty-six rows close by the carrier.  The carrier
misses 144,867 bodies on ten rows:

| q | failures | body FNV64 |
|---:|---:|:---|
| 25 | 9 | `8fe6c7431c54c6ff` |
| 50 | 99,836 | `4a8d6e89f7aa5d75` |
| 96 | 1,130 | `5766a653739cf191` |
| 100 | 14 | `f0fdda5ad80a46b7` |
| 105 | 60 | `a6f76a994fd08192` |
| 206 | 7 | `def96e7fae5db299` |
| 210 | 43,799 | `f756dc445e2790cd` |
| 256 | 7 | `392c4d67394b2e17` |
| 420 | 1 | `664e6def29193bd0` |
| 462 | 4 | `5f993deacddbf9d3` |

The aggregate failure FNV64 is `6f51fa88f3b09cdc`; aggregate exposed bodies
are 5,001,257 and carrier-hit incidences are 76,551,030.

For each failed body `B`, the direct audit partitions the circle at every
wall of the pool and pair.  On pair-safe cells let `F` be the pool-failure
mask and let `w(F)` be the cell width on the exact integer grid.  It evaluates

```text
L_q(B) = sum_{F & B = 0, popcount(F) <= 9} w(F)
M_q(B) = sum_{F & B = 0} w(F).
```

Since every omitted width is nonnegative, `L_q(B) <= M_q(B)`.  Every one of
the 144,867 carrier failures has strict surplus `63*L_q(B)-4*D_q > 0`:

| q | minimum truncated surplus | minimizing body |
|---:|---:|:---|
| 25 | 3,845,468,183,436,180 | `0d1cd000` |
| 50 | 3,359,314,066,788,540 | `003c6405` |
| 96 | 1,528,699,623,231,060 | `2d14c001` |
| 100 | 3,874,870,450,697,778 | `36704001` |
| 105 | 731,279,928,922,452 | `0458d082` |
| 206 | 81,246,765,687,627,708 | `14087c01` |
| 210 | 706,240,375,488,648 | `06b82090` |
| 256 | 13,754,370,514,849,146 | `1f106001` |
| 420 | 795,987,916,404,408 | `07584088` |
| 462 | 890,738,628,297,084 | `1e582001` |

Thus every endpoint-588 row satisfies the literal Haar target without carrier
surgery.  O2/O3 carrier outputs are byte-identical.  O2/O3 aggregated-class
audits and a separate raw-cell implementation agree byte-for-byte on all
144,867 truncated and full masses (detail SHA256
`8c2652b9cba7cfdacd4d6e908334220eb08445306fbc09db816d0b0b9c48664d`).

Consuming the layer gives typed union 2,207 (FNV64
`18d067b5614cf47f`) and residual 20,440 (FNV64 `794bd808e92e27cd`).
The next endpoint is 587 with ten rows and FNV64 `f48ca5f1904d6f52`.

This is a finite result for the frozen pool, carrier, typed universe, and
endpoint layer.  It proves no physical owner/entry map, terminating descent,
or LRC(14).
