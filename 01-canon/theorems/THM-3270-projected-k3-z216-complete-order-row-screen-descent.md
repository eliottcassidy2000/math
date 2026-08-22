---
id: THM-3270
title: "Projected-k3 z216 complete order-row screen descent"
status: "PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED"
source: root/z216-order-row-screen/2026-08-03
depends_on:
  - THM-3139-projected-k3-z225-terminal-and-z224-screen-double-layer-descent
  - THM-3264-projected-k3-z216-low-cost-gcd8-seventeen-row-terminal-descent
related:
  - THM-3261-projected-k3-z216-unique-gcd18-terminal-descent
script: 04-computation/lrc14_j7_k3_z216_order_row_screen_descent_thm3270.py
output: 05-knowledge/results/lrc14_j7_k3_z216_order_row_screen_descent_thm3270.out
script_sha256: dbcc60e3d691483e023486cdb9b935381e5172cea2b828b51ed3da99560fe2ab
output_sha256: 0fa283cc68ac831579e17b4e9d817cd5dd13a4c58363d51335e0c1911897d1f0
semantic_sha256: a8eb670dbb06bedc9c9180dc23041a77cb5c5f717f1d2ba2b1de52f8322bc602
hash_basis: LF-normalized bytes
audit: >
  An independent hostile audit reconstructed the complete 33-row order
  sublayer and all four intrinsic cost families, checked empty overlap with
  THM-3264, reran all 481 screened states and the 18 exact Farkas
  certificates with their contradiction directions, verified the empty
  residual bank and ledger arithmetic, and reproduced assertion-independent
  normal/-O/stored byte parity with empty stderr.
---

# THM-3270 -- projected-k3 z216 complete order-row screen descent

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

## 1. The complete order sublayer

After THM-3264, the occupied `z1=216` layer contains

```text
462 rows = 429 wall + 33 order.                             (1)
```

Fresh reconstruction from the pinned 6,060-row projected atlas gives the
complete order-row index set

```text
0,1,2,4,12,13,16,20,22,29,48,51,70,73,74,85,86,
99,100,114,128,130,172,241,244,253,255,267,287,289,
291,317,397.                                                (2)
```

With the inherited exact-work invoice `C(E)=L(E)r(E)`, the four intrinsic
families are

```text
gcd(216,L)  L     rows       total C
24          1176     1         25,872
24          1680    23        836,640
72           504     1          9,072
72          1008     8        157,248
-------------------------------------
total               33      1,028,832.                     (3)
```

Equivalently, gcd `24` contributes `24` rows of total cost `862,512`, and
gcd `72` contributes `9` rows of total cost `166,320`. The ordered packet
digest binding indices, bodies, gcds, rulers, component counts, and costs is

```text
ce13ef797ac9962b3a5cee0fa41450ccb6a2545e5df03463c375b620080c6060. (4)
```

## 2. Disjointness from the inherited wall closures

THM-3264 closes exactly the `17` wall indices

```text
8,23,39,57,64,66,115,142,150,152,197,272,277,300,304,306,338. (5)
```

The companion rederives `(5)` from `gcd(216,L)=8` and the theorem's pinned
cost cutoff, checks that every row in `(5)` is a wall row, checks that every
row in `(2)` is an order row, and verifies empty intersection. Thus this is
not double counting a prior closure.

## 3. Complete exact screen

The complete THM-3139 screen on all `33` rows gives

```text
481 states = 463 crude + 18 exact-status + 0 residual,
direct Farkas=0, inherited exact Farkas=18.                 (6)
```

The ordered canonical-screen digest is

```text
f35dcbe7658ea48df9cfe08c18602ec17242afaeacfcdc88f21d310bfd328ced. (7)
```

The residual bank is empty, with digest

```text
2e38e77b22c314a449e91fafed92a43826ac6aa403ae6a8acb6cf58239fbaf5d. (8)
```

The inherited screen is an upper relaxation of the projected row conditions.
Consequently `(6)` excludes every order row. Since no state reaches the
residual bank, no terminal argument is used.

## 4. Consequence and stopping boundary

Composing only after THM-3264 gives

```text
373266 - 33 = 373233.                                      (9)
```

The projected cap remains `z1<=216`, and the occupied layer becomes

```text
429 rows = 429 wall + 0 order.                             (10)
```

This theorem proves only the complete `z216` order sublayer. It asserts no
terminal probe, divisor action, free-factor or parity tower, common carrier,
owner or phase transport, physical-cover classification, final rung, or
LRC(14).

## 5. Exact evidence

The assertion-independent companion pins THM-3139 and THM-3264 by LF source,
output, and semantic hashes; reconstructs the complete `480=447+33` layer and
all gcd strata; binds `(2)--(5)`; reruns all `481` screened states; checks the
empty residual and ledger arithmetic; and prints a deterministic transcript.
Normal and optimized runs LF-normalize byte-for-byte to the stored output with
empty standard error.
