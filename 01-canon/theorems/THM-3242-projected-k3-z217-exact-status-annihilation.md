---
id: THM-3242
title: "Projected-k3 z217 exact-status annihilation"
status: "PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED"
source: root/projected-k3-z217/2026-08-03
depends_on:
  - THM-3109-projected-k3-z231-exact-screen-and-complete-cell-cardinality-descent
  - THM-3230-projected-k3-z219-common-gcd-three-terminal-descent-and-cap218
related:
  - THM-2941-critical-seven-slot-scalar-wall-and-balanced-boundary
  - THM-2596-modular-free-factor-farey-gram-owner-cocycle
script: 04-computation/lrc14_j7_k3_z217_exact_status_annihilation_thm3242.py
output: 05-knowledge/results/lrc14_j7_k3_z217_exact_status_annihilation_thm3242.out
hash_basis: LF-normalized bytes
script_sha256: 165f89184856aae9c4d784e060123c3b639d232d135ff3c744cdfefb5580454a
output_sha256: 10741cdf762edb9b370f91bbfe8739208d349d1d2cbaa37c5b68f3f1a37d34f3
semantic_sha256: a28872885af6aef23f3aa00fe026b5f4c23b6d48302924fab37279cc2a4daf50
audit: >
  An independent referee reconstructed the eight-row atlas layer, replayed all
  sixty-six raw Farkas witnesses under normal and optimized interpreters, and
  reproduced every census, modulus guard, screen digest, and empty-residual
  digest. The canonical companion sorts worker results by atlas index; its
  normal and optimized transcripts byte-match the stored output.
---

# THM-3242 -- projected-k3 z217 exact-status annihilation

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

## 1. The isolated layer

In the pinned THM-2941 projected `k=3` necessary atlas, the complete level is

```text
z1=217: 8 rows = 7 wall + 1 order,                         (1)
row SHA256 885fa20869364ee4a808324612ad373278464a339bf47e183418c1731af4ea6d.
```

For every row body `E`, with `L=14 lcm(E)`, one has

```text
gcd(217,L)=7,                    first_d=L/7.              (2)
```

Thus the ambient pointed divisor interval is isomorphic to `Div(7)`.  This is
only modulus metadata: it supplies neither seven occupied states nor an
invertible `C7` action.

## 2. Exact status annihilation

The complete exact screen gives

```text
total: 66 = 0 crude + 66 exact-status + 0 residual,
wall:  64 = 0 crude + 64 exact-status + 0 residual,
order:  2 = 0 crude +  2 exact-status + 0 residual.        (3)
```

All sixty-six exclusions use the inherited raw-verified rational Farkas path;
none uses the direct shortcut.  The canonical nineteen-field screen digest is

```text
69e3333a4b1dfbbcdd8530d933aa787b7a9b295f01308fcb613f3e65e8eed11b. (4)
```

The rowwise state counts are

```text
(1,2,6,9,10,12):16,  (1,3,4,8,10,12):2,
(1,3,4,9,10,12):10, (1,3,6,9,10,12):20,
(1,4,6,9,10,12):13, (1,6,7,9,12,14):3,
(2,4,6,9,10,12):1,  (2,6,7,9,12,14):1.                 (5)
```

Every count in `(5)` is killed by exact status.  Consequently the residual
body and mask banks are empty, with canonical empty digest

```text
2e38e77b22c314a449e91fafed92a43826ac6aa403ae6a8acb6cf58239fbaf5d. (6)
```

No terminal complete-cell step is needed.

## 3. Exact consequence and the intervening-layer boundary

Equations `(1)--(6)` prove that the entire projected `k=3`, `z1=217` layer is
empty.  This is an intentionally out-of-order closure.  THM-3230 leaves the
current cap at `218`, and the intervening `z1=218` layer remains open.  Hence
this theorem makes **no present cap or ledger change**:

```text
current projected cap = 218,              current ledger = 373411. (7)
```

If a later theorem discharges `z1=218`, this already-proved empty layer may be
used in the same contiguous descent; it must not be silently counted before
that event.

## 4. Evidence and scope

The companion reads the promoted THM-3109 engine and its pinned atlas, derives
all eight tasks afresh, enforces `(2)` on every task and returned screen row,
replays every Farkas certificate, canonicalizes parallel completion order by
atlas index, and checks `(1)--(6)`.  Normal, optimized, and stored transcripts
are byte-identical.  An independent script performed the same reconstruction
from separately frozen normal and optimized checkpoints and then reran every
worker directly.

This theorem proves only one projected necessary-atlas layer empty.  It proves
no common carrier, owner, ancestry, physical transition, `C7` torsor, modular
free-factor action, projected `k<=1` statement, final-rung statement, or
LRC(14) conclusion.
