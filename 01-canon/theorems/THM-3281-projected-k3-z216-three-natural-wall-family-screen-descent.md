---
id: THM-3281
title: "Projected-k3 z216 three natural wall-family screen descent"
status: "PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED"
source: root/z216-natural-wall-families/2026-08-03
depends_on:
  - THM-3139-projected-k3-z225-terminal-and-z224-screen-double-layer-descent
  - THM-3264-projected-k3-z216-low-cost-gcd8-seventeen-row-terminal-descent
  - THM-3270-projected-k3-z216-complete-order-row-screen-descent
related:
  - THM-3261-projected-k3-z216-unique-gcd18-terminal-descent
script: 04-computation/lrc14_j7_k3_z216_natural_wall_family_screen_descent_thm3281.py
output: 05-knowledge/results/lrc14_j7_k3_z216_natural_wall_family_screen_descent_thm3281.out
script_sha256: 430dee7ba03e0d5c9ae0df72ac512500de4f7056cb4663d1c8468bfb93a49bfe
output_sha256: b9c23b21bf9766efbbc14aa97e2bb4268ddb6abb09b54a2b7424b8744ff686b2
semantic_sha256: 2330ce22ee3b318da26e70ce2c6c01b020c8f2b78bd958d529640a599514cd32
hash_basis: LF-normalized bytes
audit: >
  An independent hostile audit rebuilt the complete z216 atlas and the three
  intrinsic ruler families, verified all 47 indices and costs, checked empty
  overlap with THM-3264 and THM-3270, reran all 1,299 states, and rechecked
  the contradiction direction of all 445 exact full-table Farkas
  certificates. It separately reproduced every family and union digest, the
  empty residual, and assertion-independent normal/-O/stored parity.
---

# THM-3281 -- projected-k3 z216 three natural wall-family screen descent

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

## 1. Three complete intrinsic families

The pinned projected atlas has

```text
z1=216: 480 rows = 447 wall + 33 order,                 (1)
```

with gcd strata

```text
gcd(216,L): 8^19, 18^1, 24^135, 36^15, 72^310.        (2)
```

After THM-3264 and THM-3270, exactly `429` wall rows remain. Inside them,
select the three complete intrinsic ruler families

```text
(gcd(216,L),L) in {(24,2352),(36,8820),(72,2520)}.     (3)
```

These are whole ruler families, not rows selected by an external work
cutoff. Their exact invoices are

| family | rows | intrinsic cost `sum Lr` | states | crude | exact status |
|---|---:|---:|---:|---:|---:|
| `gcd24/L2352` | 20 | 1,086,624 | 710 | 309 | 401 |
| `gcd36/L8820` | 4 | 829,080 | 179 | 145 | 34 |
| `gcd72/L2520` | 23 | 1,144,080 | 410 | 400 | 10 |
| **union** | **47** | **3,059,784** | **1,299** | **854** | **445** |

The ordered union indices in the original 480-row layer are

```text
6,9,10,14,15,24,25,31,32,40,41,44,45,54,55,68,69,71,
76,77,88,102,103,117,124,125,132,135,177,178,204,242,
246,257,269,270,279,284,285,294,295,322,323,346,375,
402,423.                                                 (4)
```

The packet digest binding row index, body, gcd, ruler, component count,
intrinsic cost, and wall status is

```text
967d9ef712a12d081394837460abff934e0cc81f7a356e1f6da6bd677d077560. (5)
```

## 2. Exact family screens

The complete THM-3139 upper-relaxation screen gives

```text
1299 states = 854 crude + 445 exact-status + 0 residual,
direct Farkas=0, inherited full-table exact Farkas=445. (6)
```

The family packet and screen digests are

| family | packet SHA256 | screen SHA256 |
|---|---|---|
| `gcd24/L2352` | `f4101ec3457957d91aef1fefcb2771cb08f400d8b2252e8a680519f192e252e7` | `44984dba0ca096cf43eafb6479f1977eb4c7ea9eefbb9c804c1844837008e459` |
| `gcd36/L8820` | `a03cf6fbcf98781424ccc44b0b92c55c9b22625cbe698ea1098997534459e777` | `671d9823cd1ef155cd58d30cdc0a56b7ec4023b6b13d1d9258f1083ec71b2c93` |
| `gcd72/L2520` | `3a1c6a5d45f9fb5ef8c763ff0d32540dd4b03496a99c6f17f19793faa7b5e553` | `ede315ae68c29463aff5f6e02dbcca37e0283c449b5936510bffa4209fba11a0` |

The ordered union screen digest is

```text
05ee24dbeeaec5e67e4a2ae93482453493f1d6edc226d80187c6b32e4acea8de. (7)
```

Its empty residual digest is

```text
2e38e77b22c314a449e91fafed92a43826ac6aa403ae6a8acb6cf58239fbaf5d. (8)
```

The logical direction is one-way and sufficient. The ray quotient screened
by THM-3139 is an upper relaxation of every actual projected packet. Crude
upper exclusions and the 445 verified exact dual contradictions empty that
enlarged state set, so they exclude all 47 rows. No converse from quotient
feasibility is used. The digest binds basis-invariant canonical rows and
branch counts rather than solver-selected dual representatives.

## 3. Disjoint composition and consequence

Every row in `(4)` is a wall row. The set is disjoint from both

```text
THM-3264: 17 previously closed low-cost gcd-8 wall rows,
THM-3270: all 33 order rows.                              (9)
```

Therefore the closures compose without double counting. The projected cap
remains `z1<=216`, while the live ledger and first occupied layer become

```text
373233 - 47 = 373186,
z216: 429 wall + 0 order -> 382 wall + 0 order.          (10)
```

## 4. Scope and exact evidence

This theorem closes exactly the three necessary wall families in `(3)`. It
uses no terminal probe and asserts no divisor action, parity or free-factor
tower, common carrier, owner/phase transport, physical-cover classification,
final rung, or `LRC(14)`.

The assertion-independent companion pins THM-3139, THM-3264, and THM-3270 by
LF source/output/semantic hashes; reconstructs the 480-row layer; binds every
family and union packet; reruns all 1,299 states; verifies every status branch,
the empty residual, disjointness, and ledger arithmetic; and emits a
deterministic transcript. Normal and optimized runs LF-normalize byte-for-byte
to the stored output with empty standard error.
