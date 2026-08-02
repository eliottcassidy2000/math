---
id: THM-3114
title: "Projected-k3 z227 screen and z226 terminal double-layer descent"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  All thirty rows in
  the projected k=3, z1=227 layer and all thirteen rows in the z1=226 layer
  are empty.  The projected ledger is 374172 and the cap is
  z1<=225.  No LRC(14) claim is made.
source: root/codex-thm3088-push-2026-08-02
audit: >
  An independent lower-engine reconstruction verified the 30-row and 13-row
  screens, all 504 exact Farkas certificates, all three positive two-high
  gaps, all 90 one-high masks, all three direct complete-cell carriers, and
  all six exact full-support cases.  It caught two pre-promotion evidence
  defects: raw solver witnesses were removed from semantic hashing under
  MISTAKE-331/333, and three values initially labelled mask hashes were in
  fact row[15] status-instance hashes.  The corrected mask, screen, terminal,
  carrier, and case digests below all agree with the independent build.
  Fresh immutable normal and optimized runs are byte-identical to the stored
  transcript and match all declared LF and semantic hashes.
depends_on:
  - THM-3113-projected-k3-z229-terminal-and-z228-screen-double-layer-descent
  - THM-3111-projected-k3-z230-exact-screen-and-compressed-complete-cell-descent
  - THM-3078-z234-direct-farkas-normalization-and-four-two-high-boundary
  - THM-2941-critical-seven-slot-scalar-wall-and-balanced-boundary
  - THM-1166-seven-wall-fano-gcd-discrepancy
script: 04-computation/lrc14_j7_k3_z227_screen_z226_terminal_double_layer_descent_thm3114.py
output: 05-knowledge/results/lrc14_j7_k3_z227_screen_z226_terminal_double_layer_descent_thm3114.out
script_sha256: 8de1e3d03b5070a84b040ac13a173a598646107f85e7ba0defc2ca070808f162
output_sha256: c6cfd6e9b6bd9a89c2d679d86f9fe545d949269097698f51daa86b30350a543f
semantic_sha256: ad9c2724b6468d586ae70ab120f57519350c7de440c66a3765609bd1a7880d51
hash_basis: LF-normalized bytes
---

# THM-3114 -- projected-k3 z227 screen and z226 terminal double-layer descent

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

## 1. Exact adjacent layers

In the pinned THM-2941 projected `k=3` necessary atlas,

```text
z1=227: 30 rows = 28 wall + 2 order,
z1=226: 13 rows = 13 wall + 0 order.                       (1)
```

The ordered-row digests are respectively

```text
17c023b6d703e6f2930e2a2606da83e5fe34d870d6067ec202329df210957c6c,
ec2f1218d61dd69d3811d68eb87a23254d276bb75647be2d4c883affe53a520e. (2)
```

Fresh exact screens give

```text
z1=227: 284 states = 143 crude + 141 status +  0 residual,
z1=226: 979 states = 526 crude + 363 status + 90 residual. (3)
```

The certificate-free screen-record digests are

```text
b4da9b80d12c969cf5561a0e30e36d282a955e77ff32a5f88fb9514cf78fb414,
5a64104b1ca45e7dd4fe08cfc5e0bcc5a3226a7f05740dc20f361d664a2d98a7. (4)
```

All two `z1=227` order rows close crudely (`9/9` states).  Thus the complete
`z1=227` layer dies at the screen.  The `z1=226` residual occupies exactly
three wall bodies.

## 2. The three z226 terminals

The complete residual table is

| body `E` | masks | positive duplicate-two-high gap | zero-high hostile | coarse / exact |
|---|---:|---:|---:|---:|
| `(1,5,6,9,12,14)` | 1 | `2436193/902641740` | 0 | `0 / 1` |
| `(1,5,9,11,12,14)` | 65 | `15491357285/4724552761929` | 61 | `61 / 4` |
| `(1,9,10,11,12,14)` | 24 | `1545772719847/373939773753546` | 23 | `23 / 1` |

The true residual-mask digests `sha256(repr(row[13]))` are

```text
85ba36a1535644d7d51164a512c8db1bf3ea76bc5f1ad47f5c6f6cfe7397b22a,
997cc17a4be38886e37970060f6c0d92ee6904b345e973083f9b1b52096ad882,
fe196c407eb6d7d1fc5a65f283a8ce078e534602152521e8b1081d908ed725ea. (5)
```

Every displayed gap is strictly positive.  Together with the inherited wall
gate requiring at least one high suffix slot, this forces exactly one high
slot in every actual residual assignment.  The `84` scalar zero-high passes
are hostile relaxations excluded by that gate.  All `90` one-high masks close:

```text
90 = 84 coarse cardinality + 6 exact cardinality,          (6)
```

with no max-gap fallback and no failure.  The terminal tuple digest is

```text
f343c72ad5b49b37eca86171374b274be6ef48bd9e28115107daf9ccb140c5a3. (7)
```

## 3. Direct complete-cell audit

An independent scalar and NumPy reconstruction agrees on the three carriers:

| body | low labels | cells | cell digest |
|---|---|---:|---|
| `(1,5,6,9,12,14)` | `(234,243)` | 3,782 | `b336c32dfe1bb0e3cf7f9ca3e2be9ef7b1e9017cb95747d6e419d2ec51018b5d` |
| `(1,5,9,11,12,14)` | `(234,243)` | 38,756 | `f471b80160444ce50bdee83b9ca6fe2618e6c6d8f7ea9d034688db72dc3e1812` |
| `(1,9,10,11,12,14)` | `(234,260)` | 39,320 | `190facf1ee0fdbbb650380de0ae914fb6af7553a365583d4bc5ae87da964480e` |

Every direct complete cell is an inner carrier wholly contained in the strict
open safe set.  For a high denominator `d`, the translated cap is

```text
kappa(d)=ceil(d/7).                                      (8)
```

The coarse cases exceed it by at least one.  The six exact cases have

```text
d=(2,2,3,9,15,8),
support-kappa(d)=(1,1,2,7,12,6),                        (9)
```

and every support is all of `Z/dZ`.  Hence one whole fixed safe cell survives
every translated high band; THM-2941 completes the carrier and THM-1166 gives
the contradiction.  The carrier and case digests are

```text
f5060567a3cc0341902c541ec6910161d509b48d16f5303a467095c40fbab66f,
6eaf93061816e63631a0e5180697cc1be7fe69980705d2d87f0a4a822698cfa9. (10)
```

## 4. Ledger, evidence, and scope

The two layers contain `30+13=43` rows, so promotion would give

```text
374215 - 43 = 374172,             z1 <= 225.             (11)
```

The next layer is occupied:

```text
z1=225: 78 rows = 78 wall + 0 order,
row digest 9a4869fbb3886b058bdb2b60209b1fc88ab4c82ccbc7c56eba42b2d118ee8719. (12)
```

Every returned Farkas certificate is rebuilt and checked over exact rationals.
Under MISTAKE-331/333, however, the persisted screen and semantic hashes bind
only `tuple(row[:19])` plus the basis-invariant direct/legacy counts; they never
bind the solver-selected dual representative or contradiction magnitude.

This is only a theorem about the maintained projected `k=3` necessary atlas.
It does not classify physical covers outside that projection, alter `k<=1`,
close the final rung, or prove LRC(14).
