---
id: THM-3113
title: "Projected-k3 z229 terminal and z228 screen double-layer descent"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  All forty-three
  rows in the projected k=3, z1=229 layer and all five rows in the z1=228
  layer are empty.  The projected ledger is 374215 and the cap is z1<=227.
  No physical-cover or LRC(14) claim is made.
source: root/codex-thm3088-push-2026-08-02
depends_on:
  - THM-3111-projected-k3-z230-exact-screen-and-compressed-complete-cell-descent
  - THM-3109-projected-k3-z231-exact-screen-and-complete-cell-cardinality-descent
  - THM-3078-z234-direct-farkas-normalization-and-four-two-high-boundary
  - THM-2941-critical-seven-slot-scalar-wall-and-balanced-boundary
  - THM-1166-seven-wall-fano-gcd-discrepancy
script: 04-computation/lrc14_j7_k3_z229_terminal_z228_screen_double_layer_descent_thm3113.py
output: 05-knowledge/results/lrc14_j7_k3_z229_terminal_z228_screen_double_layer_descent_thm3113.out
script_sha256: 59922de78673b94711d0cfe45c9954b8b2056b85f06e904d690ef18f389d111e
output_sha256: a872fe57e12f878b8e237ba71611a39c1f03ba893d45d3c3ae4dd6480e660595
semantic_sha256: a024972dce2c159096db81a751b49654ebc8fad18e2fb1229dbfa0f01784e7e4
hash_basis: LF-normalized bytes
---

# THM-3113 -- projected-k3 z229 terminal and z228 screen double-layer descent

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

## 1. Two exact atlas layers

In the pinned THM-2941 projected `k=3` necessary atlas, the complete adjacent
layers have

```text
z1=229: 43 rows = 36 wall + 7 order,
z1=228:  5 rows =  5 wall + 0 order.                       (1)
```

Their ordered-tuple digests are respectively

```text
7449dd7ad70cf3c76c32edb2cc509e29989ac008c2e9a9968bceaabf3e2a2587,
ba3819e24debb3708459e0293ca4f35127ac533b1a7491a30212fb3fd9242cda. (2)
```

Fresh THM-3078 screens give

```text
z1=229: 1195 states = 571 crude + 512 status + 112 residual,
z1=228:   99 states =  94 crude +   5 status +   0 residual. (3)
```

The screen-record digests are

```text
a3b822b2ba2413b498dd02b787b41171be51942e71e0392a7cbd53a0ff188059,
71cbcedb29695778b51252c711d9dbfc6ed83e5f886ed382c1fc9c8956e22928. (4)
```

All seven `z1=229` order rows close crudely (`24/24` states); every residual
belongs to a wall row.  The `z1=228` layer reaches zero before any terminal
construction.

## 2. The six `z1=229` terminals

Exactly six wall bodies retain a residual bank.  The complete terminal table
is

| body `E` | masks | positive duplicate-two-high gap | zero-high hostile | coarse / exact |
|---|---:|---:|---:|---:|
| `(1,4,5,9,11,14)` | 14 | `1346009063111/373226766637725` | 11 | `11 / 3` |
| `(1,5,6,8,9,14)` | 1 | `60679854023/18484273855320` | 1 | `1 / 0` |
| `(1,5,9,11,12,14)` | 26 | `145522417903/38298144512628` | 23 | `23 / 3` |
| `(1,8,10,11,12,14)` | 19 | `14971895384/3550561219205` | 16 | `15 / 4` |
| `(2,5,8,9,11,14)` | 42 | `26907523189147/7661147274160380` | 39 | `39 / 3` |
| `(2,8,10,11,12,14)` | 10 | `396680527987/91537658972760` | 10 | `10 / 0` |

The tuple digest is

```text
516b596b690bcb176afbbe02fee672171f2a8bbf359ba031d1983819ad87a9d0. (5)
```

Every displayed gap is strictly positive.  Therefore two high suffix slots
cannot carry the wall demand, while the inherited wall premise requires at
least one high slot.  Every actual residual assignment has exactly one high
slot.  The `100` scalar zero-high passes are hostile relaxations excluded by
that premise, not discarded physical cases.

The high-ray supremum enlarges the actual assignment set.  It produces all
`112` one-high cases, partitioned into `99` coarse-cardinality and `13`
exact-cardinality certificates, with no max-gap fallback and no failure.

## 3. Direct inner carriers and translated support

For each body, all cases have one invariant low-label pair.  Direct enumeration
of every complete grid cell agrees exactly with the scalar and vector builders.
The six carriers contain

```text
34442, 6812, 35738, 25358, 73130, 25798 cells             (6)
```

in the table order above.  Their record digest is

```text
b6ed1106b0061b86909e11559c11244f267fd5112acfe8f550b978f4dccc6c51. (7)
```

The full `112`-case record digest is

```text
e7f794f8ce3ffbfbc014899b7fc7e7b22158ff7eb56d8ad02c93c6b0d7659e33. (8)
```

Every direct complete cell is an inner carrier wholly contained in the strict
open safe set.  For a high denominator `d`, use the translated cap

```text
kappa(d)=ceil(d/7),                                      (9)
```

not a centered-band proxy.  The coarse branches exceed this cap by at least
one.  In each of the thirteen boundary branches the projected support is in
fact all of `Z/dZ`.  The only denominators are

```text
d=2,3,8,9,15,
support-kappa(d)=1,2,6,7,12.                            (10)
```

Thus every case leaves a whole fixed safe cell outside every translated high
band.  THM-2941 completes the carrier and THM-1166 supplies the final measure
contradiction exactly as in THM-3111.

## 4. Ledger and handoff

The two layers contain `43+5=48` rows.  Consequently the projected ledger and
cap become

```text
374263 - 48 = 374215,            z1 <= 227.            (11)
```

The next layer is occupied:

```text
z1=227: 30 rows = 28 wall + 2 order,
row digest 17c023b6d703e6f2930e2a2606da83e5fe34d870d6067ec202329df210957c6c. (12)
```

This is the exact handoff.  The candidate makes no assertion about that layer.

## 5. Scope and evidence

This is a theorem only about the maintained projected `k=3` necessary atlas.
It does not classify physical covers outside that projection, alter `k<=1`,
close the final rung, or prove LRC(14).

The companion recomputes all forty-eight screens, six terminal gaps, all six
direct/scalar/vector carriers, and every translated support.  Ordinary,
optimized, eight-process, and single-process executions byte-match the stored
transcript.  An independent proof-direction audit checked that the screen
enlarges the actual assignment set, the one-high ray bank remains an outer
relaxation, each direct complete cell is an inner safe carrier, and support
strictly above `ceil(d/7)` supplies the contradiction in the asserted
direction.  This is an independent hostile audit, not a second implementation.
