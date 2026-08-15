---
id: THM-3378
title: "Projected k=3 z216 gcd24/L129360 row-94 one-high torsion closure"
status: >
  PROVED analytic mechanism + FINITE-EXACT in the declared projected atlas.
  Wall row 94, with E=(1,3,8,10,11,14), has an empty necessary projected
  screen.  The inherited high gate and a positive duplicate-permitting
  two-high gap force exactly one high drift.  All 84 residual passports then
  close by independently verified located differences of effective order two
  or three.  The projected ledger falls 372914 to 372913 and wall rows 110 to
  109; the thirteen-row gcd24/L129360 family remains live, the family count
  stays 12, and the cap stays 216.  This proves no physical entry, arbitrary-k
  statement, rung, or LRC(14).
source: root/repository-archaeology-second-pass-2026-08-14
audit: >
  dependency-pinned queue reconstruction, exact-rational two-high terminal,
  independent direct order-two/three residue search, scalar hostiles, no
  assert or floating-point shortcut, and byte-matching ordinary/optimized
  replays
depends_on:
  - THM-2941-critical-seven-slot-scalar-wall-and-balanced-boundary
  - THM-2979-projected-k3-z275-ten-body-status-and-located-torsion-closure
  - THM-2984-projected-k3-signed-ray-attainment-and-unit-phase-gate
  - THM-3351-projected-k3-z216-sixteen-family-translated-two-high-closure
  - THM-3361-projected-k3-L720720-one-high-translated-residue-closure
related:
  - THM-3320-projected-k3-z216-fourth-ruler-prefix-and-affine-multicover-closure
script: 04-computation/lrc14_j7_k3_z216_gcd24_L129360_row94_one_high_torsion_thm3378.py
output: 05-knowledge/results/lrc14_j7_k3_z216_gcd24_L129360_row94_one_high_torsion_thm3378.out
script_sha256: 62769c8024ccd3f85a71b858e65f909a216e4ea85a7bca2bfb5fd0aa61a43a73
output_sha256: 26db8dc95acac0afe179aedce04543150755860215cd506b55c09cf0c872c37a
semantic_sha256: 436c5160c5e9d3dcaa0f4f3dc4104450670a67582ab5e375c8751dc0ea82c93f
hash_basis: LF-normalized bytes
---

# THM-3378 -- projected row `94` closes by located torsion

**PROVED analytic mechanism + FINITE-EXACT in the declared projected
atlas.**

## 1. Statement and inheritance

Reconstruct the necessary projected `k=3,z1=216` atlas and the intrinsic-cost
queue left by
[THM-3361](THM-3361-projected-k3-L720720-one-high-translated-residue-closure.md).
The next complete family is `gcd24/L129360`, of cost `52,778,880`, on rows

```text
94,161,174,215,237,263,319,354,369,399,430,443,472.       (1)
```

Its first row is

```text
row 94: E=(1,3,8,10,11,14),  L=129360,
        first denominator=5390, high floor=12741,
        30 components, component cost=3880800.             (2)
```

Row `94` has an empty necessary projected screen.  Therefore, inside this
declared projected universe,

```text
projected ledger       372914 -> 372913,
z1=216 wall rows          110 -> 109,
complete families          12 -> 12,
projected k=3 cap         216 unchanged.                  (3)
```

The other twelve rows in `(1)` remain live, with total family cost
`48,898,080`.  Intrinsic cost orders the queue; it is not an exclusion
inequality.

## 2. Exact inherited screen

The exact ray/common-status screen for row `94` partitions as

```text
571 = 292 crude capacity exclusions
      +195 exact common-status exclusions
      + 84 residual passports.                            (4)
```

The full-table status route supplies all `195` status witnesses.  Every
residual denominator tuple contains `5390` exactly once and has lcm `129360`.
The companion pins the queue, the row, its component count, all inherited
source/output hashes, and the digest

```text
6f1137bbe6f0541217d4c3854fd94fcc9bc5eaf106ba3510e7bf06edfae86f9f (5)
```

of the ordered residual packet.

## 3. Exactly one high drift

THM-2941's strict scalar wall forces the largest later drift to satisfy

```text
M > 13L/132 = 12740.                                    (6)
```

Since the first drift is `216<12741`, every projected completion has at least
one high later drift.  Exact scalar arithmetic alone still admits `78`
zero-high hostile records, so this inherited strict gate is load-bearing.

Across all `84` residual passports, the duplicate-permitting two-high
exclusion has the exact minimum gap

```text
268369544 / 94245446295 > 0.                            (7)
```

Thus no completion has two high drifts, even when both chosen high labels may
use the same denominator.  Equations `(6)--(7)` force exactly one high drift.
Enumerating its label and the two low labels gives exactly `84` one-high
cases, one per residual passport, in five low-pair classes.

## 4. Located order-two/three phase terminal

Fix a one-high case.  Let `d|L` be its high denominator and let `C` be the
complete cells fixed-safe for the first drift and both low labels.  THM-2984
identifies a concrete translated high drift in the normal form

```text
z=(L/d)u+hL,                  gcd(u,d)=1.                (8)
```

For two complete cells `c_1,c_2`, their phase difference is
`u(c_2-c_1)/d mod 1`: the height `h` cancels, and multiplication by the unit
`u` preserves exact additive order.  The cells dangerous at a fixed local
coordinate form a translated strict-open cyclic interval of length `1/7`.

For every one of the `84` cases, an exact residue search finds two cells in
`C` whose residue difference has effective additive order two or three in
`Z/dZ`.  Their phase separation is therefore at least `1/3`, so no translated
strict-open interval of length `1/7` can contain both.  At every local
coordinate at least one of the two cells remains safe.  The projected safe
set is consequently the whole local circle, of mass

```text
1 > 36/91,                                               (9)
```

contradicting THM-2941's necessary completion bound.

The exact census is

```text
effective order 2: 80 cases,
effective order 3:  4 cases.                            (10)
```

The inherited located-torsion engine and a second direct implementation of
the modular shift search independently return the same census `(10)`.  The
second implementation constructs the residue map itself, tests shifts `d/2`
and `d/3`, checks the effective order, and rechecks that both selected cells
are fixed-safe for all three fixed labels.  It does not invoke the inherited
torsion-pigeonhole routine.

This identifies exactly what the earlier scalar screen forgot: not more
aggregate mass, but a located phase difference between two lawful cells.

## 5. Scope and exact audit

The conclusion is one row, not the whole family.  Forgetting located phase
leaves the `84` one-high cases; closing row `94` leaves the twelve later rows
in `(1)`.  No entry, owner, physical current, arbitrary `k`, rung, or LRC(14)
claim follows.

The exact-rational companion:

1. reconstructs the post-THM-3361 queue and pins the thirteen-row family;
2. replays all `571` row states and freezes the `84` residual passports;
3. checks the strict high gate, all `78` scalar hostiles, the positive gap
   `(7)`, and all `84` one-high cases;
4. checks both the inherited and independent located-torsion witnesses `(10)`;
5. verifies that the inherited terminal returns no unresolved row-94 case;
6. uses no floating-point arithmetic or `assert`-dependent correctness; and
7. pins every imported computational dependency by LF-normalized SHA-256 and
   checks THM-3361's queue-boundary tuple explicitly.

Ordinary and optimized replays byte-match the stored output.  Reproduce with

```text
python3 04-computation/lrc14_j7_k3_z216_gcd24_L129360_row94_one_high_torsion_thm3378.py
python3 -O 04-computation/lrc14_j7_k3_z216_gcd24_L129360_row94_one_high_torsion_thm3378.py
```

The artifact hashes and semantic digest are pinned in the frontmatter.  Hence
the necessary projected screen of row `94` is empty, and `(3)` follows.

**QED.**
