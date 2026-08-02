---
id: THM-3114
title: "Projected-k3 z227 exact-screen zero-residual descent"
status: >
  PROVED + VERIFIED-EXACT + HOSTILE-AUDITED.  All thirty rows and all 284
  inherited screen states in the projected k=3, z1=227 layer close without a
  terminal construction.  The projected ledger is 374185 and the cap is
  z1<=226.  No physical-cover or LRC(14) claim is made.
source: root/codex-thm3088-push-2026-08-02
depends_on:
  - THM-3113-projected-k3-z229-terminal-and-z228-screen-double-layer-descent
  - THM-3111-projected-k3-z230-exact-screen-and-compressed-complete-cell-descent
  - THM-2941-critical-seven-slot-scalar-wall-and-balanced-boundary
script: 04-computation/lrc14_j7_k3_z227_exact_screen_zero_residual_descent_thm3114.py
output: 05-knowledge/results/lrc14_j7_k3_z227_exact_screen_zero_residual_descent_thm3114.out
script_sha256: a05dc531b8e32201e0a5dc3ef92bd0f89ecc5c7f25723275dce60d1f98703840
output_sha256: 424bf01ca72a3d5a222eb2411ae58a718bde127109ff0d64182f75342b0c2a45
semantic_sha256: 17347efff69b101d47afc7d7097fcde498d210ed25737c1ef9935d6b3a5baacb
hash_basis: LF-normalized bytes
---

# THM-3114 -- projected-k3 z227 exact-screen zero-residual descent

**PROVED + VERIFIED-EXACT + HOSTILE-AUDITED.**

## 1. Exact layer and screen

In the pinned THM-2941 projected `k=3` necessary atlas, the complete layer is

```text
z1=227: 30 rows = 28 wall + 2 order,
row digest 17c023b6d703e6f2930e2a2606da83e5fe34d870d6067ec202329df210957c6c.
```

Recomputing the full THM-3113/THM-3111 screen gives

```text
284 states = 143 crude + 141 exact-status + 0 residual.
```

Both order rows close crudely (`9/9` states).  The 141 status closures use the
inherited exact Farkas/status route; the direct normalization is retained as
a hostile control and closes no additional state at this layer.  The complete
screen-record digest is

```text
a86dc48861dbbf23d86cea24201a1b9634c118eec6e124effb5fa4e6aea688df.
```

The ray quotient and exact-status screen enlarges the universe of actual
projected assignments.  Its empty residual bank therefore closes every
actual assignment in the layer; no terminal or carrier implication is needed.

## 2. Ledger, next boundary, and scope

Removing the thirty rows gives

```text
374215 - 30 = 374185,                 z1 <= 226.
```

The next layer is occupied:

```text
z1=226: 13 rows = 13 wall + 0 order,
row digest ec2f1218d61dd69d3811d68eb87a23254d276bb75647be2d4c883affe53a520e.
```

That layer is not claimed here.  This theorem concerns only the maintained
projected `k=3` necessary atlas.  It neither classifies physical covers nor
changes `k<=1`, the six-body/seven-tail rung, or LRC(14).

Ordinary, optimized, eight-process, and single-process executions are
byte-identical to the stored transcript.  The verifier freezes the complete
row order, every state partition, the exact-status count, the zero residual,
the next-layer census, and the ledger arithmetic with explicit runtime gates.
