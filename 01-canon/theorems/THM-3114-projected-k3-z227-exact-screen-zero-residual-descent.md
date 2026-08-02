---
id: THM-3114
title: "Projected-k3 z227 screen and z226 terminal double-layer descent"
status: >
  PROVED + VERIFIED-EXACT + HOSTILE-AUDITED.  All thirty z1=227 rows close at
  the inherited screen.  All thirteen z1=226 rows then close through the same
  screen plus three exact one-high complete-cell terminals.  The projected
  ledger is 374172 and the cap is z1<=225.  No physical-cover or LRC(14)
  claim is made.
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
addendum_script: 04-computation/lrc14_j7_k3_z226_terminal_descent_thm3114.py
addendum_output: 05-knowledge/results/lrc14_j7_k3_z226_terminal_descent_thm3114.out
addendum_script_sha256: cebf964f1e7facf1de78c6afcc58bf56b620d49a01fc3788c7667d164e27c997
addendum_output_sha256: e5b2dc822dc40f5f351078f4dd54732fbea15a589a08ab64af0312017b46c1d3
addendum_semantic_sha256: cdd4f83bbd3a3244244bf14baad8935a43a4a33c0bc973acff5dcc4048a4dc52
hash_basis: LF-normalized bytes
---

# THM-3114 -- projected-k3 z227 screen and z226 terminal descent

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

## 2. The `z1=226` terminal layer

The next complete layer has thirteen wall rows.  Its inherited screen gives

```text
979 states = 526 crude + 363 exact-status + 90 residual.
```

The ninety masks occur on exactly three bodies, with residual counts

```text
(1,5,6,9,12,14):   1,
(1,5,9,11,12,14): 65,
(1,9,10,11,12,14):24.
```

Their duplicate-permitting two-high gaps are respectively

```text
2436193/902641740,
15491357285/4724552761929,
1545772719847/373939773753546,
```

all strictly positive.  The wall premise requires a high suffix slot, while
each positive gap forbids two; every actual residual assignment therefore
has exactly one high slot.  The high-ray supremum enlarges the actual bank
and yields `90` one-high cases.  Of these, `84` exceed the translated
`ceil(d/7)` cap by coarse cardinality.  The remaining six have full exact
support in `Z/dZ`.

For the three invariant low-label pairs, direct enumeration agrees with both
the scalar and vector complete-cell builders.  The inner carriers have

```text
3782, 38756, 39320
```

strict safe cells.  The minimum coarse margin and the minimum exact support
slack are both one.  Thus every case leaves a fixed complete safe cell outside
the translated high band, and the THM-3113 carrier contradiction applies.

## 3. Ledger, next boundary, and scope

Removing both layers gives

```text
374215 - 30 - 13 = 374172,            z1 <= 225.
```

The next layer is occupied:

```text
z1=225: 78 rows = 78 wall + 0 order,
row digest 9a4869fbb3886b058bdb2b60209b1fc88ab4c82ccbc7c56eba42b2d118ee8719.
```

That layer is not claimed here.  This theorem concerns only the maintained
projected `k=3` necessary atlas.  It neither classifies physical covers nor
changes `k<=1`, the six-body/seven-tail rung, or LRC(14).

For both layers, ordinary, optimized, eight-process, and single-process
executions are byte-identical to their stored transcripts.  The addendum
freezes all thirteen screens, three positive gaps, ninety cases, three
direct/scalar/vector carriers, every translated support, the next-layer
census, and the ledger arithmetic with explicit runtime gates.
