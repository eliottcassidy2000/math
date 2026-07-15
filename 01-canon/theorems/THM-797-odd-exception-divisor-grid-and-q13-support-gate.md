---
id: THM-797
title: Odd-exception divisor grids and the q=13 folded-support gate
status: CLAIM STUB — general residue-shell lemma and q=13 corollaries proved; exact replay and full integration in progress
source: codex-2026-07-14-S10 global-erosion selection analysis
depends_on:
  - THM-772   # two-sheet packet arithmetic and mandatory deep exception
  - THM-789   # global erosion target
related: [THM-769, THM-774, HYP-6800, HYP-6820]
---

# THM-797 — Odd-exception divisor grids and the q=13 support gate

> **Namespace claim stub.**  An odd divisor of one exception turns deep
> unit-grid points into an explicit acceptance shell for the other exception.
> At `q=13` this proves a new uniform support gate.  The complete proof,
> verifier, and sharp-survivor audit are being integrated now.

In THM-772's two-sheet packet `2U union {x,y}`, let odd `q|x`.  The folded
unit classes `p` for which

```text
least_absolute_residue_q(up) >= ceil(q/11)  for every u in U
```

give points `p/q in E_U`.  If `yp=qN+s` is balanced, then

```text
p/q lies in the folded diamond
iff s is odd and 1<=|s|<=floor(2q/13).
```

Thus failure of containment in this explicit `y`-shell proves the global
erosion escape required by THM-789.

For `q=13`, let `C=(Z/13)^*/{+/-1}` and let `S(U)` be the folded residue
support of `U`.  If exactly `x` is 13-divisible, hypothetical containment
forces

```text
S(U)=C  or  S(U)=C\{+/-y}.
```

If both exceptions are 13-divisible, it forces `S(U)=C`.  Hence all support
at most four cores, every misaligned five-class core, and every non-full core
in the double-13 branch are eliminated uniformly.
