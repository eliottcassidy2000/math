---
id: THM-1243
title: THE RESONANT TOOTHPICK RAY HAS A UNIFORM GLOBAL REROUTE — the local seven-crack blocker star is punctured by an explicit parity phase
status: RESERVED/PROOF IN PROGRESS. Exact residue formulas indicate a uniform witness for every m>=27, hence in particular the entire m>=42 one-blocker ray of THM-1239; referee and Lean finite arithmetic are being prepared.
source: codex-2026-07-19-S78 continuation
depends_on: [THM-1239]
related: [THM-1236, THM-1240, THM-1242, MISTAKE-185]
---

# THM-1243 — resonant toothpick ray global reroute

For `m>=27`, put `a=7m+1`,

```text
V_m={1,2,3,4} union {a,a+1,...,a+7} union {14a},
q=14m+9=2a+7,
p=3m+4+(m mod 2).
```

The reserved target is the exact all-runner certificate

```text
min_(v in V_m) ||vp/q||
  = floor((3m+5)/2)/(14m+9) > 1/14.
```

Thus THM-1239's single speed `14a` can erase all seven curvature cracks in
one selected `a`-gap for every `m>=42`, but the full thirteen-speed packet has
a uniformly deep lonely phase in another address cell.  The proof is a
parity split and direct modular reduction using `q=2a+7`; no truncated
resonance inverse or equidistribution input is intended.

