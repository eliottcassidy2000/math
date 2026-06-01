---
id: HYP-2036
status: STUB
source: codex-2026-06-01-S546
related:
  - HYP-2017
  - HYP-2019
  - HYP-2029
  - HYP-2031
  - HYP-2032
  - HYP-2033
  - HYP-2035
  - THM-369
---

# HYP-2036: p-adic zero-branch covers form a tree trienerment for LRC sieve channels

**Claim reserved.** This session will build the p-adic tree layer requested by
the user after the Gabor/t(r)ienerment work.  The intended object is the
prime-power zero-branch tree: for every prime power `q=p^d <= n`, the node
records how many speeds fall in the observer residue ball `0 mod q`.  An empty
zero branch gives the THM-369 sieve witness `t=1/q`; a covered zero branch is a
local obstruction/debt carrier.

**Why it is being claimed now.** HYP-2032 already identifies the p-adic
ultrametric as the cleanest arithmetic metric and asks for the full p-adic
tree for composite `n*`.  This HYP-2034 slot is reserved for the concrete
tree-cover core computation, its Tournament Analysis fingerprints, and its
connection to Gabor zero columns/event alphabets.  It is intentionally narrower
than the incoming S544 p-adic integration and HYP-2035 channel-rank statement:
this session measures which runners serve as singleton carriers for
prime-power zero branches and how those cover cores behave across fixed speed
families.  Evidence is still missing until
`lrc_padic_tree_trienerment_s546.py` is implemented and run.

**Files planned.** `04-computation/lrc_padic_tree_trienerment_s546.py`;
`05-knowledge/results/lrc_padic_tree_trienerment_s546.out`;
`07-reflections/lrc-padic-tree-cover-trienerment-s546.md`.
