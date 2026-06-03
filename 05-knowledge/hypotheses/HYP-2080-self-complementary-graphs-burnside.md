---
name: HYP-2080-self-complementary-graphs-burnside
description: The S565 anti-automorphism Burnside engine extends cleanly to self-complementary UNDIRECTED GRAPHS (A000171) via the edge-complement involution; formula verified exact to n=40 against OEIS (which has n=1..100). New repo capability.
metadata:
  type: hypothesis
  status: CONFIRMED
---

**WHAT:** Extend the S565 orientation-reversal Burnside engine to the parallel family:
self-complementary SIMPLE UNDIRECTED GRAPHS (OEIS A000171), with total count A000088.

This fills the last major gap in the repo's Burnside engine: the four fundamental graph
orientation types are now all covered in one framework.

**KEY DIFFERENCE FROM S565:** For undirected graphs, transposing the two vertices of a
pair {i,j} does NOT flip the edge color (no "iota" from pair-swapping). Modified rules:
- Total iso: every pair-orbit contributes C=2 choices regardless of swap → A000088
- Self-complementary: pair-orbit of length L contributes 2 if L even, 0 if L odd → A000171
  (complement involution 0↔1 applied L times returns to start iff L is even)

Contrast with S565 directed families:
- tournament (C=2, Cfix=0, swap applies iota): total=A000568, SC=A002785
- oriented (C=3, Cfix=1): total=A001174, SC=A005639
- digraph (C=4, Cfix=2): total=A000273, SC=A002499
- **undirected (C=2, no-swap override)**: total=A000088, SC=A000171 (this entry)

**VERIFIED (`self_complementary_graphs_burnside_s566.py`):**
- Brute-force independent enumeration: PASS (n=1..6, 0 mismatches)
- OEIS A000088 cross-check: PASS (n=1..14, all match)
- OEIS A000171 cross-check: PASS (n=1..40, all 0 mismatches vs OEIS b-file n=1..100)
  - Note: OEIS b-file goes to n=100; my computation (n=1..40, 34.9s) is a verification pass
  - Extension to n>100 possible with higher NMAX (roughly p(n)*n² work ≈ 400s for n=100)
- Structural sanity: PASS (sc=0 for n≢0,1 mod 4; (iso-sc) even; sc≤total for all n≤40)

**BUG CAUGHT:** My initial "known" value A000171(13) = 5765760 was WRONG.
The correct OEIS value (confirmed via b-file) is A000171(13) = 5600. The formula agrees.
This is a DATA QUALITY CATCH — do not trust hand-recalled OEIS values; always fetch.

**A000171 nonzero values, n=1..40 (verified against OEIS):**
```
n=1:  1
n=4:  1
n=5:  2
n=8:  10
n=9:  36
n=12: 720
n=13: 5600
n=16: 703760
n=17: 11220000
n=20: 9168331776
n=21: 293293716992
n=24: 1601371799340544
n=25: 102484848265030656
n=28: 3837878966366932639744
n=29: 491247277315343649710080
n=32: 128777257564337108286016980992
n=33: 32966971058719932671168222859264
n=36: 61454877497308462618188532330410573824
n=37: 31464896751148469761776612436741418123264
n=40: 422314689395950135433730499958070655419345928192
```

**HONEST:** verification + repo gap-fill. NOT a new OEIS extension (OEIS has n=1..100).
Value = completing the repo's Burnside framework for all four major graph families,
in a single unified engine. The "swap-override" trick for undirected structures is the
new methodological contribution.

**See:** `04-computation/self_complementary_graphs_burnside_s566.py` (+.out),
`05-knowledge/results/b_graphs_total_s566.txt` (A000088 n=1..40),
`05-knowledge/results/b_graphs_selfcomp_s566.txt` (A000171 n=1..40);
HYP-2078 (S565 directed families), THM-283, HYP-2074, HYP-2064.
