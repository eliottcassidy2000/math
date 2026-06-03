---
id: HYP-2081
status: SPECULATION + grounded core — worry/ignore separation & prime-gear clock grounded; covering-system bridge speculative
source: opus-2026-06-02-S564
related:
  - HYP-2080
  - HYP-2075
  - HYP-2063
  - THM-369
---

# HYP-2081: which clocks matter — worry/ignore split, the prime-gear n-clock, and the covering-system bridge

**Grounded:**
- WORRY vs IGNORE: IGNORE = positive safe-measure (orbit spends an interval in B, hits trivially) = ~all sets (300/302 random; ALL incommensurate; ALL low-resonance); certified by ~30 random clocks (each safe w.p. ~(1-2/n)^13≈0.135). WORRY = measure-0 = resonance-maximal tight family (2/302 = AP, V*). Only worry-cases can fail LRC.
- CLOCKS: complete family = O(k²) PAIR-SUM clocks t=m/(v_a+v_b) (S557/S562); for the WORRY set, ONE clock = the n-clock t=j/n (tight ⟹ binding pair sums to n, S553/S556).
- PRIME-GEAR MODEL (n=14): CRT ℤ/14=ℤ/2×ℤ/7 ⇒ the n-clock is two coupled gears (parity × mod-7); a runner is at the observer iff it FULLY aligns (n | vj). AP has no multiple of 14 ⇒ no runner ever fully aligns ⇒ every gear clear. LRC@14 = can the two gears always be co-set so no runner fully aligns? Difficulty = the gear-coupling (apex/parity, S559).

**Speculative remodels (menu):** (1) LRC ⇄ COVERING SYSTEMS — counterexample = danger arcs exactly cover; the divisibility sieve IS a covering system; sieve-modulus unboundedness (S551) = Erdős minimum-modulus; import Hough's bound — MOST PROMISING. (2) dual rhythms: pair-differences=collision beats, pair-sums=pinch clocks. (3) CRT-fold 13 runners→1 super-runner. (4) orbit as cyclic code (H-spectrum). (5) reverse: box radiates, non-illuminated integer directions = Nullstellensatz zero set (= polynomial method, S559). (6) Farey/Stern-Brocot dissection of time.

**Efficient determination:** ~30 random clocks dispose of the easy 99.x%; the single n-gear check resolves the rare resonant residual; the only object a proof must control is the n-gear on the resonance-maximal family = a covering-system of prime gears.

**Honest:** split + pair-sum/n-clock + prime-gear model grounded; covering-system bridge and other remodels are speculation (but covering-systems is a mature, un-mined toolset converging on the same object).

**See:** `07-reflections/lrc-which-clocks-matter-prime-gears-and-wild-remodels-s564.md`, `04-computation/lrc_which_clocks_gears_s564.py` (+.out); HYP-2080 (orbit/resonance), HYP-2075 (multi-sieve), HYP-2063 (apex/CRT), THM-369.
