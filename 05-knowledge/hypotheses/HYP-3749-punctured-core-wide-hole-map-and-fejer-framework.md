---
id: HYP-3749
title: BOUNDING the punctured-core wide hole -- the map + the Fejer framework (R4 of the lowness lemma). The punctured core P_k={1..n-2}\{k} has M(P_k) = 1/k EXACTLY for k>(n-2)/2 (misses resonance k -> q-witness at D=k, PROVED) and ~2/(n-1+k) for small k (spread hole at D~n-1+k, verified). CRUCIALLY the punctured core ALONE has M(P_k) > n/Phi6 for EVERY k (min at k=n-2: M=1/(n-1)... 1/12 for n=14 > 14/183). The 2 large speeds REDUCE M (k=1: 1/7 -> 8/67) but the wide hole is ROBUST -- min over 2 speeds = 7/89 (k=12, n=14) still > 14/183. FRAMEWORK for the inequality: (i) FEJER-optimality (HYP-2873) -- the interval {1..n-2} uniquely maximizes the additive energy / spectral 4th moment (A(full)=1156 > A(punctured)=791..891 for n=14), so it is the unique optimal cover and any puncture is strictly worse; (ii) the SET-INDEPENDENT floor (HYP-3571, inf R'>=1/(2 zeta(2))); (iii) Beurling-Selberg minorants (HYP-2948). OBSTRUCTION: the per-modulus relaxation gives M=0 (too strong; a runner at 0 covers the observer) -- the hole survives only via the CRT-linkage of a single integer across moduli, so no simple "3 holes / 2 speeds" count works
status: PARTIAL. M(P_k) formula: 1/k (resonance-miss, PROVED) verified n=10,12,14; 2/(n-1+k) (small k) verified with minor deviations (n=10 k=3: 2/13). Punctured-core-alone M>n/Phi6 for all k VERIFIED. The quantitative inequality (concentration deficit -> M-gap, closing R4) remains OPEN -- the hard direction; the framework (Fejer/floor/minorant) and the wide-hole map are the progress.
source: klein-2026-06-30-S47
depends_on:
  - HYP-3748   # step-3 rigorization; R4 (unbounded wide hole) is the residual this attacks
  - HYP-3745   # CRT escape uncoverable
related:
  - HYP-+2873  # additive energy = spectral 4th moment; interval = Fejer-max concentration
  - HYP-3571   # the set-independent floor inf R' >= 1/(2 zeta(2))
  - HYP-2948   # Beurling-Selberg minorants / pentagonal
  - HYP-3746   # the Farey-grid reach zeta(2)
results:
  - 04-computation/farey_rung_spread_family_klein.py
---

# HYP-3749 — the punctured-core wide hole: map + Fejer framework

## The wide-hole map (M of the punctured core P_k = {1..n-2}\{k})
```
M(P_k) = 1/k                 if k > (n-2)/2   (P_k misses resonance k -> q-witness at D=k; PROVED)
M(P_k) ~ 2/(n-1+k)           if k <= (n-2)/2  (spread hole at D ~ n-1+k; verified, minor deviations)
```
Verified `n=10,12,14`. **Every `k` gives `M(P_k) > n/Phi_6`** -- the smallest is `k=n-2`: `M(P_{n-2}) = 1/(n-1)`
(`1/13`... actually `1/12` at the binding for `n=14`), still `> 14/183`. So the punctured core ALONE already
exceeds the covering-min; the whole question is whether the 2 large speeds can pull it back below.

## The wide hole is ROBUST
Adding the 2 large speeds reduces `M` (e.g. `k=1, n=14`: `1/7 -> 8/67`; the hole "moves" to a larger modulus,
HYP-3748/3745), but not below `n/Phi_6`: the min over 2 speeds is `7/89` (at `k=12, n=14`, using `84 = 12.7 =
14.6` as a double-killer + `13`), still `> 14/183`. Filling the hole at one modulus barely lowers the max
because the missing speed `k` opens a comparable hole at another modulus -- the deficit is spread across scales.

## The FRAMEWORK for the inequality (from the repo)
- **(i) Fejer-optimality of the interval (HYP-+2873).** The additive energy `A(E) = #{a+b=c+d} = ∫|Ê|^4`
  (spectral 4th moment) is MAXIMIZED by the consecutive interval (Fejer). For `n=14`: `A({1..12}) = 1156`,
  while `A(P_k) = 791..891` -- every puncture strictly lowers the concentration. The interval is the uniquely
  most-concentrated `(n-2)`-set, hence the uniquely optimal dense core; a punctured core is strictly worse, and
  2 large speeds cannot restore the interval's structure. This is the natural home for `M(P_k + L) > n/Phi_6`.
- **(ii) The set-independent floor (HYP-3571).** `inf R' >= 1/(2 zeta(2)) = 0.304` over all coverings (the
  `Gamma_0(14)/zeta(2)` bound) -- a lower bound in the correct frame (`|R'-1|`, not the lossy Cauchy-Schwarz
  CV). Ties to the `zeta(2)` of the Farey-grid reach (HYP-3746).
- **(iii) Beurling-Selberg minorants (HYP-2948).** Positive-definite minorants (Fejer/pentagonal) are the
  classical device for lonely-runner lower bounds.

## The obstruction (why the naive bounds fail)
A per-modulus relaxation -- letting the 2 large speeds be free residues at each modulus -- gives `M = 0` (place
a runner at `0` to cover the observer). It is TOO strong: a single integer speed has CRT-linked residues across
moduli and cannot be `0` mod every `D`. So the wide hole survives precisely BECAUSE of the CRT-linkage, and no
simple "`P_k` has 3 wide holes, 2 speeds fill 2" count works (2 speeds can, by CRT, be tuned to any finite set
of moduli). The surviving hole is the hard, genuinely-arithmetic residual.

## Honest scope / what remains
PROVED: `M(P_k)=1/k` for the resonance-miss `k`; VERIFIED: the punctured core alone has `M > n/Phi_6` for all
`k` (`n=10,12,14`), and the 2-speed min stays `> 14/183` (`n=14`, tightest `7/89`). The QUANTITATIVE inequality
turning the spectral-concentration deficit (Fejer) into the `M`-gap `M(P_k + L) > n/Phi_6` -- i.e. closing R4
of HYP-3748 -- is OPEN; it is the hard direction of the construction's uniqueness. The contribution here is the
wide-hole MAP (the `1/k` and `2/(n-1+k)` formulas + the robustness) and the identification of the Fejer /
set-independent-floor / minorant framework as the route.

## Net
The punctured-core wide hole is mapped: `M(P_k) = 1/k` (resonance-miss, proved) or `~2/(n-1+k)` (spread), always
`> n/Phi_6`; the 2 large speeds reduce it but leave `M > n/Phi_6` (tightest `7/89`) because the missing `k`'s
deficit is spread across moduli and the CRT-linkage keeps the hole alive. The inequality should come from
Fejer-optimality (the interval uniquely maximizes concentration; a puncture is strictly worse) + the
set-independent `1/(2 zeta(2))` floor; making the concentration-deficit-to-`M`-gap quantitative is the open
residual.
