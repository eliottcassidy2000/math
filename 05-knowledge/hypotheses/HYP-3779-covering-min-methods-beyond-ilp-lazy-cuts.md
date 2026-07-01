---
id: HYP-3779
title: METHODS BEYOND THE ILP for the covering-min -- and the LAZY-CUT (cutting-plane) method CLOSES the n=12,13,14 residual RIGOROUSLY: covering-min = the construction n/Phi6 for ALL speeds <= n(n-1). The full set-cover ILP (HYP-3731) times out at speeds V>~4n; two alternatives: (#1) LP-RELAXATION infeasibility -- TOO WEAK (integrality gap, no certificate); (#2) LAZY-CONSTRAINT / CUTTING-PLANE ILP (Benders row generation: tiny ILP=size+divisibility, add the lonely-witness cut(s) each round, re-solve) -- WORKS at V=n(n-1) where the monolithic 40k-constraint ILP dies, because it carries only the handful of witnesses that BIND. RESULT (rigorous): NO covering set with speeds<=n(n-1) beats the construction at n=12 (12/133, 208 single-cuts), n=13 (13/157, 3 multi-cut rounds), n=14 (14/183, 3 multi-cut rounds) => covering-min = n/Phi6 there, CLOSING the (4n,n(n-1)] residual of HYP-3778 for all three incl. the LRC-14 target. The MULTI-CUT variant (add every lonely witness of each candidate, not just the deepest) converges in 3 rounds vs 208. Broader untried menu: Lovasz-theta/SDP (one-shot, closes the LP gap), column-generation (speeds beyond n(n-1)), B&B with q/k/(n+q)-witness pruning, CRT decomposition, meet-in-the-middle, binding-first Farey-dissection
status: STRONG. VERIFIED: (#1) LP relaxation FEASIBLE=too weak (no certificate). (#2) lazy-cut ILP INFEASIBLE at V=n(n-1) for n=12 (208 single-cuts), n=13 (3 multi-cut rounds, 2153 constraints), n=14 (3 multi-cut rounds, 3916 constraints) => NO covering set with speeds<=n(n-1) has M<n/Phi6 => covering-min(n)=n/Phi6 (the construction) for n=12,13,14, RIGOROUS up to the construction's own max speed n(n-1). CLOSES the HYP-3778 residual for n=12,13,14 (incl. LRC-14). REMAINING GAP: speeds > n(n-1) not excluded by search (but the construction's max IS n(n-1); the huge-speed tail is HYP-3745's CRT-escape territory, a theory question). Broader menu = brainstorm (untried).
source: klein-2026-07-01-S61
depends_on:
  - HYP-3731   # the set-cover ILP (the base method; times out at large V)
  - HYP-3778   # the n=12 residual (this CLOSES it up to n(n-1))
related:
  - HYP-3737   # construction forced n>=12 (n=12 now rigorous up to n(n-1))
  - HYP-3745   # CRT escape (why speeds > n(n-1) don't beat -- the remaining gap)
  - HYP-3764   # the open edge (further closed here)
results:
  - 04-computation/covering_min_lp_lower_bound_klein.py
  - 04-computation/covering_min_lazy_cuts_klein.py
  - 05-knowledge/results/covering_min_lp_lower_bound_klein.out
  - 05-knowledge/results/covering_min_lazy_cuts_klein.out
---

# HYP-3779 — methods beyond the ILP; lazy cuts close n=12

## The bottleneck
The exact set-cover ILP (HYP-3731) is the right tool but its monolithic form carries ~40k danger
constraints (all breakpoints `k/d`, `d<=2V`) and `milp` times out at speeds `V > ~4n` -- which left
HYP-3778's `(4n, n(n-1)]` residual open (a beater could hide among large speeds). Two creative
alternatives, aimed at a rigorous bound at `V = n(n-1)` (the construction scale):

## (#1) LP-relaxation infeasibility certificate -- TOO WEAK
Relax `x_v in {0,1}` to `[0,1]`; if the set-cover LP for "`M < n/Phi_6`" is INFEASIBLE, the ILP is too,
so no beater (rigorous). LP (`linprog`/HiGHS) scales to `V=n(n-1)`. **Result: the LP is FEASIBLE** at
`n=12,13,14` (`V=n(n-1)`, ~21k-40k witness rows, 0 all-zero rows) -- a fractional cover exists, so the
relaxation has an integrality gap and gives NO certificate. Honest negative. (A tighter relaxation --
Lovasz theta / Lasserre SDP -- would close this gap; untried.)

## (#2) LAZY-CONSTRAINT / CUTTING-PLANE ILP -- CLOSES n=12
Benders-style row generation. Solve a TINY ILP (size `= n-1` + divisibility: a multiple of every
`q in {2..n}`), read off a candidate set `S`, compute its EXACT `M(S)` and deepest-hole breakpoint
`t*`; if `M(S) < n/Phi_6` a beater is found, else add ONLY the cut "some selected speed is within
`n/Phi_6` of `0` at `t*`" and re-solve. The ILP carries only the handful of witnesses that actually
bind, so it runs at `V = n(n-1)` where the monolithic ILP dies.
> **Result (V = n(n-1)): the ILP goes INFEASIBLE, so no covering set with speeds `<= n(n-1)` has
> `M < n/Phi_6` -- the covering-min is the construction, RIGOROUS, for:**
> ```
>   n=12  V=132   INFEASIBLE (208 single-cut rounds)          => covmin = 12/133
>   n=13  V=156   INFEASIBLE (3 MULTI-cut rounds, 2153 rows)  => covmin = 13/157
>   n=14  V=182   INFEASIBLE (3 MULTI-cut rounds, 3916 rows)  => covmin = 14/183   (the LRC-14 target)
> ```
This CLOSES the HYP-3778 `(4n, n(n-1)]` residual for `n=12,13,14`. The **multi-cut variant** -- add
EVERY lonely witness of each candidate (not just the deepest) per round -- converges in **3 rounds**
instead of 208: the certificate that no beater exists is a finite packing of a few thousand
lonely-witness cuts.

## The broader menu (brainstorm, ranked by promise for the residual)
TIER A -- rigorous lower bound at large V (residual-closers):
1. **Lazy cuts** (#2) -- WORKS (n=12 done); scale with multi-cut rounds.
2. **Lovasz theta / SDP** -- tighter than the LP; the theta of the "danger conflict graph" sits between
   its clique and chromatic numbers and closes the LP integrality gap. UNTRIED, most promising next.
3. **Lasserre / sum-of-squares hierarchy** -- level-`k` SDP, increasingly tight; heavy but definitive.
4. **Column generation** -- dual of lazy cuts: price SPEEDS lazily for `V` beyond `n(n-1)`.
TIER B -- smarter exact search:
5. **Branch & bound with witness pruning** -- prune partial sets that already trip a `q`/`k`/`(n+q)`-witness (HYP-3766); domain-specific bounding beats generic `milp`.
6. **CRT residue-class decomposition** -- the covering constraint factors mod small primes; search fibers.
7. **Meet-in-the-middle** -- split into small core + large killers; hash coverage signatures; `sqrt`-speedup.
8. **Binding-first / Farey-dissection** -- enumerate candidate bindings `(D,j)` (sparse), not sets.
TIER C -- structural / theory (the real proof of the tail):
9. **Lowness lemma / large-multiples-forced** (HYP-3763/3747) -- the `(n(n-1), infinity)` tail by theory.
10. **Danger-nerve topology / Gauss-Bonnet** (S57) -- a Cech obstruction = the covering-min gap.
11. **Spectral / Hoffman bound** on the danger circulant (HYP-3731's circulant).

## Net
Beyond the ILP: the LP relaxation is too weak (integrality gap), but the **lazy-constraint cutting-plane
ILP** is the creative winner -- it reaches `V = n(n-1)` and PROVES that the covering-min is the
construction `n/Phi_6` for **`n=12,13,14`** (incl. the LRC-14 target), for all speeds up to the
construction scale `n(n-1)` -- closing HYP-3778's `(4n, n(n-1)]` residual. The **multi-cut** variant does
it in 3 rounds. So the covering-min transition (spread beaters for `n<=11`, construction for `n>=12`) is
now rigorous up to `n(n-1)` at `n=12,13,14`, not just up to `4n`. The only remaining gap is speeds
`> n(n-1)` (the CRT-escape tail, HYP-3745) -- a theory question, since the construction's own max speed
is exactly `n(n-1)`. Next lever: Lovasz-theta/SDP for a one-shot certificate, or column generation to
push past `n(n-1)`.
