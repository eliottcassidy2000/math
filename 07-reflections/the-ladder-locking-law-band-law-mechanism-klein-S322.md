# The ladder-locking law: why the first gap populates exactly at N ≡ 1 (mod 6)

> **⚠ SCOPE BANNER (klein-S322, MISTAKE-194):** every EXACT-M fact in this note
> (ladder rung values, plateau law instances, trichotomy values, K-rungs) is
> referee-exact and UNAFFECTED. The census-completeness claims ("complete sub-1/13
> table", "all on-ladder", "zero off-ladder", h_min exactness) rest on the census
> harness whose mask prune was unsound; they are UNDER-REVERIFICATION by the
> patched suite and this banner will be updated with the outcome.


**Instance:** klein-2026-07-19-S322 (owner: "keep clarifying structure like this; see
how insights weave into the remaining LRC(14) frontier").
**Evidence:** exact-M sweep N = 6..13 speeds, both ladders, all rungs m,s ≤ 6
(`lrc_crossN_two_ladders_klein_S322.py` + frozen out; referee engine gate-validated).
**Status of claims:** the sweep facts are EXACT; the mechanism identification rests on
opus-S409's proved trichotomy; the "law" at general N is VERIFIED 6..13, not proved.

## 1. The sweep fact (exact, complete for N = 6..13, m = 1..6)

Define, at N speeds (LRC(N+1)), the two ladders rooted at the AP:

- **L-ladder** L_N(m) = {1..N−2, N} ∪ {(N−1)m} — L_N(1) = {1..N} = the AP.
- **K-ladder** K_N(s) = {1..N−1} ∪ {Ns} — K_N(1) = the AP; Kravitz rungs.

Then:

| N | L-profile (m = 1..6) | locked? |
|---|---|---|
| 6 | P ? ? ? ? ? | no |
| **7** | **P P F F F F** | **YES** |
| 8–12 | P ? ? ? ? ? | no |
| **13** | **P P F F F F** | **YES** |

where P = floor plateau (M = 1/(N+1)), F = the exact formula **M = m/((N−1)m + 5)**
(active pair (5, (N−1)m), D = m), and ? = neither. The K-ladder is **unconditional**:
M(K_N(s)) = s/(Ns+1) exactly, all 16 tested (N, s) — the universal spine.

**The L-ladder locks (full P P F F F… profile) at exactly N = 7 and N = 13 — the two
N whose first gap is populated in the single-far stratum (THM-1284) — and rung 3 of a
locked ladder IS the mediant 3/(3N+2), inside the window.**

## 2. The m = 2 trichotomy = opus-S409's proved GW law (independent convergence, same day)

The derailed L_N(2) values are not noise; they follow exactly:

> M(L_N(2)) = 1/N (N even) · 2/(2N+1) (N odd, N ≢ 1 mod 6) · **1/(N+1) (N ≡ 1 mod 6 — TIGHT)**

and [odd ∧ 3 | 2N+1] ⟺ N ≡ 1 (mod 6). This is precisely the trichotomy opus-S409
proved today (16/16 exact + the gate mechanism: **3a ≡ ±2 phase-denial blocks the D=2
escape at 3 | 2N+1**): their "D=2 escape" is the 2/(2N+1) branch, their gate firing is
the collapse to the floor. **L_N(2) is the GW family**: the tight locus {AP, GW} is
rungs m = 1, 2 of the L-ladder, and GW is tight exactly at lock-N.

## 3. What this explains (the weave into the frontier)

1. **THM-1284's band law gets a mechanism.** "First gap populated ⟺ N ∈ {6} ∪
   {N ≡ 1 mod 6}, arithmetic not width" was a validated mystery. Now: mediant
   attained (single-far stratum) ⟺ the L-ladder locks ⟺ GW is tight ⟺ opus's
   gate fires ⟺ N ≡ 1 (mod 6). The N = 6 member (5/33, a k=3 stratum value) is
   off-both-ladders — matching THM-1284's own classification that only N = 7, 13
   populate the single-far stratum. A proof of the band law in this stratum =
   opus-S409's gate (proved) + the locked-formula continuation for m ≥ 3
   (verified; the (5, (N−1)m) witness surviving the base — one lemma to write).
2. **The S-T k=12 vs k=13 asymmetry is retrodicted.** k=12 (N=12 ≢ 1 mod 6):
   tight locus {AP} alone (S-T's uniqueness); k=13 (N=13 ≡ 1 mod 6): {AP, GW}
   (kps HYP-7940's terminal table T=2). The "extra" tight family at n=14 exists
   because 6 | 12 = N−1. **The next lock is N = 19 (LRC(20))** — a falsifiable
   cross-N prediction: GW₁₉ = {1..17,19,36} should be tight and 3/59 attained.
3. **Wall A's rigidity core in ladder language.** "Only {AP, GW} are tight at
   n=14" ⟺ "no third family plateaus at the floor" — plateau-uniqueness. The
   bottom-census (exhaustive to h = 45: eight families, all on-ladder) is the
   height-bounded version of exactly this statement.
4. **The 5.** The universal active leg 5 = numerator 3 + slack 2 of the mediant
   (s_pair = 3N + 2 = 3(N−1) + 5) — the same 5 in m/(12m+5), in THM-1284's
   binding pair {5, 3(N−1)}, and in the L-ladder's formula at every N. The
   frontier's (D, s)/slack/k coordinates all meet in this one integer: slack of
   the mediant = 2 = its Fan–Sun k, and the leg 5 = D + k = 3 + 2. Conjecturally
   the locked ladder's rung-m active leg is always D + k = m + 2 = 5 at m = 3.

## 4. Honest scope and named next steps

- The sweep is N ≤ 13, m ≤ 6; the law is a verified pattern + a proved m=2
  mechanism (opus-S409), NOT yet a theorem for m ≥ 3 or general N. The missing
  lemma is clean: *at lock-N, the (5, (N−1)m) witness t\* = (a₁+a₂)/((N−1)m+5)
  clears the base {1..N−2, N} for all m ≥ 3.* One session, likely Lean-able.
- **N = 19 verification: CONFIRMED (same session, referee-exact).** M({1..17,19,36})
  = 1/20 (GW₁₉ TIGHT), M({1..17,19,54}) = 3/59 = the mediant, inside LRC(20)'s
  first window, and L₁₉(4) = 4/77 on-formula — the lock law's first out-of-sample
  prediction lands; THM-1284's populated band extends {6,7,13} → {6,7,13,19}.
- Cross-link kps HYP-8040/8050 (escape classes): the ladder families are the
  ON-shell bottom; the escape atlas (tiny-anchor + mid-band clusters,
  q − v_max ∈ [4, 8]) is the OFF-shell boundary — the two structures partition
  the near-floor landscape the small-witness law must govern.

## 5. Cross-links

opus-S409 (the proved trichotomy + gate — the mechanism) · THM-1284 (the band law
this explains; HYP-4516 gate solutions = the locked witnesses) · THM-1230/1290/1291
(ladders, census, CF leg) · HYP-7940 (terminal table {AP,GW} = the m ≤ 2 plateau) ·
S-T k=12 uniqueness (retrodicted) · the S321 reflection (the N = 13 instance this
generalizes) · kps HYP-8040/8050 (the off-shell complement).
