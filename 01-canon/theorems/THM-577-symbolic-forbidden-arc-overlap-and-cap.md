---
id: THM-577
title: Closed-form forbidden-arc overlap (apex-14 threshold) ⟹ symbolic evaluation of the coverage-extremal caps cap_11=66/91 and cap_10=55/91 (closing the value half of THM-576's j=3 search-only case); the binding dip at k=8,9 is a finite computable higher-order remainder
status: PROVED (overlap lemma: rigorous derivation + verified 0 mismatches over all 120 pairs p,q<=16; symbolic cap_11, cap_10 by assembly). The OPTIMALITY of the minimizer (that {1,12,13} attains the min) remains from THM-576's exact search; this theorem makes its VALUE symbolic.
source: mac-mini-2026-06-27-S64
depends_on:
  - THM-576   # cap_k = min meas(lonely(P)) = pairwise avoidance; identifies the minimizers (by search)
  - THM-534   # the moment-LP / inclusion-exclusion frame
related:
  - HYP-3092  # cap_k = C(k+1,2)/91 = pair-normalized Pascal mass (the values)
  - HYP-3090  # codex: triangular caps
  - OPEN-Q-108
external: Lonely Runner Conjecture (first open case = 13 speeds).
---

# THM-577 — Closed-form forbidden-arc overlap and the symbolic caps

## Setup
For a speed `p`, the **forbidden set** is `Fb(p) = {x∈[0,1): ‖p x‖ < 1/14}` — `p` arcs of half-width
`1/(14p)` centered at `a/p`, total measure `1/7`. The coverage-extremal cap (THM-576) is
`cap_k = min_{|P|=j} meas(lonely(P))`, `j = 13−k`, `lonely(P) = ⋂_{p∈P} Fb(p)^c`, evaluated by
inclusion–exclusion over the joint forbidden measures `F(T) = meas(⋂_{p∈T} Fb(p))`.

## Lemma (closed-form pairwise overlap; the apex-14 threshold)
For positive integers `p < q` with `gcd(p,q)=1`,
```
o(p,q) := meas(Fb(p) ∩ Fb(q)) = 1/(7q) + (1/(7pq)) · Σ_{m≥1, 14m<p+q} (p + q − 14m),
```
where the offset sum is **empty unless `p+q > 14`** (the apex `14 = 2·7`). For general `gcd(p,q)=g`,
each residue is realized `g` times and `m` runs over multiples of `g`. (Caveat: when `p=1, q≥15` the
`m=1` term must be capped at the concentric value `1/(7q)`; for all `p,q ≤ 13` — the minimizer regime —
no cap is needed.)

**Proof.** By CRT the map `(a,b) ∈ ℤ/p × ℤ/q ↦ (aq−bp) mod pq` is a bijection (coprime case). The arc of
`Fb(p)` at `a/p` and the arc of `Fb(q)` at `b/q` are centered at distance `‖(aq−bp)/(pq)‖`; writing the
offset as the integer `m = aq−bp` (centered representative), two arcs of half-widths `1/(14p),1/(14q)`
overlap in length `2/(14q)` if `m=0` (concentric) and `(p+q−14|m|)/(14pq)` if `0<|m|` and that is positive.
Summing the bijective family over `m` (with `±m` giving the factor in the sum) yields the stated formula. ∎
**Verification:** `04-computation/lrc_symbolic_overlap_closedform_macmini_S64.py` — 0 mismatches against the
exact rational measure over all 120 pairs `p,q ≤ 16` (coprime and non-coprime).

Examples (the minimizer top-pairs):
`o(1,13)=1/91` (p+q=14, offset OFF), `o(1,12)=1/84` (off), `o(12,13)=1/91+11/1092=23/1092` (p+q=25, ON).

## Theorem (symbolic evaluation of cap_11, cap_10)
Assembling the closed forms over the THM-576 minimizers (all pairwise coprime):
```
cap_11 = meas(lonely{1,13})    = 1 − 2/7 + o(1,13)                          = 5/7 + 1/91     = 66/91 = C(12,2)/91.
cap_10 = meas(lonely{1,12,13}) = 1 − 3/7 + [o(1,12)+o(1,13)+o(12,13)] − F(1,12,13)
       = 4/7 + (1/84 + 1/91 + 23/1092) − 1/91 = 4/7 + 3/91                  = 55/91 = C(11,2)/91.
```
The triple `F(1,12,13) = 1/91 = 1/(7·13)` is the **narrowest nested arc** inside the speed-1 container
`I_0=(−1/14,1/14)` (the speed-13 central arc, contained in the speed-12 central arc). Both caps thus have a
**fully symbolic** value via the overlap lemma — closing the VALUE half of THM-576's `j=3` (k=10) case, which
was previously only exact-by-search. (`cap_12=6/7`, `cap_13=1` are the trivial `j=1,0` rows.)

## The binding dip (k=8,9) as a finite higher-order remainder
Write the inclusion–exclusion as `meas(lonely(P)) = (1 − j/7) + Σ_{r≥2} O_r`, `O_r=(−1)^r Σ_{|T|=r}F(T)`.
Using the algebraic identity `C(14−j,2)/91 = 1 − j/7 + C(j,2)/91`, the cap equals the pair-Pascal mass minus
a **dip**:
```
cap_k = C(k+1,2)/91 − dip_j,   dip_j = C(j,2)/91 − Σ_{r≥2} O_r.
dip_j = 0 (j≤3, k≥10);   dip_4 = 1/4004 (k=9);   dip_5 = 1081/76440 (k=8).
```
The orders (`lrc_symbolic_coverage_inclexcl_macmini_S64.py`): at `j=3`, `O_2=4/91, O_3=−1/91` (net `3/91`,
dip 0); at `j=4` the net of `O_2,O_3,O_4` falls `1/4004` short of `C(4,2)/91`; at `j=5` (the minimizer
**breaks** to the middle-spread, 3-correlated `{1,5,7,8,9}`) it falls `1081/76440` short. **Each dip is a
finite sum of closed-form `o`/`F` terms** (the lemma extends to triples as nested/offset arcs), so the
binding rows reduce to a finite symbolic remainder, not a search.

## Net / what remains
- **Symbolic VALUE of the coverage extremal: PROVED for k≥10** (cap_10..13 = `C(k+1,2)/91`, closed form),
  and the k=8,9 dips are finite computable remainders of the same closed-form overlaps.
- **OPTIMALITY** (that the named config attains the min over all `P`) is still THM-576's exact search; the
  symbolic optimality (a packing/spreading argument that the top-cluster minimizes the inclusion–exclusion)
  is the remaining gap = the concentration extremality (OPEN-Q-108 / gK8 HYP-3085). The apex-14 threshold
  in `o(p,q)` (offset ON iff `p+q>14`) is the structural lever for that argument: overlaps grow only across
  the apex, so spreading speeds to the cluster `{1, …, 13}` extremes minimizes total overlap.
