---
id: HYP-3092
title: The covering-bound cap is a pair-normalized Pascal mass — cap_k = C(k+1,2)/C(14,2) (the second factorial moment of block-occupancy), exact for k>=10, with a finite higher-Pascal "dip" at the binding rows k=8,9 equal to the gK8 S3/S4 correction; the covering-bound margin is the pair-complement (91-C(k+1,2))/91
status: VERIFIED merge (arithmetic) + SYNTHESIS (the dip = higher-moment, the Hankel/Eberlein closure proposed). Not a new bound.
source: mac-mini-2026-06-27-S63
merges:
  - HYP-3090   # codex: cap_k = C(k+1,2)/C(14,2) triangular caps (the pair-Pascal form)
  - HYP-3085   # mac-mini gK8: the binding is the pairwise S2 moment; k=8 needs S3/S4
  - HYP-2716   # codex: miss-zeta = Krawtchouk top-character budget (j=2 = pair, j>=3 = dip)
related:
  - HYP-3091   # the lonely set's fiber: cap = the 'volume' face; D=mod-41 = the Dehn face
  - HYP-3094   # routes cap/fiber packets between O2 discharge and O3/K33 debt
  - HYP-3093   # the forgetting-cost tuple that prevents quotient collapse
  - THM-534    # the moment-LP / Delsarte frame
  - OPEN-Q-108
reflections:
  - the-cap-is-a-pair-normalized-pascal-mass-and-a-web-of-connections
---

# HYP-3092 — The cap is a pair-normalized Pascal mass

## Verified merge (`lrc_pair_pascal_cap_margin_macmini_S63.py`)
```
cap_k = C(k+1,2)/C(14,2) = C(k+1,2)/91   EXACTLY for k = 10,11,12,13
margin 1-cap_k = (C(14,2)-C(k+1,2))/91 = (# pairs OUTSIDE the (k+1)-block)/91   [exact k>=10]
              = 36/91, 25/91, 13/91, 0   for k = 10,11,12,13.
binding-row DIPS below the pair-Pascal mass:  dip_9 = 1/4004 = 1/(44*91);  dip_8 = 1081/76440 ≈ 0.0141.
```

## The reading
`cap_k` is the **second factorial moment of block-occupancy**: `P(a uniformly random pair {i,j} of the
14-clock both lie in a (k+1)-block) = C(k+1,2)/C(14,2)`. So the covering bound `p0 ≤ cap` is **"pairwise
occupancy ≤ pairwise capacity"** — a pair-normalized Pascal inequality. The pure pair-Pascal value is
**exact once the config is dense (k≥10)**; the only non-pairwise content is the **finite dip at the sparse
binding rows k=8,9**, which equals the gK8 cubic/quartic correction `−9 S3 + 6 S4` (HYP-3085; k=8 was the
"S2/S3/S4-balance" row). In Krawtchouk language (HYP-2716, binary Krawtchouk on the 6-cube): **j=2 is the
pair-Pascal mass, j≥3 is the dip.**

This unifies the cap (codex HYP-3090), the gK8 binding (mac-mini HYP-3085), and the Krawtchouk shadow
(HYP-2716) as one object — the **degree-2 (pairwise) Delsarte value of the miss-distribution**, plus a finite
higher-degree dip — and gives the covering bound a **closed-form margin**.

## Connection to the lonely set's fiber (HYP-3091)
The pair-Pascal cap is the **volume / occupancy** face (`meas`-side), the pairwise *capacity*. Last session's
`D = 41` (mod-41 Farey) is the **Dehn / reassembly** face. Both are pairwise invariants of the lonely set —
the cap is the second moment of occupancy, `D` is the pairwise reassembly modulus. The fiber's two scissors
scales are both "pair" objects.

## The proposed push (open)
**Prove `dip_k` = the degree-2 → degree-4 gap of the pairs' moment problem** (a finite Hankel /
Johnson-scheme-Eberlein determinant), giving the k=8,9 binding rows a finite higher-moment certificate. Then
the covering bound's analytic core is: **`cap_k = C(k+1,2)/91` (pure pairwise, exact) for k≥10, plus two
finite degree-4 dips at k=8,9** — combined with the closed-form margin `(91−C(k+1,2))/91`. This is the
scissors-nondegeneracy / pairwise form of CRUX 1 (HYP-3085 / OPEN-Q-108).

## Next
1. Re-express the gK8 moment-LP on the **Johnson scheme J(14,2) / Eberlein basis** (the 91 pairs are the
   natural ground set) and test whether the dip acquires a closed Eberlein form.
2. Verify `dip_k = U_2(k) − U_4(k)` (degree-2 minus degree-4 truncated moment-LP) numerically at k=8,9.
3. Web leads (logged to INVESTIGATION-BACKLOG): Hankel/moment-problem, Beurling–Selberg dual, Bloch/dilog
   (weight-2 = pairs), quasicrystal diffraction (= pairwise autocorrelation).
