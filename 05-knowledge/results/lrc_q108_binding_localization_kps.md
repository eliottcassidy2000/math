# OPEN-Q-108 binding case localized to a FEASIBLE finite window (kps-S23)

The sole LRC(14) residual (span>14 => p0<=cap) has its `sup p0` localized exactly:

## sup p0 over wide k=9 configs, by span band (cap_9=0.4943):
```text
span (15,20]: 0.369   worst = [0,1,2,3,4,5,6,7,18]   (consec_8 + ONE far)
span (20,28]: 0.372   worst = [0,1,2,3,4,5,6,7,21]
span (28,40]: 0.368   worst = [0,1,2,3,4,5,6,7,29]
span (40,70]: 0.157   <- DROPS sharply (decorrelation onset ~ 5*max(base))
span (70,140]:0.108
```

## Findings
- **The binding wide case = "bounded base + ONE far element `w`", with `w` up to ~40** (the SINGLE-FAR case,
  THM-546/557 regime). p0 ~ 0.37 < cap_9 (margin 0.12). Multi-far (far_count>=2) is COMFORTABLE (p0 ~ 0.13).
- **The actual `Delta_w` is tiny** (Plat(consec_8)=0.362, p0(consec_8+18)=0.369 => Delta=0.007 << margin 0.13);
  the comb bound `(6/49)V/w` is just LOOSE at small w (it gives ~1.8 at w=18). So this is NOT an analytic gap —
  it is the looseness of the single-far bound at small w.
- **Decorrelation onset ~ `5*max(base)`**: for bounded base (max<=14) the far element decorrelates (p0 drops)
  once `w > ~70`; the comb bound `(6/49)V/w` then rigorously gives `|Delta_w|<margin` for `w` past the cutoff.

## The route (now with a FEASIBLE finite window)
1. `span<=14`: finite check. DONE.
2. bounded base (`E\{w}` span<=14) + single far `w in (14, W*)`: FINITE check. `W*` set by the single-far cutoff
   (`(6/49)V/W* = margin`); with the actual small `V` (and a SHARP `O(m)` bound from the transfer-operator /
   expanding-map spectral gap) `W* ~ 70-80`, giving ~`C(14,8)*65 ~ 2e5` sets — FEASIBLE.
3. bounded base + far `w > W*`: comb bound (THM-546/557), RIGOROUS.
4. wide base / multi-far: p0 comfortable (<0.16<<cap), via decorrelation (carrier product, splitting lowers cover).

## The ONE thing that closes it
A **sharp single-far bound** `|Delta_w| <= C_sharp(m)/w` with `C_sharp = O(m)` (not the lossy `O(m^2)`), so the
cutoff `W*` drops to ~tens and the finite window (step 2) is feasible. The transfer-operator/expanding-map framing
(the anchor map x->frac(wx) is the w-fold expanding map; its spectral gap is the `1/w` rate) is the route to it.
-> OPEN-Q-108, THM-546/547/557, HYP-2705 (cat-map), the S23 closure workflow.
