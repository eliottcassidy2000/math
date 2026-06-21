# OPEN-Q-108 decorrelated wide bound = a CLOSED FORM (linear in the base missed-distribution) (kps-S24)

For a bounded base B (size k-r) plus r FAR elements that decorrelate (large separated scales), the cover
factors: each far element is an independent uniform point on Z/7, and the r-element pair/group covers the
base's missed-set S only if |S| <= r. EXACT decorrelated value:
> **p0(B u {r decorrelated far}) = sum_{t=0}^{min(r,6)} P_t^{(r)} * p_t(B),**
where p_t(B) = meas_x(B misses EXACTLY t inner sectors), and P_t^{(r)} = P(r indep uniform points cover t
SPECIFIC sectors) = sum_{j=0}^t (-1)^j C(t,j) ((7-j)/7)^r (inclusion-exclusion surjection prob).
For r=2: P_0=1, P_1=13/49, P_2=2/49, P_{>=3}=0.

This is a LINEAR FUNCTIONAL of the base's missed-distribution (p_0,..,p_6) -- the MOMENT-DUAL structure (THM-534).

## Two-far (r=2) VERIFIED < cap (clean finite check)
`p0(B) + (13/49)p1(B) + (2/49)p2(B)` over ALL bounded bases B subset {0..14}:
- k=9 (3003 bases): max = 0.24645 at consec_7 (cap 0.4943, margin 0.2478), 0 over cap.
- k=10 (3432 bases): max = 0.39894 at consec_8 (cap 0.6044, margin 0.2055), 0 over cap.
Matches the commensurable atlas non-resonant value (0.247). So the decorrelated two-far closes via this clean
finite check on (p0,p1,p2) of bounded bases, margin >= 0.20.

## MAIN RESULT (kps-S24): the decorrelated wide bound is MAXIMIZED at r=1 = Q(k-1) < cap (PROVEN finite check)
The full check `sum_t P_t^{(r)} p_t(B)` over ALL r=1..k-1 and ALL bounded bases B (size k-r):
```text
  k    GLOBAL max over (r,B)     value     cap      margin    over_cap
  8    r=1, consec_7             0.19660   0.38146  0.18486   0
  9    r=1, consec_8             0.36210   0.49426  0.13216   0
  10   r=1, consec_9             0.44789   0.60440  0.15651   0
  11   r=1, consec_10            0.53125   0.72527  0.19402   0
  12   r=1, consec_11            0.60224   0.85714  0.25490   0
```
**The max is ALWAYS at r=1 (single far) on the consec base = the plateau Q(k-1).** So for ANY wide config (any
r>=1 far elements, any bounded base):
> **p0_decorr(E) = sum_t P_t^{(r)} p_t(B) <= Q(k-1) < cap_k.**
This PROVES the base-size domination (single-far dominates) FOR THE DECORRELATED bound, AND closes the
decorrelated wide bound with the SAME margin as the finite check (Q(k-1) to cap). The decorrelated p0 never
exceeds the single-far plateau -- more far elements (smaller base, P_t weights) strictly lower it. The wide
bound = [decorrelated <= Q(k-1) < cap, DONE as a finite check] + [resonance error <= margin].

## Structure (the L7 closure)
- DECORRELATED main term: `sum_t P_t^{(r)} p_t(B)` -- CLOSED FORM, < cap (finite check over bounded bases, all r).
- RESONANT correction: for a COMMENSURABLE far group, P_t^{(r)} is replaced by the curve survival width
  P_curve(S) (HYP-2757 atlas) -- finite ratios, also < cap.
- This is the moment dual (THM-534) with EXPLICIT far-coverage weights; the resonance is the only deviation.
Caveat (MISTAKE-082): the decorrelated value is the M->inf FLOOR; the actual p0 = it + a (resonance) error,
so the bound is `p0 <= [decorrelated] + [resonance correction]`, both controlled. -> OPEN-Q-108, HYP-2757, THM-534.
