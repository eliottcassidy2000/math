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

## Structure (the L7 closure)
- DECORRELATED main term: `sum_t P_t^{(r)} p_t(B)` -- CLOSED FORM, < cap (finite check over bounded bases, all r).
- RESONANT correction: for a COMMENSURABLE far group, P_t^{(r)} is replaced by the curve survival width
  P_curve(S) (HYP-2757 atlas) -- finite ratios, also < cap.
- This is the moment dual (THM-534) with EXPLICIT far-coverage weights; the resonance is the only deviation.
Caveat (MISTAKE-082): the decorrelated value is the M->inf FLOOR; the actual p0 = it + a (resonance) error,
so the bound is `p0 <= [decorrelated] + [resonance correction]`, both controlled. -> OPEN-Q-108, HYP-2757, THM-534.
