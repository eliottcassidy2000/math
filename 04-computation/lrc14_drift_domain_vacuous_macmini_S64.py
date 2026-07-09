"""
mac-mini-2026-07-09-S64 -- the 1/7-vs-2/7 tension RESOLVED, and the local drift-absorbed hembed
discharge (klein-S205 `minReach_ge_of_driftGap`) is shown VACUOUS for the proof.

THE TRUE THRESHOLD IS RATIO-DEPENDENT, NOT 1/7 OR 2/7.
klein-S205's hypothesis (sorry-free, correct) is
    hgap :  1/7 + 2 * spread * (a + g/2) / Vmax  <  g
with a tooth-free gap (a, a+g) subset [0,1] (so 0 <= a, a+g <= 1, hence g <= 1), the binding
v_i = Vmax - e_i, and |e_i| <= spread.  The tightest admissible bound is spread = Vmax - Vmin
(over the same index set the binding quantifies).  Since a >= 0 we have phi := a + g/2 >= g/2, so

    g  >  1/7 + 2*spread*(g/2)/Vmax  =  1/7 + spread*g/Vmax
    =>   g * (1 - spread/Vmax)  >  1/7
    =>   g * (Vmin/Vmax)        >  1/7                      [since 1 - spread/Vmax = Vmin/Vmax]
    =>   g  >  (1/7) * (Vmax/Vmin)  =  r/7,    r := Vmax/Vmin.

A gap on the circle has g <= 1, therefore the hypothesis is SATISFIABLE only when  r < 7.

BUT kps-S28's `spread13_lonely` (LRCSpread13.lean, kernel-pure, one line) already gives
    a <= |v_i| <= b  and  b <= 13a   ==>   Lonely 14 v (1/(a+b)),
i.e. loneliness for EVERY ratio r <= 13.  And {r < 7} is a STRICT SUBSET of {r <= 13}.

=> `minReach_ge_of_driftGap` can never fire on a case `spread13_lonely` does not already close.
   The drift-absorbed hembed discharge is VACUOUS for the proof.

AND IN THE ONLY REGIME THAT NEEDS PROOF (r > 13) THE LOCAL ROUTE IS PROVABLY IMPOSSIBLE:
there it demands g > r/7 >= 13/7 = 1.857 > 1, which no gap can supply.

WHY (the structural reason): the slow-fast frame freezes the teeth t_i = frac(e_i * tau) while the
fast phase phi = frac(Vmax*tau) sweeps one ruler period.  Over that period the teeth move by
e_i/Vmax <= spread/Vmax = 1 - 1/r.  For r > 13 that is > 12/13: the teeth traverse nearly a full
turn while phi sweeps once.  There is NO slow-fast separation, so no local (single j, single phi)
embedding can exist.  The 1/7 and 2/7 thresholds are red herrings; the obstruction is the RATIO.
"""
from fractions import Fraction as F

print("DOMAIN of klein-S205 minReach_ge_of_driftGap  (necessary condition g > r/7, and g <= 1)")
print(f"{'ratio r':>10}{'min gap g needed = r/7':>26}{'g <= 1 feasible?':>19}{'spread13_lonely covers?':>25}")
for r in [F(1), F(2), F(4), F(6), F(69, 10), F(7), F(13), F(14), F(20), F(50)]:
    need = r / 7
    print(f"{str(r):>10}{str(need):>26}{str(need < 1):>19}{str(r <= 13):>25}")

print()
print("CONCLUSION")
print("  drift lemma fires  <=>  r < 7.")
print("  spread13_lonely closes  all r <= 13   (kps-S28, kernel-pure, one line).")
print("  {r < 7} is a STRICT SUBSET of {r <= 13}  =>  the drift lemma is VACUOUS for the proof:")
print("  it never fires on a case spread13_lonely leaves open.")
print()
print("THE OPEN REGIME r > 13:  local drift-absorbed embedding is IMPOSSIBLE")
for r in [F(13), F(14), F(20), F(50)]:
    sp = 1 - 1 / r          # spread/Vmax
    print(f"  r={str(r):>3}: spread/Vmax = {float(sp):.4f} (teeth move {float(sp)*100:.1f}% of a turn per ruler period);"
          f"  need g > {float(r/7):.3f} > 1  IMPOSSIBLE")

print()
print("RESOLUTION OF THE 1/7-vs-2/7 TENSION")
print("  The correct threshold is neither 1/7 nor 2/7: it is  g > (1/7) * (Vmax/Vmin),  ratio-dependent.")
print("  g=1/7 needs r<1 (impossible);  g=2/7 needs r<2;  and g<=1 caps the whole local route at r<7.")
print("  So the local slow-fast bridge dies at r>=7, strictly inside the regime spread13_lonely")
print("  already handles.  hembed must be NON-LOCAL in the open regime r>13.")
print()
print("CONSISTENT with mac-mini-S64's exact counterexample (worst7Struct, Vmax=91, Vmin=9, r=10.1):")
print("  good periods j=5,10,11 have NO lonely phi (exact max_phi minReach = 3/43, 2/31, 1/23 < 1/14),")
print("  while a DISTANT j=25 is lonely (0.2306).  The lonely time is not at the good period.")
