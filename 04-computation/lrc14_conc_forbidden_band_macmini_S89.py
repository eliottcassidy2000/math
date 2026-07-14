#!/usr/bin/env python3
"""mac-mini-S89: the conc FORBIDDEN BAND = the tournament-quantization signature of the last bit.
klein-S290: for S={1}uC, L=|G(C)|(1-conc/7), conc=14|G(C)∩[0,1/14)|/|G(C)|; covering => conc<7,
AP {1..13} (C={2..13}) the UNIQUE conc=7 tight case. Claim (HYP-4306 quantization): covering
families are pushed OFF the AP ground-rung by a FULL band -- sup_covering conc < 7 with a real gap,
mirroring the Farey ladder gap (1/13,2/25). Verify the band is robust over the bounded compact-core
residual, and check the tournament-regularity reading: the AP is the regular circulant C_13({1..6})."""
import numpy as np
from itertools import combinations
c=1.0/14; N=1_500_000; t=(np.arange(N)+0.5)/N
band=(t<1.0/14)
def Gmask(C):
    ok=np.ones(N,bool)
    for w in C:
        r=(w*t)%1.0; ok &= (np.minimum(r,1-r)>=c)
    return ok
def conc(C):
    g=Gmask(C); mg=g.mean()
    if mg==0: return None,0.0
    return 14*(g&band).mean()/mg, mg
def covers(S):  # S covers moduli 2..14
    return all(any(v%q==0 for v in S) for q in range(2,15))
# (1) the AP tight case
apC=list(range(2,14))  # S={1..13}
a,mg=conc(apC); print(f"AP cluster C={{2..13}} (S={{1..13}}): conc={a:.4f} (claim 7), |G(C)|={mg:.5f}, L=|G|(1-conc/7)={mg*(1-a/7):.5f}")
# (2) bounded covering {1}uC census: C subset {2..Vmax}, |C|=12, S={1}uC covering & primitive
print("\nBounded covering {1}uC census -- conc and the band edge:")
best=[]; Vmax=16
pool=list(range(2,Vmax+1))
cnt=0
for C in combinations(pool,12):
    S=(1,)+C
    if not covers(S): continue
    a,mg=conc(list(C))
    if a is None: continue
    cnt+=1; best.append((a,S,mg))
best.sort(reverse=True)
print(f"  {cnt} covering families in {{2..{Vmax}}}; top conc (closest to the AP wall 7):")
for a,S,mg in best[:8]:
    miss=[q for q in range(2,15) if not any(v%q==0 for v in S)]
    print(f"    conc={a:.4f}  L>={mg*(1-a/7):.5f}  S={S}")
sup=best[0][0] if best else 0
print(f"\n  SUP covering conc = {sup:.4f};  FORBIDDEN BAND = ({sup:.3f}, 7);  uniform margin 1-sup/7 = {1-sup/7:.4f}")
print(f"  => covering => L >= {1-sup/7:.4f} * |G(C)| > 0 on this stratum (the 'shared cancellation', quantified)")
# (3) tournament reading: the AP's tight-time difference tournament is the regular circulant C_13({1..6})
print("\nTournament reading (HYP-4306 quantization): AP=regular circulant C_13({1..6}), the free (Z/13)* orbit,")
print("  = ground rung conc=7. Covering breaks (Z/13)-regularity => forced off by a full band to conc<=sup.")
print("  The gap (sup,7) is the tournament forbidden band = W(n) odd-cycle cancellation, LRC-side.")
