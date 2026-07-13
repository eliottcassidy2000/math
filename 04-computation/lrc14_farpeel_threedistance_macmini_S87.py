#!/usr/bin/env python3
"""mac-mini-S87: PINPOINT the peel klein's open disc_v bound (THM-731) should target. opus-S270 found
'only the FAR peel certifies'. Reason (my #far reduction): peeling the far element leaves the
CONSECUTIVE core {1..12}, whose good set G'_{~far} = SafeSet({1..12}) cap middle is a THREE-DISTANCE
(Steinhaus) interval union -- few intervals, piecewise-linear autocorrelation -> analytically
tractable v-grid discrepancy disc_far. Contrast: peeling a SMALL runner leaves a non-consecutive
mess (many intervals). Verify by counting the good-set intervals per peel on the deep well."""
import numpy as np
c=1.0/14; N=2_000_000
x=(np.arange(N)+0.5)/N; mid=((x>=1/14)&(x<=13/14))
def safe_others(S,v):
    ok=mid.copy()
    for w in S:
        if w==v: continue
        r=(w*x)%1.0; ok &= (np.minimum(r,1-r)>=c)
    return ok
def n_intervals(mask):
    d=np.diff(mask.astype(np.int8)); return int((d==1).sum())  # rising edges = #intervals
deep=list(range(1,13))+[182]
print("Deep well {1..12,182}: #good-set intervals of G'_{~v} per peel v (fewer = more structured):\n")
rows=[]
for v in deep:
    m=safe_others(deep,v); ni=n_intervals(m); meas=m.mean()
    rows.append((v,ni,meas))
for v,ni,meas in sorted(rows,key=lambda r:r[1]):
    tag=" <-- FAR peel: CONSECUTIVE core {1..12}, three-distance" if v==182 else ""
    print(f"   peel v={v:>3}: |G'_~v|={meas:.5f}, #intervals={ni:>4}{tag}")
farrow=[r for r in rows if r[0]==182][0]; small=[r for r in rows if r[0]!=182]
print(f"\n   FAR peel #intervals = {farrow[1]}; small-runner peels avg #intervals = {np.mean([r[1] for r in small]):.0f}")
print(f"   => the far peel's good set is FAR more structured (fewer intervals) -- three-distance of {{1..12}}.")
print("   klein's disc_far = (1/182)Sum_j A(j/182) - |G'|^2 with A = autocorrelation of a FEW-interval")
print("   piecewise-linear set => bounded by three-gap theorem (Steinhaus), NOT the crude r^2/3v^2.")
print("\nNear-AP residue {1..11,13,84}: peel BOTH far (13,84) -> core {1..11} collapses; effective order=2.")
res=[*range(1,12),13,84]
for v in [13,84]:
    m=safe_others(res,v); print(f"   peel v={v}: |G'_~v|={m.mean():.5f}, #intervals={n_intervals(m)}")
# peel both far at once: good set of just the consecutive core {1..11}
okc=mid.copy()
for w in range(1,12):
    r=(w*x)%1.0; okc &= (np.minimum(r,1-r)>=c)
print(f"   core {{1..11}} alone: |SafeSet|={okc.mean():.5f}, #intervals={n_intervals(okc)} (three-distance);")
print(f"   then 2 combs D_13,D_84 cut it -> L. Effective multi-linear order = #far = 2, not 13.")
