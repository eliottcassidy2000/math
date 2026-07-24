from fractions import Fraction as Fr
import itertools
h=Fr(3,41)
def lonely_arcs(V,h):
    arcs=[]
    for v in V:
        for m in range(0,v+1):
            lo=(Fr(m)-h)/v; hi=(Fr(m)+h)/v
            arcs.append((lo,hi))
    norm=[]
    for lo,hi in arcs:
        if lo<0: norm.append((Fr(0),hi)); norm.append((lo%1,Fr(1)))
        elif hi>1: norm.append((lo,Fr(1))); norm.append((Fr(0),hi-1))
        else: norm.append((lo,hi))
    norm=[(a,b) for a,b in norm if a<b]; norm.sort()
    merged=[]
    for a,b in norm:
        if merged and a<=merged[-1][1]: merged[-1]=(merged[-1][0],max(merged[-1][1],b))
        else: merged.append((a,b))
    gaps=[]
    for i in range(len(merged)):
        e=merged[i][1]
        nxt=merged[i+1][0] if i+1<len(merged) else merged[0][0]+1
        if nxt>e: gaps.append(nxt-e)
    return max(gaps) if gaps else Fr(0)
TWOH=2*h
print("DEFECT-2 CLOSURE TEST at h=3/41.")
print("Criterion (exact): gap(C u {s,r}) <= h  =>  r <= 2h / L_max(C u {s})")
print("  because D_r has bands of width 2h/r separated by gaps, so an interval longer")
print("  than 2h/r cannot fit inside a band.  If all r_max <= 300, my exhaustive")
print("  adds<=300 scan already covered every possibility => DEFECT-2 CLOSED.\n")
worst=[]; zeros=0
for i,j in itertools.combinations(range(1,14),2):
    C=[v for v in range(1,14) if v not in (i,j)]
    for s in range(14,71):                 # klein's lemma: min(far) <= 70
        L=lonely_arcs(C+[s],h)
        if L==0: zeros+=1; continue
        rmax=int(TWOH/L)
        worst.append((rmax,(i,j),s,float(L)))
worst.sort(reverse=True)
print(f"  (C,s) pairs tested: {len(worst)+zeros}   (L_max=0 cases: {zeros})")
print("  LARGEST r_max values:")
for rmax,(i,j),s,L in worst[:12]:
    print(f"    drop({i:2d},{j:2d}) s={s:3d}: L_max(Cu{{s}})={L:.6f}  =>  r <= {rmax}")
MX=worst[0][0]
print(f"\n  MAX r_max over ALL (core, s<=70) = {MX}")
print(f"  my exhaustive scan covered adds <= 300  ->  {'CLOSED: every defect-2 near-tight config was scanned' if MX<=300 else 'RESIDUAL remains for r in (300, %d]'%MX}")
if MX>300:
    rem=sum(1 for rmax,_,_,_ in worst if rmax>300)
    tot=sum(max(0,rmax-300) for rmax,_,_,_ in worst if rmax>300)
    print(f"  pairs with r_max>300: {rem}; total unscanned (s,r) configs = {tot:,}  <- finite, scannable")
