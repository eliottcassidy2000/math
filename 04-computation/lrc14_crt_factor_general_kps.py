#!/usr/bin/env python3
"""Is the (2/7)-type clean 7-adic x archimedean factorization SPECIAL to drop-6, or general?
Check #surviving cells and R_a-constancy for ALL single-deletion cores. EXACT."""
import sys
from fractions import Fraction as F
sys.stdout.reconfigure(line_buffering=True) if hasattr(sys.stdout,'reconfigure') else None
def cell_mass(E,a):
    lo=F(a,7)-F(1,14); hi=F(a,7)+F(1,14); bands=[]
    for v in E:
        if v==0: continue
        w=F(1,14*v); 
        for k in range(int(lo*v)-2,int(hi*v)+3):
            c=F(k,v); blo=max(c-w,lo); bhi=min(c+w,hi)
            if blo<bhi: bands.append((blo,bhi))
    bands.sort(); merged=[]
    for s,e in bands:
        if merged and s<=merged[-1][1]: merged[-1]=(merged[-1][0],max(merged[-1][1],e))
        else: merged.append((s,e))
    return (hi-lo)-sum(e-s for s,e in merged)
print(f"{'drop':>4} {'L':>14} {'#surv':>6} {'surv cells':>16} {'R-constant?':>11}")
for e in range(1,14):
    core=[v for v in range(1,14) if v!=e]
    cells=[cell_mass(core,a) for a in range(7)]
    surv=[a for a in range(7) if cells[a]>0]
    Rs=[7*cells[a] for a in surv]
    const = (len(set(Rs))==1)
    L=sum(cells,F(0))
    print(f"{e:>4} {str(L):>14} {len(surv):>6} {str(surv):>16} {str(const):>11}")
print("\n=> '#surv cells' is the 7-adic factor; R-constant=True means CLEAN CRT product L=(#surv/7)*R.")
