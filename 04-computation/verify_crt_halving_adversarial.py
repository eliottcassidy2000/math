#!/usr/bin/env python3
"""
INDEPENDENT adversarial verification of the CRT/halving claim for LRC(14).
Target gap 1/14. Lonely measure on [0,1) = measure of x where every speed v in E
keeps frac(v x) at distance > 1/14 from 0 (equivalently ||v x|| > 1/14).
We compute the COMPLEMENT: union of bands [k/v - 1/(14v), k/v + 1/(14v)] over all
integer k, intersected with [0,1). Lonely = 1 - measure(union).

This is my own from-scratch implementation. Exact fractions.
"""
from fractions import Fraction as F

def lonely(E):
    """Exact lonely measure on [0,1). x is lonely iff for all v in E, ||v x||>1/14."""
    bands=[]
    for v in E:
        if v==0: continue
        w=F(1,14*v)
        # centers k/v for k=0..v intersected with [0,1]; clamp
        for k in range(0,v+1):
            c=F(k,v)
            lo=c-w; hi=c+w
            if hi<=0 or lo>=1: continue
            bands.append((max(lo,F(0)),min(hi,F(1))))
    bands.sort()
    merged=[]
    for s,e in bands:
        if merged and s<=merged[-1][1]:
            if e>merged[-1][1]: merged[-1]=(merged[-1][0],e)
        else:
            merged.append((s,e))
    cov=sum(e-s for s,e in merged)
    return F(1)-cov

def odd_part(v):
    while v%2==0: v//=2
    return v

# ---------- FACT (1): RIGID ODD BASE ----------
print("="*70)
print("FACT 1: odd speeds of [1,13]\\{e} for even e are always {1,3,5,7,9,11,13}")
print("        and L_odd is constant.")
ODD=[1,3,5,7,9,11,13]
L_odd=lonely(ODD)
print(f"  L(odd base {ODD}) = {L_odd} = {float(L_odd):.8f}")
print(f"  claimed 75454/315315 = {F(75454,315315)} = {float(F(75454,315315)):.8f}")
print(f"  MATCH: {L_odd==F(75454,315315)}")
for e in [2,4,6,8,10,12]:
    core=[v for v in range(1,14) if v!=e]
    odds=sorted(v for v in core if v%2==1)
    print(f"  drop {e:>2}: odd speeds = {odds}  L_odd={float(lonely(odds)):.8f}")

# ---------- FACT (2): DYADIC HALVING CASCADE ----------
print("\n"+"="*70)
print("FACT 2: dyadic halving cascade on CONSEC {1..13}")
print("  odd base -> +level1{2,6,10} -> +level2{4,12} -> +level3{8}")
cur=list(ODD); L=lonely(cur); print(f"  odd base:        L={float(L):.8f}")
for lvl,evs in [(1,[2,6,10]),(2,[4,12]),(3,[8])]:
    cur=sorted(cur+evs); newL=lonely(cur)
    ratio=float(newL/L) if L!=0 else float('nan')
    print(f"  +level{lvl} {evs}: L={float(newL):.8f}  ratio={ratio:.4f}")
    L=newL
print(f"  final (=consec) L={float(L):.8f}  (should be 0)")

# ---------- FACT (3): CRT CLEAN PRODUCT for drop-6 ----------
print("\n"+"="*70)
print("FACT 3: drop-6 CRT clean product L = (2/7)*(49/1716) = 7/858")
C0=[v for v in range(1,14) if v!=6]
Lc0=lonely(C0)
print(f"  L(drop-6 core {C0}) = {Lc0} = {float(Lc0):.8f}")
print(f"  claimed 7/858 = {F(7,858)} = {float(F(7,858)):.8f}")
print(f"  MATCH: {Lc0==F(7,858)}")

# per-cell decomposition a/7 + xi/7
def lonely_in_cell(E,a):
    lo=F(a,7)-F(1,14); hi=F(a,7)+F(1,14)
    bands=[]
    for v in E:
        if v==0: continue
        w=F(1,14*v)
        kmin=int(lo*v)-2; kmax=int(hi*v)+2
        for k in range(kmin,kmax+1):
            c=F(k,v); blo=max(c-w,lo); bhi=min(c+w,hi)
            if blo<bhi: bands.append((blo,bhi))
    bands.sort(); merged=[]
    for s,e in bands:
        if merged and s<=merged[-1][1]:
            if e>merged[-1][1]: merged[-1]=(merged[-1][0],e)
        else: merged.append((s,e))
    return (hi-lo)-sum(e-s for s,e in merged)

cells=[lonely_in_cell(C0,a) for a in range(7)]
print(f"  per-cell masses: {[float(c) for c in cells]}")
print(f"  sum cells = {sum(cells,F(0))}  (should equal L(drop-6)? {sum(cells,F(0))==Lc0})")
surv=[a for a in range(7) if cells[a]>0]
print(f"  surviving cells: {surv}")
Rs=[7*cells[a] for a in surv]
print(f"  R_a = 7*m_a on surviving: {[str(r) for r in Rs]}")
print(f"  all equal (clean product)? {len(set(Rs))==1}")
if len(set(Rs))==1:
    print(f"  product = (#surv/7)*R = ({len(surv)}/7)*{Rs[0]} = {F(len(surv),7)*Rs[0]}")

# NOTE: cell decomposition sums over k in a window, but full [0,1) lonely uses k=0..v.
# Check the cell sum reconstructs the full measure (the 7-cell partition is exhaustive).

# ---------- FACT (4): all 5 even deletions beat all 7 odd deletions ----------
print("\n"+"="*70)
print("FACT 4: ranking all single deletions; do all evens beat all odds?")
rows=[]
for e in range(1,14):
    core=[v for v in range(1,14) if v!=e]
    rows.append((lonely(core),e))
rows.sort()
even_Ls=[float(L) for L,e in rows if e%2==0]
odd_Ls=[float(L) for L,e in rows if e%2==1]
for L,e in rows:
    par='EVEN' if e%2==0 else 'ODD'
    print(f"  drop {e:>2} [{par}]: L={float(L):.8f}")
print(f"  max even-deletion L = {max(even_Ls):.8f}")
print(f"  min odd-deletion  L = {min(odd_Ls):.8f}")
print(f"  ALL evens < ALL odds? {max(even_Ls) < min(odd_Ls)}")

# ---------- ADVERSARIAL: is drop-6 really the UNIQUE minimizer? ----------
print("\n"+"="*70)
print("ADVERSARIAL: unique minimizer among single deletions?")
minL=rows[0][0]
mins=[e for L,e in rows if L==minL]
print(f"  min L = {minL} = {float(minL):.8f}, achieved by deletions: {mins}")

# ---------- ADVERSARIAL: CRT product on OTHER cores ----------
print("\n"+"="*70)
print("ADVERSARIAL: does the clean R-constant product hold for ALL cores or just some?")
for e in range(1,14):
    core=[v for v in range(1,14) if v!=e]
    cells=[lonely_in_cell(core,a) for a in range(7)]
    surv=[a for a in range(7) if cells[a]>0]
    Rs=[7*cells[a] for a in surv]
    clean=(len(set(Rs))<=1)
    total=sum(cells,F(0))
    recon = (F(len(surv),7)*Rs[0]) if (clean and surv) else None
    ok = (recon==total) if recon is not None else False
    print(f"  drop {e:>2}: #surv={len(surv)} R-constant={clean} cellsum={float(total):.6f} "
          f"{'product=cellsum:'+str(ok) if clean else ''}")
