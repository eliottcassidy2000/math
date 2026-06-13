"""
monad-explorer-2026-06-07-S4  (exploration of HYP-2298 Q1)
=========================================================
The Moser-angle LADDER.  Each unit direction

    w_t = exp(i arccos(1 - 1/2t)) = ((2t-1) + i sqrt(4t-1)) / (2t),   |w_t| = 1,

has CM discriminant sqrt(4t-1):  t=1 sqrt3, t=2 sqrt7, t=3 sqrt11, t=4 sqrt15, ...
The Moser lattice (Engel et al.) is the rank-4 lattice L_t = Z[zeta6, w_t] at t=3
(zeta6 = w_1, the 60-deg generator; w_3 = the Moser-spindle closure angle cos=5/6).

QUESTION (HYP-2298): how many UNIT VECTORS does L_t have for each t?  Is t=3's 18
special?  EXACT integer arithmetic: scale every coordinate by S=4t to clear
denominators, work in Z[1, sqrt3, sqrt(4t-1), sqrt(3(4t-1))]; unit vector iff
|S z|^2 == (S^2, 0, 0, 0).
"""
from itertools import product
from math import isqrt

def is_square(n):
    if n < 0: return False
    r = isqrt(n); return r*r == n

def count_units(t, R):
    r1, r2 = 3, 4*t-1
    r12 = r1*r2
    S = 4*t
    SS = S*S
    # integer mul table for basis 0:1 1:sqrt(r1) 2:sqrt(r2) 3:sqrt(r12)
    MT = {
      (0,0):(1,0),(0,1):(1,1),(0,2):(1,2),(0,3):(1,3),
      (1,1):(r1,0),(1,2):(1,3),(1,3):(r1,2),
      (2,2):(r2,0),(2,3):(r2,1),
      (3,3):(r12,0),
    }
    def mul(x,y):
        r=[0,0,0,0]
        for i in range(4):
            xi=x[i]
            if xi==0: continue
            for j in range(4):
                yj=y[j]
                if yj==0: continue
                a,b=(i,j) if i<=j else (j,i)
                c,idx=MT[(a,b)]; r[idx]+=xi*yj*c
        return r
    # S * (basis complex coords)  (re, im) integer 4-tuples:
    m=2*t-1
    RE=[(S,0,0,0),(2*t,0,0,0),(2*m,0,0,0),(m,0,0,-1)]      # 1, w1, w_t, w1 w_t
    IM=[(0,0,0,0),(0,2*t,0,0),(0,0,2,0),(0,m,1,0)]
    def coord(v):
        re=[0,0,0,0]; im=[0,0,0,0]
        for k,(rr,ii) in zip(v,zip(RE,IM)):
            if k==0: continue
            for q in range(4):
                re[q]+=k*rr[q]; im[q]+=k*ii[q]
        return re,im
    TARGET=[SS,0,0,0]
    units=[]
    for v in product(range(-R,R+1),repeat=4):
        if v==(0,0,0,0): continue
        re,im=coord(v)
        n=[a+b for a,b in zip(mul(re,re),mul(im,im))]
        if n==TARGET: units.append(v)
    return units

print("="*74)
print("MOSER-ANGLE LADDER: unit-vector counts of L_t = Z[zeta6, w_t]")
print("  w_t = ((2t-1)+i sqrt(4t-1))/(2t),  disc = sqrt(4t-1).  EXACT integer arithmetic.")
print("="*74)
print(f"  {'t':>2} {'4t-1':>5} {'disc':>8} {'#units':>7} {'stable':>7}  {'transv':>6}   notes")
summary=[]
for t in range(1,15):
    r2=4*t-1
    degen = is_square(r2) or is_square(3*r2) or (r2 % 3 == 0 and is_square(r2//3))
    if degen:
        note=("rank-2: = Z[zeta6] (t=1 triangular rosette, 6 units)" if t==1
              else f"rank-2 collapse: 4t-1={r2} keeps w_t in Q(sqrt-3) [count N/A in 4D]")
        print(f"  {t:>2} {r2:>5} {'sqrt'+str(r2):>8} {'rank2':>7} {'--':>7}  {'--':>6}   {note}")
        summary.append((t,r2,None,note)); continue
    u5=count_units(t,5); u8=count_units(t,8)
    stable = len(u5)==len(u8)
    cnt=len(u8)
    transv=sum(1 for u in u8 if u[2]!=0 or u[3]!=0)
    note = ("<== MOSER lattice (Engel et al.); spindle cos=5/6" if t==3 else
            "sqrt7 rung (THM-421/HYP-2262 radius); Q(sqrt-3,sqrt-7)" if t==2 else "")
    print(f"  {t:>2} {r2:>5} {'sqrt'+str(r2):>8} {cnt:>7} {('yes' if stable else 'NO!'):>7}  {transv:>6}   {note}")
    summary.append((t,r2,cnt,note))

print()
print("READING:")
print(" - t=1: rank-2 Z[zeta6], the triangular rosette = 6 unit vectors.")
print(" - rank-4 t: EXACT unit-vector count (box-stable R=5 vs R=8).")
print("   18 = 6 triangular + 12 transverse w_t-units;  12 = 6 triangular + 6 transverse.")
print(" - which discriminants give 18 vs 12:")
c18=[r2 for (t,r2,c,n) in summary if c==18]; c12=[r2 for (t,r2,c,n) in summary if c==12]
cother=[(r2,c) for (t,r2,c,n) in summary if c not in (12,18,None)]
print(f"     18-unit (4t-1): {c18}")
print(f"     12-unit (4t-1): {c12}")
if cother: print(f"     other counts  : {cother}")
print(" - degenerate t (4t-1 square or 3*square) collapse into Q(sqrt-3); not rank-4.")
print("   (4t-1=27 at t=7 => sqrt27=3 sqrt3 => w_7 in Q(sqrt-3): triangular-compatible angle.)")
