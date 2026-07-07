from fractions import Fraction as F
def distZ(num,den):
    r=num%den; return F(min(r,den-r),den)
def reach_M_wit(V):
    V=[abs(v) for v in V if v!=0]; n=len(V); dens=set()
    for i in range(n):
        dens.add(2*V[i])
        for j in range(i+1,n):
            dens.add(V[i]+V[j])
            if V[i]!=V[j]: dens.add(abs(V[i]-V[j]))
    best=F(0); bt=None
    for d in dens:
        if d==0: continue
        for m in range(1,d):
            mn=min(distZ(v*m,d) for v in V)
            if mn>best: best=mn; bt=F(m,d)
    return best,bt
def gcount(V,t):  # # distinct gap lengths of {0} U {v_i t mod 1} on the circle
    pts=sorted(set([F(0)]+[ (F(v)*t) - int(F(v)*t) + (1 if (F(v)*t)-int(F(v)*t)<0 else 0) for v in V]))
    pts=sorted(set((F(v)*t)%1 for v in V) | {F(0)})
    gaps=[]
    for a,b in zip(pts,pts[1:]): gaps.append(b-a)
    gaps.append(F(1)-pts[-1]+pts[0])
    return len(set(gaps)), sorted(set(gaps))

print("=== LRC(14) three-gap rigidity g-count (extends mac-mini-S15 from 12-speed to 13-speed), opus-S132 ===")
print("    g = # distinct gap lengths of phases at the M-witness. S15: near-tight => g<=3, loose => g>=7\n")
fams = {
  "AP {1..13} (tight)": list(range(1,14)),
  "GW {1..11,13,24} (2nd tight)": [1,2,3,4,5,6,7,8,9,10,11,13,24],
  "2*AP {2..26 even} (sat tight)": [2*i for i in range(1,14)],
  "{1..12,182} (deep well, sat)": list(range(1,13))+[182],
  "{1..4,10..18} (sat, min-M ss)": [1,2,3,4,10,11,12,13,14,15,16,17,18],
  "{2..14} consec (sat, loose)": list(range(2,15)),
  "{1,2,3,4,5,7..14} (sat)": [1,2,3,4,5,7,8,9,10,11,12,13,14],
  "generic loose {3,5,7,11,13,17,19,23,29,31,37,41,43}": [3,5,7,11,13,17,19,23,29,31,37,41,43],
}
print(f"{'family':<44} {'M':>9} {'~M':>7} {'g':>3}")
for name,V in fams.items():
    M,t=reach_M_wit(V); g,gaps=gcount(V,t)
    print(f"  {name:<44} {str(M):>9} {float(M):>7.4f} {g:>3}   t={t}")
print("\n  => if near-tight LRC14 families have g<=3-4 and loose g>=7, S15's three-gap rigidity")
print("     EXTENDS to 13 speeds; the moat/(G) = converse-three-gap: M near 1/14 => g<=3 => {k*alpha} orbit => M on Farey ladder.")
