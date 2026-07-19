from fractions import Fraction as F
P=[1,2,4,7,9,11,12]
def danger(v):
    w=F(1,14*v); out=[]
    for j in range(v):
        c=F(j,v); lo=(c-w)%1; hi=(c+w)%1
        if lo<hi: out.append((lo,hi))
        else: out.append((lo,F(1))); out.append((F(0),hi))
    return out
def sub(safe,arcs):
    for clo,chi in sorted(arcs):
        new=[]
        for lo,hi in safe:
            if chi<=lo or clo>=hi: new.append((lo,hi)); continue
            if clo>lo: new.append((lo,clo))
            if chi<hi: new.append((chi,hi))
        safe=new
    return safe
S0=[(F(0),F(1))]
for v in P: S0=sub(S0,danger(v))  # core-safe arcs (exact)
def Lb(b):  # max gap within core arcs only (gap is always inside a core arc)
    best=F(0)
    for (a,bb) in S0:
        cuts=[]
        for k in [b,b+2,b+4,b+6,b+8]:
            w=F(1,14*k); jlo=int(a*k)-1; jhi=int(bb*k)+2
            for j in range(jlo,jhi+1):
                c=F(j,k); lo=c-w; hi=c+w
                if hi>a and lo<bb: cuts.append((max(lo,a),min(hi,bb)))
        cuts.sort(); cur=a
        for lo,hi in cuts:
            if lo>cur and lo-cur>best: best=lo-cur
            if hi>cur: cur=hi
        if bb-cur>best: best=bb-cur
    return best
worst=F(0); wb=None; allok=True
for b in range(157,400):
    L=Lb(b); R=F(1,7*(b+8))/L
    if R>=1: allok=False; print("FAIL b=%d R=%.5f"%(b,float(R)))
    if R>worst: worst=R; wb=b
print("EXACT (arc-restricted) finite check b in [157,399]:")
print("  ALL R_sharp<1:",allok,"; max R_sharp=%.6f = %s at b=%d"%(float(worst),worst,wb))
