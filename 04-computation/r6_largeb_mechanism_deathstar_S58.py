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
for v in P: S0=sub(S0,danger(v))
S0f=sorted([(float(a),float(b),a,b) for a,b in S0],key=lambda x:-(x[1]-x[0]))
for b in [1000,10000]:
    # find max gap per core arc, record which arc
    best=(0.0,None,None,None)
    for (af,bf,a,bb) in S0f:
        cuts=[]
        for k in [b,b+2,b+4,b+6,b+8]:
            w=1.0/(14*k); jlo=int(af*k)-1; jhi=int(bf*k)+1
            for j in range(jlo,jhi+1):
                c=j/k; lo=c-w; hi=c+w
                if hi>af and lo<bf: cuts.append((max(lo,af),min(hi,bf)))
        cuts.sort(); cur=af
        for lo,hi in cuts:
            if lo>cur and lo-cur>best[0]:
                best=(lo-cur,(af,bf),cur,lo)
            if hi>cur: cur=hi
        if bf-cur>best[0]: best=(bf-cur,(af,bf),cur,bf)
    L=best[0]
    print("b=%d: max gap L=%.7f  L*b=%.5f  in core arc (%.4f,%.4f) width %.4f"%(
        b,L,L*b,best[1][0],best[1][1],best[1][1]-best[1][0]))
    print("   gap spans (%.6f,%.6f), center %.5f, sigma=2t=%.5f (mod1)"%(best[2],best[3],(best[2]+best[3])/2,(2*(best[2]+best[3])/2)%1))
    print("   killer period 1/b=%.6f, gap/period=%.3f, arc/period=%.2f"%(1.0/b,L*b,(best[1][1]-best[1][0])*b))
