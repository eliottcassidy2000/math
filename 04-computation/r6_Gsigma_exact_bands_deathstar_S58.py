from fractions import Fraction as F
def G(s):
    w=F(1,14); arcs=[]
    for m in range(5):
        c=(m*s)%1; lo=(c-w)%1; hi=(c+w)%1
        if lo<hi: arcs.append((lo,hi))
        else: arcs.append((lo,F(1))); arcs.append((F(0),hi))
    arcs.sort(); merged=[]
    for lo,hi in arcs:
        if merged and lo<=merged[-1][1]: merged[-1]=(merged[-1][0],max(merged[-1][1],hi))
        else: merged.append((lo,hi))
    mg=F(0)
    for i in range(len(merged)):
        a=merged[i][1]; b=merged[(i+1)%len(merged)][0]+(1 if i+1==len(merged) else 0)
        if b-a>mg: mg=b-a
    return mg
# candidate breakpoints: sigma = (k +/- 1/7)/d, d=1..4  -- where arc endpoints coincide
bp=set([F(0),F(1)])
for d in range(1,5):
    for k in range(0,5):
        for e in (F(1,7),-F(1,7),F(0)):
            v=F(k,d)+F(e,d)
            if 0<=v<=1: bp.add(v)
bp=sorted(bp)
# G is piecewise-linear between breakpoints; evaluate at midpoints & endpoints, test >1/7
thr=F(1,7); bands=[]
for i in range(len(bp)-1):
    a,b=bp[i],bp[i+1]; mid=(a+b)/2
    Ga,Gm,Gb=G(a),G(mid),G(b)
    # >1/7 region within [a,b]: since piecewise linear, check where it exceeds; report if midpoint qualifies
    if Gm>thr: bands.append((a,b,Gm))
print("EXACT sigma-bands where G(sigma)>1/7 (piecewise, by breakpoint intervals):")
# merge adjacent qualifying intervals
merged=[]
for a,b,gm in bands:
    if merged and a==merged[-1][1]: merged[-1]=(merged[-1][0],b)
    else: merged.append([a,b])
for a,b in merged:
    print("  (%s, %s) = (%.5f, %.5f)"%(a,b,float(a),float(b)))
print("min over all breakpoint-midpoints:",min(float(G((bp[i]+bp[i+1])/2)) for i in range(len(bp)-1)))
print("G(1/5) =",G(F(1,5)),"(the even-spread minimum)")
print("G(1/3) =",G(F(1,3)),"= %.5f (only 3 distinct arcs, big margin)"%float(G(F(1,3))))
print()
# --- verify worst core [1,2,4,7,9,11,12]: do its safe arcs (x2) hit a good band? ---
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
S=[(F(0),F(1))]
for v in P: S=sub(S,danger(v))
goodbands=merged
hits=0
for lo,hi in sorted(S,key=lambda x:-(x[1]-x[0]))[:8]:
    sig_lo=(2*lo)%1; sig_hi=(2*hi)%1
    # does 2*[lo,hi] intersect a good band (mod 1)? width in sigma = 2*(hi-lo)
    inband=any(not(2*hi<=a or 2*lo>=b) for a,b in goodbands) or any(not((2*hi-1)<=a or (2*lo-1)>=b) for a,b in goodbands)
    if inband: hits+=1
    print("  core arc (%.4f,%.4f) w=%.4f -> sigma=2t in (%.4f,%.4f) hits good band: %s"%(float(lo),float(hi),float(hi-lo),float(sig_lo),float(sig_hi),inband))
print("=> widest core arcs whose 2t lands in a good G-band:",hits)
