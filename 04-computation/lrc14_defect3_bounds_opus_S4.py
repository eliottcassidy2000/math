import numpy as np, itertools, time
h=3/41.0
def Lmax(V):
    """largest gap of the complement of union D_v (floats, fast)."""
    segs=[]
    for v in V:
        m=np.arange(0,v+1)
        lo=(m-h)/v; hi=(m+h)/v
        segs.append(np.stack([lo,hi],axis=1))
    A=np.concatenate(segs,axis=0)
    out=[]
    for lo,hi in A:
        if lo<0: out.append((0.0,hi)); out.append((lo%1.0,1.0))
        elif hi>1: out.append((lo,1.0)); out.append((0.0,hi-1.0))
        else: out.append((lo,hi))
    out=[(a,b) for a,b in out if b>a]
    out.sort()
    merged=[]
    for a,b in out:
        if merged and a<=merged[-1][1]+1e-15: merged[-1]=(merged[-1][0],max(merged[-1][1],b))
        else: merged.append((a,b))
    best=0.0
    for i in range(len(merged)):
        e=merged[i][1]
        nxt=merged[i+1][0] if i+1<len(merged) else merged[0][0]+1.0
        if nxt-e>best: best=nxt-e
    return best
F3=(1-6*h)/(2*h); F2=(1-4*h)/(2*h); TWOH=2*h
print(f"DEFECT-3 closure attempt. factors: k=3 -> {F3:.4f} (23/6), k=2 -> {F2:.4f} (29/6)")
t0=time.time()
s1cap=0; s2cap=0; s3cap=0; worst=None
cores=list(itertools.combinations(range(1,14),3))
for C3 in cores:
    C=[v for v in range(1,14) if v not in C3]
    L=Lmax(C)
    if L<=0: continue
    B3=L*F3; s1max=int(3/B3)
    s1cap=max(s1cap,s1max)
for C3 in cores:
    C=[v for v in range(1,14) if v not in C3]
    L=Lmax(C); B3=L*F3; s1max=int(3/B3)
    for s1 in range(14,min(s1max,400)+1):
        L1=Lmax(C+[s1])
        if L1<=0: continue
        B2=L1*F2; s2max=int(2/B2)
        if s2max<s1: continue
        s2cap=max(s2cap,s2max)
        for s2 in range(s1,min(s2max,400)+1):
            L2=Lmax(C+[s1,s2])
            if L2<=0: continue
            s3max=int(TWOH/L2)
            if s3max<s2: continue
            if s3max>s3cap:
                s3cap=s3max; worst=(C3,s1,s2,L2)
el=time.time()-t0
print(f"  max s1 (smallest far)  <= {s1cap}")
print(f"  max s2 (middle far)    <= {s2cap}")
print(f"  max s3 (largest far)   <= {s3cap}     worst={worst}")
M=max(s1cap,s2cap,s3cap)
print(f"  => ALL three far speeds <= {M}    [{el:.0f}s]")
print(f"  my exhaustive defect-3 scan covered adds <= 55 -> "
      f"{'CLOSED' if M<=55 else 'need scan to '+str(M)}")
