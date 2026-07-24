import itertools, time
h=3/41.0
def lonely(V):
    segs=[]
    for v in V:
        for m in range(0,v+1):
            lo=(m-h)/v; hi=(m+h)/v
            if lo<0: segs.append((0.0,hi)); segs.append((lo+1.0,1.0))
            elif hi>1: segs.append((lo,1.0)); segs.append((0.0,hi-1.0))
            else: segs.append((lo,hi))
    segs=[s for s in segs if s[1]>s[0]]; segs.sort()
    mg=[]
    for a,b in segs:
        if mg and a<=mg[-1][1]+1e-15: mg[-1]=(mg[-1][0],max(mg[-1][1],b))
        else: mg.append((a,b))
    gaps=[]
    for i in range(len(mg)):
        e=mg[i][1]
        nxt=mg[i+1][0] if i+1<len(mg) else mg[0][0]+1.0
        if nxt-e>1e-12: gaps.append((e,nxt))
    return gaps
def covered(a,b,S):
    """is [a,b] covered by union of D_s for s in S ?"""
    bands=[]
    for s in S:
        m0=int((a*s)-1); m1=int((b*s)+2)
        for m in range(m0,m1+1):
            lo=(m-h)/s; hi=(m+h)/s
            if hi>a and lo<b: bands.append((lo,hi))
    if not bands: return False
    bands.sort(); cur=a
    for lo,hi in bands:
        if lo>cur+1e-12: return False
        if hi>cur: cur=hi
        if cur>=b-1e-12: return True
    return cur>=b-1e-12
t0=time.time(); n=0; hits=[]
cores=list(itertools.combinations(range(1,14),3))
for C3 in cores:
    C=[v for v in range(1,14) if v not in C3]
    G=lonely(C)
    G.sort(key=lambda x:-(x[1]-x[0]))     # longest first -> fastest rejection
    if not G: continue
    for s1,s2,s3 in itertools.combinations(range(14,83),3):
        n+=1
        S=(s1,s2,s3); ok=True
        for (a,b) in G:
            if not covered(a,b,S): ok=False; break
        if ok: hits.append((C3,s1,s2,s3))
el=time.time()-t0
print(f"DEFECT-3 EXHAUSTIVE over the PROVED region (all far speeds <= 82):")
print(f"  cores={len(cores)}  triples/core={len(list(itertools.combinations(range(14,83),3)))}  scanned={n:,}  [{el:.0f}s]")
print(f"  near-tight (gap<=3/41) found: {len(hits)}")
for x in hits[:10]: print("   ",x)
print(f"\n  => {'DEFECT-3 CLOSED: no defect-3 config has gap<=3/41' if not hits else 'HITS FOUND - investigate'}")
