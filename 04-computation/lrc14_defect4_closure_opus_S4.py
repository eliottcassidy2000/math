import itertools, time, sys
h=3/41.0
def lonely_intervals(V):
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
    out=[]
    for i in range(len(mg)):
        e=mg[i][1]; nxt=mg[i+1][0] if i+1<len(mg) else mg[0][0]+1.0
        if nxt-e>1e-13: out.append((e,min(nxt,1.0)))
    return [(a,b) for a,b in out if b>a]
def carve(lon,s):
    out=[]
    for a,b in lon:
        cur=a; m0=int(a*s)-1; m1=int(b*s)+2
        for m in range(m0,m1+1):
            lo=(m-h)/s; hi=(m+h)/s
            if hi<=cur or lo>=b: continue
            if lo>cur: out.append((cur,min(lo,b)))
            if hi>cur: cur=hi
            if cur>=b: break
        if cur<b: out.append((cur,b))
    return [(a,b) for a,b in out if b-a>1e-13]
NODES=0
def search(lon,lo,krem,chosen):
    global NODES
    if not lon: return list(chosen)
    if krem==0: return None
    L=max(b-a for a,b in lon)
    if krem==1:
        smax=int(2*h/L)                      # band criterion: one speed must cover all
    else:
        R=L*(1-2*krem*h)/(2*h)
        if R<=0: return None
        smax=int(krem/R)                     # klein lemma: krem speeds all >= s
    if smax<lo: return None
    for s in range(lo,smax+1):
        NODES+=1
        r=search(carve(lon,s),s+1,krem-1,chosen+[s])
        if r: return r
    return None
K=4
t0=time.time(); found=[]
cores=list(itertools.combinations(range(1,14),K))
for idx,drop in enumerate(cores):
    C=[v for v in range(1,14) if v not in drop]
    lon=lonely_intervals(C)
    if not lon: continue
    r=search(lon,14,K,[])
    if r: found.append((drop,r))
    if time.time()-t0>500: print(f"  TIME CAP at core {idx+1}/{len(cores)}"); break
el=time.time()-t0
print(f"DEFECT-4 closure search (klein-lemma pruning at every node + band criterion at the last)")
print(f"  9-cores processed: {idx+1}/{len(cores)}   nodes={NODES:,}   [{el:.0f}s]")
print(f"  near-tight defect-4 configs found: {len(found)}")
for d,r in found[:10]: print("    drop",d,"far",r)
print(f"  => {'DEFECT-4 CLOSED (no config has gap<=3/41)' if not found and idx+1==len(cores) else 'incomplete or hits found'}")
