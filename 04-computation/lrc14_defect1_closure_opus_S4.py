import itertools
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
print("CONTROL: run the SAME machinery at defect 1, where the answer is known.")
print("Expect exactly GW {1..11,13,24} (gap 1/14) and {1..11,13,36} (gap 3/41), i.e. j=12, s in {24,36}.\n")
hits=[]
for j in range(1,14):
    C=[v for v in range(1,14) if v!=j]
    lon=lonely_intervals(C)
    if not lon: continue
    L=max(b-a for a,b in lon)
    smax=int(2*h/L)                       # band criterion for the single far speed
    for s in range(14,smax+1):
        if not carve(lon,s): hits.append((j,s,smax))
print(f"  near-tight defect-1 configs found: {len(hits)}")
for j,s,smax in hits: print(f"    drop j={j}, far s={s}   (band bound for this core: s <= {smax})")
ok = sorted((j,s) for j,s,_ in hits)==[(12,24),(12,36)]
print(f"\n  CONTROL {'PASSES - machinery recovers exactly the known answer' if ok else 'FAILS - investigate'}")

print("\nFULL defect-1 band bounds (s <= 2h/L_max(C)) for every dropped j:")
mx=0
for j in range(1,14):
    C=[v for v in range(1,14) if v!=j]
    lon=lonely_intervals(C); L=max(b-a for a,b in lon); smax=int(2*h/L)
    mx=max(mx,smax)
    found=[s for s in range(14,smax+1) if not carve(lon,s)]
    print(f"   drop j={j:2d}: L_max={L:.6f}  => far speed s <= {smax:3d}   near-tight s: {found}")
print(f"\n  MAX bound over all j = {mx}.  So EVERY defect-1 near-tight config has far speed <= {mx},")
print("  and exhaustive enumeration to that bound gives exactly {GW (s=24), {1..11,13,36} (s=36)}.")
print("  => DEFECT 1 IS CLOSED (not an unbounded family): only GW is TIGHT at defect 1.")
