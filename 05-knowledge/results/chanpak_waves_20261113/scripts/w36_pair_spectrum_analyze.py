import sys
from collections import defaultdict
for fn in sys.argv[1:]:
    pts=set(); H=set()
    for line in open(fn):
        h,p=map(int,line.split()); pts.add((h,p)); H.add(h)
    print(f'{fn}: classes -> distinct pairs {len(pts)}, distinct H {len(H)}, max H {max(H)}, max P2 {max(p for _,p in pts)}')
    # box density: lattice points (h odd, P even) in boxes around the bulk
    hs=sorted(h for h,_ in pts); hmed=hs[len(hs)//2]
    byh=defaultdict(set)
    for h,p in pts: byh[h].add(p)
    # for each attained h, the attained P2 values: fraction of even P in [min,max] hit
    fr=[]
    for h,ps in byh.items():
        lo,hi=min(ps),max(ps); slots=(hi-lo)//2+1
        fr.append((h,len(ps),slots,len(ps)/slots))
    fr.sort()
    mid=[x for x in fr if x[2]>=20]
    avg=sum(x[3] for x in mid)/max(len(mid),1)
    print(f'   per-H fibre fill (even P2 values hit / slots in [min,max]), fibres with >= 20 slots: mean {avg:.4f}; best {max((x[3] for x in mid),default=0):.4f}')
    # lines alpha*H + beta*P2 for small kernels: coverage of odd N in [lo, hi] for (alpha,beta) in a few cases
    for a,b in [(3,1),(5,3),(1,1),(7,2)]:
        vals=set(a*h+b*p for h,p in pts)
        ov=sorted(v for v in vals if v%2==1)
        if not ov: continue
        # longest run of consecutive odd values from the minimum
        s=set(ov); start=ov[0]; run=0; v=start
        while v in s: run+=1; v+=2
        lo=sorted(ov)[len(ov)//10]; hi=sorted(ov)[len(ov)*9//10]
        inwin=sum(1 for x in ov if lo<=x<=hi); slots=(hi-lo)//2+1
        print(f'   line {a}H+{b}P2: {len(vals)} values; odd fill of central window [{lo},{hi}]: {inwin/slots:.4f}')
