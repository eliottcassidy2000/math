from pathlib import Path
from collections import deque
import json

def check(c,m):
    if not c: raise RuntimeError(m)

base=Path(__file__).resolve().parent
bank={tuple(map(int,x.split())) for x in (base/'bank6.txt').read_text().splitlines()}
base67={tuple(map(int,x.split())) for x in (base.parent/'global-induction'/'bank67.txt').read_text().splitlines()}
check(base67 <= bank,'bank67 not represented in independent order-six bank')
for m in range(1,21):
    for j in range(m):
        bank.add((1,int(j==0),2**j,int(j==m-1),2**(m-1-j),2**(m-1),int(m==1),1))
bank=sorted(bank)

def next_states(v):
    h,e,r,d=v
    # Source singleton represented by the independently generated bank.
    yield h+r-e,d,2*r-d,d
    yield h,0,2*r,d
    for g,s,l,E,q,b,delta,dg in bank:
        yield g*h-s*e+l*r,E*h-delta*e,q*h+b*r,dg*h

def close(h,d,mul):
    cap=mul*h
    seed=(h,d,h,d)
    queue=deque([seed]); seen={seed}; values={h}
    while queue:
        for out in next_states(queue.popleft()):
            H,E,R,D=out
            check(0<=E<=D<=H<=R and R+E>=H+D,('nonphysical',seed,out))
            if H<=cap: values.add(H)
            if R<=cap and out not in seen:
                seen.add(out);queue.append(out)
    return values,len(seen)

physical=[(3,1),(5,1),(9,1),(9,3),(15,3),(15,5),(17,1),(17,5),(19,5),(23,5),(27,3),(37,5),(49,13)]
results=[]
for h,d in physical:
    vals,states=close(h,d,5)
    want={h,2*h-d,3*h,4*h-3*d,5*h}
    check(vals==want,('physical failure',h,d,vals,want))
    results.append({'h':h,'d':d,'values':sorted(vals),'retained_states':states})
count=0
for d in range(1,102,2):
    for gap in range(0,202,2):
        h=3*d+gap
        vals,states=close(h,d,5)
        check(vals=={h,2*h-d,3*h,4*h-3*d,5*h},('parameter failure',h,d,vals))
        count+=1
for h in (1,3,5,9,11,13,15,101,10**50+1):
    vals,states=close(h,h,5)
    check(vals=={h,3*h,5*h},('boundary failure',h,vals))
print(json.dumps({'status':'PASS','independent_gadget_signatures_through_six':345,'signature_bank_with_transitives_through_twenty':len(bank),'original_bank67_contained':True,'parameter_pairs_checked':count,'boundary_pairs_checked':9,'physical_seed_results':results},indent=2))
