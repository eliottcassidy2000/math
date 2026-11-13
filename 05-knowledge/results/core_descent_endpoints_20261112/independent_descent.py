"""Independent one-path gap compiler, plus literal endpoint controls.

The producer convolves two-path covers. This audit instead counts Hamilton
paths in the actual two-vertex tail and uses H_tail=D+2U.
"""
from context import *
from itertools import combinations
from collections import defaultdict,Counter
from math import comb
from synchrony_probe import expand
from tensor import kernel as signed_kernel
import time

def gaps_hamilton(rows,ports):
    fixed=[v for v in range(len(rows)) if v not in ports];n=len(fixed);base=n+2
    ordered=sorted(ports,key=lambda v:-sum(rows[v]>>w&1 for w in ports))
    need(all(rows[u]>>v&1 for i,u in enumerate(ordered) for v in ordered[i+1:]),'acyclic gap alphabet')
    chains=[tuple(ordered[i] for i in range(3) if mask>>i&1) for mask in range(8)]
    def choices(left,right):
        out=[]
        for chain in chains:
            walk=(() if left is None else (left,))+chain+(() if right is None else (right,))
            if all(rows[u]>>v&1 for u,v in zip(walk,walk[1:])):
                out.append(sum(base**ports.index(v) for v in chain))
        return out
    starts=[choices(None,v) for v in fixed];ends=[choices(v,None) for v in fixed]
    bridges=[[choices(u,v) if u!=v else [] for v in fixed] for u in fixed]
    dp=[{} for _ in range(1<<n)]
    for v in range(n):dp[1<<v][v]=dict(Counter(starts[v]))
    for mask in range(1,1<<n):
        for u,poly in dp[mask].items():
            for v in range(n):
                if mask>>v&1 or not bridges[u][v]:continue
                target=dp[mask|1<<v].setdefault(v,defaultdict(int))
                for c,number in poly.items():
                    for gap in bridges[u][v]:target[c+gap]+=number
        if mask!=(1<<n)-1:dp[mask].clear()
    out=defaultdict(int)
    for v,poly in dp[-1].items():
        for c,number in poly.items():
            for gap in ends[v]:
                x=c+gap;r=[]
                for _ in range(3):r.append(x%base);x//=base
                need(x==0,'no run-code overflow')
                if all(r):out[tuple(r)]+=number
    return dict(out)

def tail(rows,O):
    n=len(rows);out=[r|(1<<n)|((1<<(n+1)) if not O>>i&1 else 0) for i,r in enumerate(rows)]
    return out+[1<<(n+1),O]

def evaluate_runs(K,lengths):
    # Inclusion-exclusion counts surjections, independent of producer expo.
    from itertools import product
    def F(r,a):return sum((-1)**(r-j)*comb(r,j)*j**a for j in range(r+1))
    return sum(c*__import__('math').prod(F(r,a) for r,a in zip(rs,lengths)) for rs,c in K.items())

def run():
    data=load('synchronized_descent.json.gz',HERE);stats=Counter();receipts=[];start=time.time()
    for r in data['records']:
        ports=r['roots']+[r['clock']];D={tuple(x[:3]):x[3] for x in r['defect_runs']};U={tuple(x[:3]):x[3] for x in r['slope_runs']}
        dd=signed_kernel(r['rows'],r['selected'],ports)
        need(dd==D,'independent signed-cut full tensor')
        H=gaps_hamilton(tail(r['rows'],r['selected']),ports)
        expected={k:D.get(k,0)+2*U.get(k,0) for k in D.keys()|U.keys()};expected={k:v for k,v in expected.items() if v}
        need(H==expected,'independent entire Hamilton-tail gap polynomial')
        stats['full_tensor_identities']+=1
        for lengths in [(1,1,1),(2,2,1),(1,1,2),(2,1,2)]:
            rr,OO,_=expand(r['rows'],r['selected'],dict(zip(ports,lengths)))
            need(len(rr)<=18,'literal int64 bound')
            h,u,w,_=tables(np.array(rr,dtype=np.int64));sel=[v for v in range(len(rr)) if OO>>v&1]
            literalU=sum(int(u[v]) for v in sel);literalD=defect(h,w,sel)
            need(literalU==evaluate_runs(U,lengths) and literalD==evaluate_runs(D,lengths),'literal full endpoint readouts')
            stats['literal_controls']+=1
        vals=[]
        for m,k in [(1,1),(1,2),(2,1),(2,2),(3,1)]:
            z=2*m+r['z_offset']+k;raw=evaluate_runs(U,[m,m,z]);j=2*m+r['A_valuation_offset'];c=r['anchor']+2**k*sum(b**m*d for b,d in r['B_odd'])
            need(raw==2**j*c and c%2==1 and c>0,'exact normalized seed coefficient')
            vals.append({'m':m,'k':k,'z':z,'j':j,'c':c})
        receipts.append({'family':r['family'],'clock':r['clock'],'values':vals})
        print('audited',r['family'],r['clock'],'seconds',round(time.time()-start,1),flush=True)
    out={'scope':'All nine monomial-base4 tail charts; full three-port D by the independent one-cut walk DP, full U by independent gap DP on the two-vertex tail; four literal controls per chart include an off-diagonal hostile.','statistics':dict(stats),'records':receipts}
    save('independent_descent.json.gz',out);print(out['statistics'],flush=True)
if __name__=='__main__':run()
