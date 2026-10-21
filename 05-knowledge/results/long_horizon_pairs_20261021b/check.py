#!/usr/bin/env python3
"""Exact four-node cancellation flow with its indispensable integer defect."""
from collections import defaultdict, deque
from dataclasses import dataclass, replace
from hashlib import sha256
from importlib.util import spec_from_file_location, module_from_spec
from itertools import combinations, permutations, product
from math import gcd
from pathlib import Path
import argparse
import json

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
POINTS=((0,0),(0,1),(-1,-1),(1,0))
EDGES=tuple((a,b) for a in range(4) for b in range(4) if a!=b)
V=tuple((POINTS[b][0]-POINTS[a][0],POINTS[b][1]-POINTS[a][1]) for a,b in EDGES)
Q=lambda v:v[0]*v[0]-v[0]*v[1]+v[1]*v[1]
COST=tuple(map(Q,V))
LAMBDA=(-3,1,1,1)
CHECKS=0

def ck(ok,msg):
    global CHECKS
    CHECKS+=1
    if not ok:raise ValueError(msg)

def integer(x,lo=0):
    if type(x) is not int or x<lo:raise ValueError('literal integer')
    return x

def counts(c):
    if type(c) is not tuple or len(c)!=12:raise ValueError('immutable twelve-type inventory')
    for x in c:integer(x)
    return c

def divergence(x):
    b=[0]*4
    for n,(a,z) in zip(x,EDGES):b[a]-=n;b[z]+=n
    return tuple(b)

def vector(x):
    return tuple(sum(n*v[k] for n,v in zip(x,V)) for k in range(2))

def weight(x):return sum(a*b for a,b in zip(x,COST))

def bounds(c):
    counts(c);lo=hi=0;first=True
    for size in range(1,4):
        for S in combinations((1,2,3),size):
            incoming=sum(n for n,(a,b) in zip(c,EDGES) if a not in S and b in S)
            outgoing=sum(n for n,(a,b) in zip(c,EDGES) if a in S and b not in S)
            l=-(outgoing//size);h=incoming//size
            if first:lo,hi=l,h;first=False
            else:lo=max(lo,l);hi=min(hi,h)
    ck(lo<=0<=hi,'empty cancellation feasible')
    return lo,hi

def residual(c,x):
    arcs=[]
    for i,((a,b),cap,flow,w) in enumerate(zip(EDGES,c,x,COST)):
        if flow<cap:arcs.append((a,b,w,cap-flow,i,1))
        if flow:arcs.append((b,a,-w,flow,i,-1))
    return arcs

CYCLES=tuple(p for k in range(2,5) for ss in combinations(range(4),k)
             for p in permutations(ss) if p[0]==min(p))

def feasible(c,t):
    """Integral source/sink maxflow, retaining original edge identities."""
    b=tuple(t*a for a in LAMBDA);x=[0]*12;remaining=list(b)
    while any(n<0 for n in remaining):
        arcs=residual(c,x);seen={a:None for a in range(4) if remaining[a]<0}
        queue=deque(seen);end=None
        while queue and end is None:
            a=queue.popleft()
            for e in arcs:
                if e[0]!=a or e[1] in seen:continue
                z=e[1];seen[z]=e
                if remaining[z]>0:end=z;break
                queue.append(z)
        if end is None:raise ValueError('infeasible prescribed divergence')
        path=[];a=end
        while seen[a] is not None:path.append(seen[a]);a=seen[a][0]
        amount=min(-remaining[a],remaining[end],*(e[3] for e in path))
        for e in path:x[e[4]]+=e[5]*amount
        remaining[a]+=amount;remaining[end]-=amount
    ck(not any(remaining),'all prescribed flow demands met')
    return x

@dataclass(frozen=True)
class Fixed:
    inventory:tuple
    defect:int
    selected:tuple
    value:int
    potential:tuple

def fixed(c,t):
    counts(c)
    if type(t) is not int:raise ValueError('literal signed defect')
    lo,hi=bounds(c)
    if not lo<=t<=hi:raise ValueError('infeasible defect')
    x=feasible(c,t)
    while True:
        best={}
        for e in residual(c,x):
            if e[:2] not in best or e[2]>best[e[:2]][2]:best[e[:2]]=e
        chosen=None;gain=0
        for p in CYCLES:
            es=[best.get((p[k],p[(k+1)%len(p)])) for k in range(len(p))]
            if any(e is None for e in es):continue
            g=sum(e[2] for e in es)
            if g>gain:chosen=es;gain=g
        if chosen is None:break
        amount=min(e[3] for e in chosen)
        for e in chosen:x[e[4]]+=amount*e[5]
    pi=[0]*4
    for _ in range(3):
        old=pi.copy()
        for a,b,w,_,_,_ in residual(c,x):pi[b]=max(pi[b],old[a]+w)
    out=Fixed(c,t,tuple(x),weight(x),tuple(pi));audit_fixed(out)
    return out

def audit_fixed(r):
    if type(r) is not Fixed:raise ValueError('literal flow certificate')
    counts(r.inventory);counts(r.selected);integer(r.value)
    if type(r.defect) is not int:raise ValueError('literal signed defect')
    if type(r.potential) is not tuple or len(r.potential)!=4 or any(type(v) is not int for v in r.potential):
        raise ValueError('literal integral potential')
    ck(all(x<=c for x,c in zip(r.selected,r.inventory)),'capacity/source ownership')
    ck(divergence(r.selected)==tuple(r.defect*a for a in LAMBDA),'exact divergence sidecar')
    ck(vector(r.selected)==(0,0),'literal zero charge')
    ck(r.value==weight(r.selected),'owned diagonal norm')
    for a,b,w,_,_,_ in residual(r.inventory,r.selected):
        ck(w<=r.potential[b]-r.potential[a],'residual potential optimality')
    return r

@dataclass(frozen=True)
class Optimum:
    best:Fixed
    neighbours:tuple

def optimize(c):
    lo,hi=bounds(c);cache={}
    def f(t):
        if t not in cache:cache[t]=fixed(c,t)
        return cache[t]
    while lo<hi:
        mid=(lo+hi)//2
        if f(mid).value<f(mid+1).value:lo=mid+1
        else:hi=mid
    lower,upper=bounds(c)
    ans=Optimum(f(lo),tuple(f(t) for t in (lo-1,lo+1) if lower<=t<=upper))
    audit_optimum(ans);return ans

def audit_optimum(r):
    if type(r) is not Optimum or type(r.neighbours) is not tuple:raise ValueError('literal optimum certificate')
    audit_fixed(r.best);c=r.best.inventory;t=r.best.defect;lo,hi=bounds(c)
    ck(lo<=t<=hi,'feasible maximizing defect')
    expected=tuple(j for j in (t-1,t+1) if lo<=j<=hi)
    ck(tuple(x.defect for x in r.neighbours)==expected,'complete adjacent defect witnesses')
    for x in r.neighbours:
        audit_fixed(x);ck(x.inventory==c and x.value<=r.best.value,'concave local maximum')
    return r

def brute(c):
    return max(weight(x) for x in product(*(range(n+1) for n in c)) if vector(x)==(0,0))

def short_circuits():
    out=[]
    for i,j in combinations(range(12),2):
        if tuple(a+b for a,b in zip(V[i],V[j]))==(0,0):
            x=[0]*12;x[i]=x[j]=1;out.append(tuple(x))
    det=lambda a,b:a[0]*b[1]-a[1]*b[0]
    for ids in combinations(range(12),3):
        a,b,c=(V[i] for i in ids);v=(det(b,c),det(c,a),det(a,b))
        if all(n<0 for n in v):v=tuple(-n for n in v)
        if not all(n>0 for n in v):continue
        g=gcd(gcd(v[0],v[1]),v[2]);x=[0]*12
        for i,n in zip(ids,v):x[i]=n//g
        out.append(tuple(x))
    return tuple(sorted(out))

def greedy_short(c):
    counts(c);left=list(c)
    for v in short_circuits():
        m=min(left[i]//n for i,n in enumerate(v) if n)
        left=[a-m*b for a,b in zip(left,v)]
    return tuple(left)

def encode(obj):return (json.dumps(obj,sort_keys=True,separators=(',',':'))+'\n').encode()

def payload(r):
    def one(x):return {'defect':x.defect,'selected':x.selected,'value':x.value,'potential':x.potential}
    return {'counts':r.best.inventory,'best':one(r.best),'neighbours':[one(x) for x in r.neighbours]}

def core_controls():
    ck(len(set(V))==12 and sorted(COST)==[1]*6+[3]*6,'exact signed alphabet')
    tested=0
    for mask in range(1<<12):
        c=tuple((mask>>i)&1 for i in range(12));o=optimize(c)
        ck(o.best.value==brute(c),'complete binary-inventory independent enumeration');tested+=1
    multiple=0
    for ids in combinations(range(12),3):
        for ns in product(range(1,4),repeat=3):
            c=[0]*12
            for i,n in zip(ids,ns):c[i]=n
            c=tuple(c);ck(optimize(c).best.value==brute(c),'three-type multiplicity universe');multiple+=1
    # Three outward spokes have zero charge but are not graph circulation.
    star=tuple(int(a==0) for a,b in EDGES)
    ck(vector(star)==(0,0) and divergence(star)==(-3,1,1,1),'indispensable defect hostile')
    ck(fixed(star,0).value==0 and optimize(star).best.value==3,'ordinary circulation misses all star norm')
    # No irrevocably increasing optimal packet bank exists for the full alphabet.
    early=tuple(int(i in (0,1,3)) for i in range(12))
    late=tuple(int(i in (0,1,2,3)) for i in range(12))
    eo,lo=optimize(early),optimize(late)
    es=[x for x in product(*(range(n+1) for n in early)) if vector(x)==(0,0) and weight(x)==eo.best.value]
    ls=[x for x in product(*(range(n+1) for n in late)) if vector(x)==(0,0) and weight(x)==lo.best.value]
    ck(eo.best.value==2 and lo.best.value==3 and len(es)==len(ls)==1,'unique temporal optima')
    ck(not all(a<=b for a,b in zip(es[0],ls[0])),'optimal cancellation requires ownership reassignment')
    # A zero-free integer inventory need not occupy one strict angular half-plane.
    holes=(1,1,0,0,0,1,0,1,0,0,0,0)
    ck(optimize(holes).best.value==0 and weight(holes)==8 and Q(vector(holes))==1,'integer residue hostile')
    quartet=tuple(int(674>>i&1) for i in range(12))
    ck(vector(quartet)==(0,0) and weight(quartet)==8,'irreducible four-type zero packet')
    ck(greedy_short(quartet)==quartet and optimize(quartet).best.value==8,'strict short-circuit-bank extension')
    ck(all(vector(x)!=(0,0) for x in product(*(range(n+1) for n in quartet)) if any(x) and x!=quartet),
       'no proper nonempty zero subpacket')
    ck(len(short_circuits())==46,'six pairs and forty primitive three-type circuits')
    def rejects(f):
        try:f()
        except (TypeError,ValueError):return True
        return False
    for bad in ([0]*12,(False,)+(0,)*11,(0.0,)+(0,)*11,(-1,)+(0,)*11):
        ck(rejects(lambda:optimize(bad)),'literal inventory rejection')
    r=optimize(star)
    for bad in (replace(r.best,defect=True),replace(r.best,value=2),replace(r.best,potential=(-99,0,0,0))):
        ck(rejects(lambda:audit_fixed(bad)),'forged fixed-flow receipt rejection')
    ck(rejects(lambda:audit_optimum(replace(r,neighbours=()))),'missing adjacent optimum rejection')
    return {'binary_inventories':tested,'three_type_multiplicity_inventories':multiple,'star':payload(r),
            'non_nested':{'early':payload(eo),'late':payload(lo)},
            'zero_free_integer_hostile':{'counts':holes,'norm':8,'vector':vector(holes),'energy':1},
            'irreducible_quartet':payload(optimize(quartet))}

FRONTIER_SHA='391f6d3380bcae4afb9a1bbd9522b74b1e9006012fa10a92be5bf0bbf938754f'

def source_rows(N,terminals):
    integer(N,1)
    if type(terminals) is not tuple or len(set(terminals))!=len(terminals):raise ValueError('immutable unique terminals')
    for u in terminals:integer(u,1)
    path=HERE.parent/'session18_native_annulus/verify_annulus.py'
    sp=spec_from_file_location('long_pair_native',path);native=module_from_spec(sp);sp.loader.exec_module(native)
    raw=(path.parent/'frontier.json').read_bytes().replace(b'\r\n',b'\n')
    ck(sha256(raw).hexdigest()==FRONTIER_SHA,'frozen grounded frontier')
    data=json.loads(raw);routes=native.Routes(data['frontier'],dict(data['frontier_ranks']))
    ck(set(terminals)<=routes.V,'actual terminal universe')
    rows=[];digest=sha256();total_norm=0;events=defaultdict(lambda:[0,0,0,0,0])
    for s in range(1,256*N,2):
        if s%8==5:continue
        u,t,d=routes.get(s);v=native.charge(64*N,s)
        ck(v==native.independent_charge(64*N,s),'independent original source charge')
        row=(s,s%15,s%35,u,t,d,*v,Q(v));digest.update(encode(row))
        signed_v=tuple((1-2*(t%2))*x for x in v);event=events[d,u]
        k=2 if s>=4*N else 0
        event[k]+=signed_v[0];event[k+1]+=signed_v[1]
        if k==2:event[4]+=Q(v)
        if s>=4*N:total_norm+=Q(v)
        else:ck(v==native.charge(N,s),'unchanged old-source charge')
        if u in terminals:rows.append(row)
    ck(total_norm==153*N,'exact complete new-source norm')
    W=defaultdict(lambda:(0,0));H=defaultdict(lambda:(0,0));A=delta=0;ledger=[]
    for d in range(257):
        for u in routes.V:
            wa,wb,ha,hb,m=events.get((d,u),(0,0,0,0,0))
            if not (wa or wb or ha or hb or m):continue
            w,h=W[u],H[u];before=Q((w[0]+h[0],w[1]+h[1]))-Q(w)
            w=(w[0]+wa,w[1]+wb);h=(h[0]+ha,h[1]+hb);W[u],H[u]=w,h
            delta+=Q((w[0]+h[0],w[1]+h[1]))-Q(w)-before;A+=m
        ledger.append([d,A,153*N-A,delta-A,delta,delta-153*N])
    if N==7:ck(ledger[8][1:]==[1054,17,30,1084,13],'unchanged native early-arrival hostile')
    return tuple(rows),digest.hexdigest(),routes.audit(),ledger

def typed_row(r):
    if type(r) is not tuple or len(r)!=9 or any(type(x) is not int for x in r):raise ValueError('literal source row')
    s,r15,r35,u,t,d,a,b,q=r
    if s<=0 or s%2!=1 or s%8==5 or u<=0 or u%2!=1 or u%8==5 or min(t,d)<0:
        raise ValueError('native row domain')
    if (r15,r35)!=(s%15,s%35) or q!=Q((a,b)) or ((a,b)!=(0,0) and (a,b) not in V):
        raise ValueError('source register/charge identity')
    return r

def signed(r):return tuple((1-2*(r[4]%2))*x for x in r[6:8])

def owned_selection(N,rows,d,o):
    """All rows must separately be source-authenticated; this checks ownership."""
    integer(N,1);integer(d);audit_optimum(o)
    if type(rows) is not tuple:raise ValueError('immutable source rows')
    seen=set();buckets=[[] for _ in V];terminals=set()
    for r in rows:
        typed_row(r)
        if r[0] in seen:raise ValueError('duplicate physical source')
        seen.add(r[0]);terminals.add(r[3])
        if r[0]<4*N or r[5]>d:continue
        v=signed(r)
        if v!=(0,0):buckets[V.index(v)].append(r)
    ck(len(terminals)==1,'no cross-terminal vector cancellation')
    ck(tuple(map(len,buckets))==o.best.inventory,'complete arrived new-source inventory')
    chosen=tuple(sorted(r[0] for k,b in zip(o.best.selected,buckets) for r in sorted(b,key=lambda x:(x[5],x[0]))[:k]))
    ck(len(chosen)==len(set(chosen)),'single source ownership')
    table={r[0]:r for r in rows};vec=tuple(sum(signed(table[s])[i] for s in chosen) for i in range(2))
    ck(vec==(0,0) and sum(table[s][8] for s in chosen)==o.best.value,'whole mature zero packet')
    ck(all(table[s][5]<=d and table[s][0]>=4*N for s in chosen),'no future credit or old-source norm')
    return chosen

def native_controls(N,terminals):
    rows,source_hash,dag,ledger=source_rows(N,terminals);results=[]
    for u in terminals:
        local=tuple(r for r in rows if r[3]==u);batches=defaultdict(list)
        for r in local:batches[r[5]].append(r)
        c=[0]*12;W=[0,0];A=0;profiles=[];best=None;certs={};last=None
        for d in range(257):
            for r in batches[d]:
                v=signed(r)
                if r[0]<4*N:W=[a+b for a,b in zip(W,v)]
                else:
                    A+=r[8]
                    if v!=(0,0):c[V.index(v)]+=1
            cc=tuple(c)
            if cc!=last:o=optimize(cc);last=cc
            z=o.best.value;opposite=sum(min(c[i],c[V.index(tuple(-a for a in v))])*COST[i] for i,v in enumerate(V))
            H=vector(c);energy=Q((W[0]+H[0],W[1]+H[1]))-Q(W)
            # energy is Q(W+H)-Q(W); the full pair ledger is energy-A.
            rem=tuple(n-x for n,x in zip(c,o.best.selected))
            ck(vector(rem)==H and weight(rem)==A-z,'source-preserving residual')
            coherence=energy-weight(rem);interaction=energy-A
            ck(coherence-z==interaction,'physical pair ledger unchanged')
            ck(z>=opposite,'all-charge optimum contains antipodal feasible choice')
            profiles.append([d,A,z,opposite,z-opposite,weight(rem),*H,*W,interaction,coherence])
            if best is None or z-opposite>best[0]:best=(z-opposite,d,o)
        gain,d,o=best;selected=owned_selection(N,local,d,o)
        greedy_rem=greedy_short(o.best.inventory)
        ck(vector(greedy_rem)==vector(o.best.inventory),'short-circuit comparison source vector')
        ck(weight(o.best.inventory)-weight(greedy_rem)<=o.best.value,'exact optimum dominates greedy short circuits')
        results.append({'terminal':u,'profile':profiles,'best_gain':gain,'best_depth':d,
                        'optimum':payload(o),'selected_source_ids':selected,'source_rows':local,
                        'greedy_short_comparison':{'residual':greedy_rem,'cancelled_norm':weight(o.best.inventory)-weight(greedy_rem)}})
    return {'N':N,'terminals':terminals,'source_rows_sha256':source_hash,'route_dag':dag,'results':results,
            'full_ledger_schema':['depth','arrived_norm','moving_norm','I','Delta_E','excess_over_153N'],
            'full_ledger':ledger}

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--write',action='store_true');args=parser.parse_args()
    out={'status':'PROVED exact finite-inventory cancellation; native long-horizon uniform energy OPEN',
         'core':core_controls(),'native':[native_controls(7,(1,3,7)),native_controls(4096,(1,87,711))],
         'profile_schema':['depth','arrived_norm','optimal_cancelled_norm','antipodal_cancelled_norm','gain',
                           'residual_norm','new_vector_a','new_vector_b','old_vector_a','old_vector_b','I_terminal','C_terminal'],
         'checks':CHECKS}
    data=encode(out);path=HERE/'receipt.json'
    if args.write:path.write_bytes(data)
    else:ck(path.read_bytes().replace(b'\r\n',b'\n')==data,'saved receipt equality')
    print(json.dumps({'checks':CHECKS,'sha256_lf':sha256(data).hexdigest()}))

if __name__=='__main__':main()
