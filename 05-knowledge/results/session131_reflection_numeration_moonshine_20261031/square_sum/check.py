#!/usr/bin/env python3
"""Portable exact verifier. No optimizer, external inputs, or assert statements."""
import argparse
from collections import Counter
from itertools import product
from math import isqrt
from pathlib import Path
import hashlib
import json

HERE = Path(__file__).resolve().parent
CHECKS = 0
BOUNDS = (4, 7, 8, 7, 8, 7, 8, 7)
BLOCK = (0, 0, 1, 2, 3, 3, 2, 1)

def require(test, message):
    global CHECKS
    CHECKS += 1
    if not test:
        raise ValueError(message)

def is_square(n):
    return n >= 0 and isqrt(n) ** 2 == n

def cycle_edges(cycle):
    return list(zip(cycle, cycle[1:] + cycle[:1]))

def valid_cycle(cycle, n):
    return (n >= 3 and len(cycle) == n and set(cycle) == set(range(1, n + 1))
            and all(a != b and is_square(a+b) for a, b in cycle_edges(cycle)))

def odd_forest(n, edges):
    adj = {v: [] for v in range(1, n+1)}
    odd, even = [], []
    for a, b in edges:
        require(a != b and is_square(a+b), 'invalid supplied square edge')
        if isqrt(a+b) % 2:
            odd.append((a,b)); adj[a].append(b); adj[b].append(a)
        else:
            even.append((a,b))
    remaining = set(adj)
    components = []
    k = [0]*4
    endpoints = Counter()
    for v in adj:
        require(len(adj[v]) <= 2, 'odd subgraph degree exceeds two')
        endpoints[v%8] += 2-len(adj[v])
    while remaining:
        start = min(remaining); remaining.remove(start); component = [start]
        for v in component:
            for w in adj[v]:
                if w in remaining:
                    remaining.remove(w); component.append(w)
        require(sum(len(adj[v]) for v in component) == 2*(len(component)-1),
                'cut subgraph contains a cycle, not a path')
        types = {BLOCK[v%8] for v in component}
        require(len(types) == 1, 'odd path crosses a residue block')
        k[next(iter(types))] += 1
        require(sum(2-len(adj[v]) for v in component) == 2,
                'lost singleton endpoint incidence')
        components.append(component)
    nr = Counter(v%8 for v in adj)
    delta = [nr[1]-nr[0],nr[7]-nr[2],nr[3]-nr[6],nr[5]-nr[4]]
    for j,(odd_r,even_r) in enumerate(((1,0),(7,2),(3,6),(5,4))):
        require(endpoints[odd_r] == k[j]+delta[j], 'odd endpoint formula')
        require(endpoints[even_r] == k[j]-delta[j], 'even endpoint formula')
    return odd,even,components,k,endpoints,delta

def audit_cycle(cycle, n):
    require(valid_cycle(cycle,n), 'not a Hamiltonian cycle')
    odd,even,components,k,E,delta = odd_forest(n,cycle_edges(cycle))
    t=len(even); X=k[0]+k[3];Y=k[1]+k[2]
    dl=delta[0]+delta[3];dr=delta[1]+delta[2]
    require(t == len(components) == sum(k), 'cut forest component count')
    require(all(x>=1 for x in k), 'empty residue block')
    require(X-Y == dr-dl and (X-dl)%2 == (Y-dr)%2 == 0, 'flow identities')
    require(t >= BOUNDS[n%8], 'uniform bound violated')
    oo=sum(a%2==b%2==1 for a,b in even)
    ee=sum(a%2==b%2==0 for a,b in even)
    require(oo-ee == n%2, 'label parity edge identity')
    # Each removed edge supplies two named halfedge slots. Each path pairs
    # its two slots, even when they belong to the same singleton vertex.
    slots={v:[] for v in range(1,n+1)};seam={};pair={}
    for a,b in even:
        slots[a].append((a,b));slots[b].append((b,a))
        seam[(a,b)]=(b,a);seam[(b,a)]=(a,b)
    for component in components:
        boundary=[s for v in component for s in slots[v]]
        require(len(boundary)==2,'wrong named boundary size')
        a,b=boundary;pair[a]=b;pair[b]=a
    require(set(pair)==set(seam) and len(pair)==2*t,'incidence slot domain')
    unseen=set(pair);boundary_components=0
    while unseen:
        first=unseen.pop();queue=[first];boundary_components+=1
        for s in queue:
            for z in (pair[s],seam[s]):
                if z in unseen:unseen.remove(z);queue.append(z)
    require(boundary_components==1,'sewn boundary is disconnected')
    return {'n':n,'even_square_edges':t,'bound':BOUNDS[n%8],
            'odd_root_path_counts':k,'endpoint_counts_mod8':[E[r] for r in range(8)],
            'even_root_seams':[list(x) for x in even],
            'boundary_components':boundary_components,
            'singleton_paths':sum(len(c)==1 for c in components)}

def audit_path(path):
    n=len(path);edges=list(zip(path,path[1:]))
    odd,even,components,k,E,delta=odd_forest(n,edges)
    b=Counter((path[0]%8,path[-1]%8));X=k[0]+k[3];Y=k[1]+k[2]
    dl=delta[0]+delta[3];dr=delta[1]+delta[2]
    require(len(even)+1==len(components)==X+Y,'marked path count')
    require(X-Y==dr-dl+b[1]+b[5]-b[3]-b[7],'marked odd flow')
    require((X-dl-b[0]-b[4])%2==0,'marked even flow left')
    require((Y-dr-b[2]-b[6])%2==0,'marked even flow right')
    if BLOCK[path[0]%8]==BLOCK[path[-1]%8]:
        require(k[BLOCK[path[0]%8]]>=2,'same-block endpoint path count')

def compiled_starts(roots,n):
    m=len(roots)
    if m<3:return []
    cs=[0]
    for s in roots:cs.append(s*s-cs[-1])
    if m%2==0 and cs[-1]!=0:return []
    if m%2 and cs[-1]%2:return []
    lo=1;hi=n
    for i,c in enumerate(cs[:-1]):
        a,b=(1-c,n-c) if i%2==0 else (c-n,c-1)
        lo=max(lo,a);hi=min(hi,b)
    forbidden=set()
    for i in range(m):
        for j in range(i):
            if i%2==j%2:
                if cs[i]==cs[j]:return []
            else:
                even_c,odd_c=(cs[i],cs[j]) if i%2==0 else (cs[j],cs[i])
                if (odd_c-even_c)%2==0:forbidden.add((odd_c-even_c)//2)
    if m%2:
        x=cs[-1]//2
        return [x] if lo<=x<=hi and x not in forbidden else []
    return [x for x in range(lo,hi+1) if x not in forbidden]

def direct_starts(roots,n):
    result=[]
    for x in range(1,n+1):
        walk=[x]
        for s in roots:walk.append(s*s-walk[-1])
        if (walk[-1]==x and len(set(walk[:-1]))==len(roots)
                and all(1<=v<=n for v in walk[:-1])):
            result.append(x)
    return result

def audit_reflections():
    words=0
    for n in range(1,17):
        alphabet=range(2,isqrt(2*n-1)+1)
        for length in range(3,7):
            for roots in product(alphabet,repeat=length):
                require(compiled_starts(roots,n)==direct_starts(roots,n),
                        'reflection compiler disagrees with direct walking')
                words+=1
    require(compiled_starts((3,3,3,3),20)==[], 'closure-only hostile accepted')
    require(is_square(1+3) and not is_square((25-1)+(25-3)),
            'partial reflection automorphism hostile')
    for r in range(1,15):
        for s in range(1,30):
            for x in range(1,s*s):
                require((r*s)**2-r*r*x==r*r*(s*s-x),'scaling conjugacy')
    return {'root_words':words,'max_graph_order':16,'lengths':[3,6]}

def audit_components():
    edges40={(a,b) for a in range(1,41) for b in range(a+1,41)
             if is_square(a+b) and isqrt(a+b)%2}
    supplied=set();vertices=[]
    for a in range(1,5):
        for p in ([24+a,25-a,a,9-a,16+a,33-a],
                  [32+a,17-a,8+a,41-a]):
            vertices.extend(p)
            supplied.update(tuple(sorted(e)) for e in zip(p,p[1:]))
    require(set(vertices)==set(range(1,41)) and len(vertices)==40,'H40 partition')
    require(supplied==edges40,'H40 complete edge template')
    parent=list(range(2001));count=0;changes=[];previous=None
    def root(v):
        while parent[v]!=v:
            parent[v]=parent[parent[v]];v=parent[v]
        return v
    for n in range(1,2001):
        count+=1;neighbors=[]
        for s in range(1,isqrt(2*n-1)+1,2):
            v=s*s-n
            if 1<=v<n:
                neighbors.append(v)
                require(BLOCK[v%8]==BLOCK[n%8],'odd-edge channel')
                a,b=root(v),root(n)
                if a!=b:parent[a]=b;count-=1
        if n>=13:require(bool(neighbors),'odd-square arrival gap')
        if n>=44:require(count==4,'four-component theorem')
        if count!=previous:changes.append([n,count]);previous=count
    return {'maximum_order':2000,'component_count_changes':changes}

def fibonacci_weights(length):
    fs=[0,1]
    for i in range(length+1):fs.append(fs[-1]+fs[-2])
    return fs[2:2+length]

def zeck(n):
    if n==0:return []
    ws=[1,2]
    while ws[-1]<=n:ws.append(ws[-1]+ws[-2])
    if ws[-1]>n:ws.pop()
    out=[0]*len(ws)
    for i in range(len(ws)-1,-1,-1):
        if ws[i]<=n:out[i]=1;n-=ws[i]
    require(n==0,'greedy remainder')
    return out

def audit_digits():
    words=legal_count=0
    for length in range(19):
        ws=fibonacci_weights(length)
        for z in product((0,1),repeat=length):
            words+=1
            legal=all(not(z[i] and z[i+1]) for i in range(length-1))
            e=z[::2];o=z[1::2]
            row_legal=all(not(e[j] and o[j]) and
                          not(o[j] and (e[j+1] if j+1<len(e) else 0))
                          for j in range(len(o)))
            require(legal==row_legal,'interweaving legality equivalence')
            if not legal:continue
            legal_count+=1
            v=sum(x*w for x,w in zip(z,ws))
            maxj=max(len(e),len(o)+1)
            As=fibonacci_weights(2*maxj)[::2]
            signed=[]
            for j in range(maxj):
                ej=e[j] if j<len(e) else 0
                oj=o[j] if j<len(o) else 0
                prev=o[j-1] if 0<=j-1<len(o) else 0
                signed.append(ej+prev-oj)
            require(sum(c*a for c,a in zip(signed,As))==v,'signed even-weight map')
            require(v%2==sum(x for i,x in enumerate(z) if (i+2)%3)%2,
                    'Fibonacci weight parity')
            normalized=list(z)
            while normalized and normalized[-1]==0:normalized.pop()
            require(zeck(v)==normalized,'canonical digit roundtrip')
    reflected=0
    for n in range(1,201):
        for s in range(2,isqrt(2*n-1)+1):
            for x in range(1,n+1):
                y=s*s-x
                if not(1<=y<=n and x!=y):continue
                zx,zy=zeck(x),zeck(y)
                require(sum(a*b for a,b in zip(zy,fibonacci_weights(len(zy))))==y,
                        'conjugated reflection value')
                require(zeck(s*s-y)==zx,'conjugated reflection involution')
                require((x+y)%2==s%2,'conjugated reflection parity')
                reflected+=1
    return {'all_binary_words':words,'legal_words':legal_count,
            'maximum_digit_length':18,'conjugated_reflections':reflected}

def audit_flow_floor():
    computed=[]
    for r in range(8):
        nr=Counter(v%8 for v in range(1,80+r+1))
        ds=(nr[1]-nr[0],nr[7]-nr[2],nr[3]-nr[6],nr[5]-nr[4])
        dl=ds[0]+ds[3];dr=ds[1]+ds[2]
        feasible=[]
        for ks in product(range(1,9),repeat=4):
            x=ks[0]+ks[3];y=ks[1]+ks[2]
            if x-y==dr-dl and (x-dl)%2==0 and (y-dr)%2==0:
                feasible.append(sum(ks))
        raw=min(feasible)
        require(raw==(4,7,8,7,4,7,8,7)[r],'flow algebra floor')
        if r==4:
            # t=4 forces k_i=1; residue2 is the only allowed endpoint pool
            # in the {2,6} even-seam sector, hence its path closes alone.
            E2=1-ds[1];E6=1-ds[2]
            require(E2==2 and E6==0,'residue4 disconnected closure obstruction')
            raw=min(t for t in feasible if t>4)
        require(raw==BOUNDS[r],'connected flow floor')
        computed.append(raw)
    return computed

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path)
    args=parser.parse_args()
    data=json.loads((HERE/'witnesses.json').read_text())
    rows=[];path_cases=0
    require([w['n'] for w in data['cycles']]==list(range(80,88)),'residue witnesses incomplete')
    for w in data['cycles']:
        row=audit_cycle(w['cycle'],w['n'])
        require(row['even_square_edges']==w['even_square_edges']==BOUNDS[w['n']%8],
                'witness does not attain elementary bound')
        rows.append(row)
        roots=tuple(isqrt(a+b) for a,b in cycle_edges(w['cycle']))
        require(w['cycle'][0] in compiled_starts(roots,w['n']),
                'long Hamiltonian reflection word compiler')
        for i in range(w['n']):
            path=w['cycle'][i+1:]+w['cycle'][:i+1]
            audit_path(path);path_cases+=1
    c32=[1,8,28,21,4,32,17,19,30,6,3,13,12,24,25,11,5,31,18,7,29,20,16,9,27,22,14,2,23,26,10,15]
    inherited=audit_cycle(c32,32)
    require(inherited['even_square_edges']==16,'G32 inherited edge statistic')
    require(inherited['singleton_paths']>0,'singleton-incidence positive control')
    bad=c32.copy();bad[-1]=bad[0]
    require(not valid_cycle(bad,32),'duplicate-label hostile accepted')
    bad=c32.copy();bad[0],bad[1]=bad[1],bad[0]
    require(not valid_cycle(bad,32),'nonsquare-edge hostile accepted')
    floors=audit_flow_floor()
    reflection=audit_reflections();components=audit_components();digits=audit_digits()
    result={'status':'PASS','scope':'Universal elementary bounds plus finite exact verification; no all-order attainment claim.',
            'uniform_residue_floor':floors,'optimal_witnesses':rows,'marked_path_cases':path_cases,
            'inherited_G32_even_root_edges':16,'reflection_compiler':reflection,
            'odd_square_components':components,'digit_interweaving':digits,
            'checked_conditions':CHECKS,
            'input_sha256':{name:hashlib.sha256((HERE/name).read_bytes()).hexdigest()
                            for name in ('check.py','witnesses.json','source_pins.json')}}
    out=json.dumps(result,sort_keys=True,indent=2)+'\n'
    if args.output:args.output.write_text(out)
    else:print(out,end='')
if __name__=='__main__':main()
