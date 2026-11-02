#!/usr/bin/env python3
"""Independent Session140 finite controls. Written proofs establish infinite claims."""
import collections
import hashlib
import itertools
import json
import pathlib

ROOT = pathlib.Path(__file__).resolve().parent.parent
CHECKS = 0

def require(condition, label):
    global CHECKS
    CHECKS += 1
    if not condition:
        raise RuntimeError(label)

def ordinary(n):
    return 3*n+1 if n & 1 else n//2

def shortcut(n):
    return (3*n+1)//2 if n & 1 else n//2

def iterate(fn,n,d):
    for _ in range(d):
        n=fn(n)
    return n

# Derive each safe gate by searching inverse paths, independently of packet table.
def gate(n):
    return next(k for k in range(1,5) if (2**k*n)%9 in (4,7))

def selected(n,fn):
    k=gate(n)
    numerator=16*n-2**(4-k)
    require(numerator%3==0, 'integrality')
    return ((32 if fn is ordinary else 16)*n,numerator//3)

residue_gates={r:gate(r) for r in (1,2,4,5,7,8)}
require(residue_gates=={1:2,2:1,4:2,5:3,7:4,8:1},'exhaustive gate table')

# Forward-only exhaustive enumeration of *complete* fibres: any d-step
# predecessor of m is <=2**d*m because both maps send n to >=n/2.
fibres={}
for fn,depth,name in ((ordinary,5,'C'),(shortcut,4,'T')):
    all_fibres=collections.defaultdict(set)
    for p in range(1,2**depth*200+1):
        all_fibres[iterate(fn,p,depth)].add(p)
    for n in range(1,201):
        if n%3:
            a,b=selected(n,fn)
            require({a,b}<=all_fibres[n],(name,'full forward fibre',n))
            require(a!=b and min(a,b)>n and a%3 and b%3,(name,'two increasing units',n))
    fibres[name]={'roots':200,'forward_enumeration_ceiling':2**depth*200}

for fn,depth,name,expected in ((ordinary,4,'C',{112,18}),(shortcut,3,'T',{56,18})):
    actual={p for p in range(1,2**depth*7+1) if iterate(fn,p,depth)==7}
    require(actual==expected,('sharpness complete fibre',name))
    require(len({p for p in actual if p%3})==1,('sharpness survives',name))
    for d in range(depth):
        earlier={p for p in range(1,2**d*7+1) if iterate(fn,p,d)==7}
        require(len({p for p in earlier if p%3})==1,('earlier clock',name,d))

large_inputs=set(range(1,100001))
for k in range(1,1025):
    large_inputs.update({2**k-1,2**k,2**k+1})
for n in large_inputs:
    exponent=n.bit_length()-1
    nxt=ordinary(n)
    delta=nxt.bit_length()-1-exponent
    require(delta in ((1,2) if n&1 else (-1,)),('dyadic exponent',n))
    require(exponent%3 != (nxt.bit_length()-1)%3,('proper explicit coloring',n))
    if n%3:
        for fn,depth in ((ordinary,5),(shortcut,4)):
            a,b=selected(n,fn)
            require(iterate(fn,a,depth)==n and iterate(fn,b,depth)==n,('large forward clock',n))
            require(min(a,b)>n and a!=b and a%3 and b%3,('large unit growth',n))
require({(n.bit_length()-1)%3 for n in (1,2,4)}=={0,1,2},'triangle needs three colors')

macro_receipts=[]
for fn,name in ((ordinary,'C'),(shortcut,'T')):
    for root in (1,5,7,11):
        frontier={root}; seen={root}
        for depth in range(1,13):
            children=[]
            for p in frontier:
                children.extend(selected(p,fn))
            nxt=set(children)
            require(len(nxt)==len(children)==2**depth,('same-level separation',name,root,depth))
            require(not nxt&seen,('cross-level separation',name,root,depth))
            seen.update(nxt);frontier=nxt
        macro_receipts.append({'map':name,'root':root,'blocks':12,'leaves':len(frontier),'all_levels':len(seen)})

# Root transfer uses only a finite arithmetic calculation, not convergence.
for x in range(1,10001):
    y=x;t=0;v2=0
    while y%2==0:
        y//=2;v2+=1
    if x%3==0:
        z=3*y+1;t=v2+1
    else:
        z=x
    require(z%3 and iterate(ordinary,x,t)==z and t<=v2+1,('root transfer',x))

# Read actual honeycomb edges, contract matching, then compare quotient neighbors.
def neighbors(p):
    x,y=p
    return {(x-1,y),(x+1,y),(x,y+(1 if (x+y)%2==0 else -1))}

def match_lower(p):
    x,y=p
    return (x,y if (x+y)%2==0 else y-1)

def square_label(p):
    x,y=match_lower(p)
    return ((x+y)//2,(x-y)//2)

for x in range(-8,9):
    for y in range(-8,9):
        if (x+y)%2:
            continue
        branch={(x,y),(x,y+1)}
        actual={square_label(v) for p in branch for v in neighbors(p) if v not in branch}
        a,b=square_label((x,y))
        expected={(a-1,b),(a+1,b),(a,b-1),(a,b+1)}
        require(actual==expected,('actual matching quotient',x,y))
hexagon=((0,0),(1,0),(2,0),(2,1),(1,1),(0,1))
require(all(hexagon[(i+1)%6] in neighbors(p) for i,p in enumerate(hexagon)),'six-cycle host edges')
require(all((sum(hexagon[i])%2)!=(sum(hexagon[(i+1)%6])%2) for i in range(6)),'subdivision is two-colorable')
require(all((sum(hexagon[i])%2)==0 for i in (0,2,4)),'three original cycle labels now share host color')

# Median lemma on every labeled tree up to order six and all repeated triples.
def tree_edges(code,n):
    degrees=[1]*n
    for v in code:
        degrees[v]+=1
    edges=[]
    for v in code:
        leaf=next(i for i,d in enumerate(degrees) if d==1)
        edges.append((leaf,v));degrees[leaf]-=1;degrees[v]-=1
    remaining=[i for i,d in enumerate(degrees) if d==1]
    edges.append(tuple(remaining))
    return edges

def path(adj,a,b):
    queue=collections.deque([a]);prev={a:None}
    while b not in prev:
        u=queue.popleft()
        for v in adj[u]:
            if v not in prev:
                prev[v]=u;queue.append(v)
    ans={b}
    while b!=a:
        b=prev[b];ans.add(b)
    return ans

trees=triples=0
for n in range(2,7):
    for code in itertools.product(range(n),repeat=n-2):
        adj=[set() for _ in range(n)]
        for a,b in tree_edges(code,n):
            adj[a].add(b);adj[b].add(a)
        paths={(a,b):path(adj,a,b) for a in range(n) for b in range(n)}
        for a,b,c in itertools.combinations_with_replacement(range(n),3):
            med=paths[a,b]&paths[a,c]&paths[b,c]
            require(len(med)==1,('unique median',n,code,a,b,c))
            m=next(iter(med));arms=[paths[m,v] for v in (a,b,c)]
            require(all(arms[i]&arms[j]=={m} for i,j in ((0,1),(0,2),(1,2))),('independent arms',n,code,a,b,c))
            triples+=1
        trees+=1

# Eisenstein unit quotient: each local step has exactly one lift, but cycles may fail.
units={(1,0),(-1,0),(0,1),(0,-1),(1,1),(-1,-1)}
q=lambda p:(p[0]+2*p[1])%7
by_residue={q(u):u for u in units}
require(set(by_residue)==set(range(1,7)),'all K7 steps lifted uniquely')
require(all((a+b)%3!=0 and a*a-a*b+b*b==1 for a,b in units),'unit graph proper three-coloring')
counts=collections.Counter()
for a,b in itertools.permutations(range(1,7),2):
    steps=[by_residue[a],by_residue[(b-a)%7],by_residue[-b%7]]
    displacement=tuple(sum(v[i] for v in steps) for i in (0,1))
    require(q(displacement)==0,'quotient closes')
    counts['exactly_closed' if displacement==(0,0) else 'only_modularly_closed']+=1
require(dict(counts)=={'only_modularly_closed':18,'exactly_closed':12},'full rooted triangle census')
require(tuple(sum(by_residue[r][i] for r in (1,1,5)) for i in (0,1))==(2,-1),'explicit failed triangle')
# Multiplication by 2-omega has matrix [[2,1],[-1,3]]: determinant 7,
# its two columns lie in ker q, proving equality of index-seven subgroups.
require(2*3-(-1)*1==7 and q((2,-1))==q((1,3))==0,'kernel ideal index')

# Endpoint equation alone does not preserve reachability.
swap=lambda v:1-v
require(iterate(swap,0,2)==0 and swap(0)==1,'endpoint-only hostile')

# Canonical suspensions exhaustively preserve first-return times and periods.
suspensions=0
for n in range(1,5):
    for f in itertools.product(range(n),repeat=n):
        for clocks in itertools.product((1,2,3),repeat=n):
            def route(state):
                x,j=state
                return (x,j+1) if j+1<clocks[x] else (f[x],0)
            for x in range(n):
                state=(x,0)
                for t in range(1,clocks[x]+1):
                    state=route(state)
                    require((state[1]==0)==(t==clocks[x]),'first-return route has no early marker')
                require(state==(f[x],0),'first-return endpoint')
                orbit=[];y=x
                while y not in orbit:
                    orbit.append(y);y=f[y]
                cycle_start=orbit.index(y)
                expected_period=sum(clocks[z] for z in orbit[cycle_start:])
                r_orbit=[];v=(x,0)
                while v not in r_orbit:
                    r_orbit.append(v);v=route(v)
                require(len(r_orbit)-r_orbit.index(v)==expected_period,'least eventual period cocycle')
                require((cycle_start==0)==(r_orbit.index(v)==0),'periodic iff periodic')
            suspensions+=1

# Critical source identities are explicit and checked, not merely labels.
pins=json.loads((ROOT/'independent/SOURCE_PINS.json').read_text())
for entry in pins['local_packets']:
    data=(ROOT/entry['path']).read_bytes()
    require(len(data)==entry['bytes'] and hashlib.sha256(data).hexdigest()==entry['sha256'],('source pin',entry['path']))

print(json.dumps({'status':'PASS','checks':CHECKS,'scope':'Independent finite exact controls only; written audit handles infinite claims; no convergence conclusion or external HN/Lean audit',
 'baseline':pins['repository_baseline'],'gate_table':residue_gates,'full_fibres':fibres,'sharp_clock_full_fibres':{'C_inverse4_of7':[18,112],'T_inverse3_of7':[18,56]},
 'arithmetic_color_samples':len(large_inputs),'large_integer_boundary_max_exponent':1024,
 'selected_trees':macro_receipts,'tree_median_control':{'labeled_trees':trees,'attachment_triples':triples},
 'rooted_K7_triangles':dict(counts),'canonical_suspensions':suspensions,
 'verified_local_source_pins':len(pins['local_packets'])},sort_keys=True,indent=2))
