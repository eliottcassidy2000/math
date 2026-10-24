#!/usr/bin/env python3
"""Independent matrix/canonical-orbit audit; does not import producer code."""
import argparse
from collections import Counter
from itertools import combinations, permutations, product
from math import comb, factorial
from pathlib import Path
import json

CHECKS=Counter()
def require(value, label):
    CHECKS[label]+=1
    if not value: raise RuntimeError(label)

def decode(n, code):
    a=[[0]*n for _ in range(n)]
    for k,(i,j) in enumerate(combinations(range(n),2)):
        a[i][j]=(code>>k)&1; a[j][i]=1-a[i][j]
    return a

def code(a, p=None):
    if p is None: p=range(len(a))
    return sum(a[p[i]][p[j]]<<k for k,(i,j) in enumerate(combinations(range(len(a)),2)))

def stable(a):
    cells=[tuple(range(len(a)))] if a else []
    while True:
        new=[]
        for cell in cells:
            bins={}
            for v in cell:
                key=tuple(sum(a[v][w] for w in c) for c in cells)
                bins.setdefault(key,[]).append(v)
            new.extend(tuple(bins[k]) for k in sorted(bins))
        if len(new)==len(cells): return tuple(new)
        cells=new

def literal_group(a):
    """Full factorial permutations, without degree or refinement restrictions."""
    target=code(a)
    return tuple(p for p in permutations(range(len(a))) if code(a,p)==target)

def vertex_orbits(group,n):
    return tuple(sorted(set(tuple(sorted({p[v] for p in group})) for v in range(n))))

def excess(a, group):
    cs={v:len(c) for c in stable(a) for v in c}
    return sum(cs[o[0]]-len(o) for o in vertex_orbits(group,len(a)))

def strong_components(a):
    n=len(a); reachable=[]
    for root in range(n):
        seen={root}; todo=[root]
        while todo:
            x=todo.pop()
            for y in range(n):
                if a[x][y] and y not in seen: seen.add(y); todo.append(y)
        reachable.append(seen)
    return tuple(sorted(set(tuple(v for v in range(n) if v in reachable[u] and u in reachable[v]) for u in range(n))))

def induced(a,vertices): return [[a[v][w] for w in vertices] for v in vertices]

def check_refinement(a,group):
    cells=stable(a); comps=strong_components(a)
    for c in cells:
        require(len(c)%2==1,'stable cell odd')
        require(all(sum(a[v][w] for w in c)==(len(c)-1)//2 for v in c),'stable cell regular')
        require(all({p[v] for v in c}==set(c) for p in group),'stable cells setwise fixed')
    rebuilt=[]; local_error=0
    for component in comps:
        b=induced(a,component)
        rebuilt.extend(tuple(component[v] for v in c) for c in stable(b))
        local_error+=excess(b,literal_group(b))
    require(set(cells)==set(rebuilt),'SCC refinement restriction')
    require(excess(a,group)==local_error,'SCC excess additivity')

def full_classes(n, universe=None):
    remaining=set(range(1<<comb(n,2))) if universe is None else set(universe)
    initial=len(remaining); result=[]; ps=tuple(permutations(range(n)))
    while remaining:
        representative=min(remaining); a=decode(n,representative)
        orbit={code(a,p) for p in ps}
        require(orbit<=remaining,'disjoint full permutation class removal')
        group=tuple(p for p in ps if code(a,p)==representative)
        require(len(orbit)*len(group)==factorial(n),'literal orbit stabilizer')
        result.append((a,group)); remaining.difference_update(orbit)
    return initial,result

def regular_codes(n):
    target=(n-1)//2; found=set()
    def visit(v,a):
        if v==n:
            require(all(sum(row)==target for row in a),'regular leaf')
            found.add(code(a)); return
        need=target-sum(a[v]); candidates=range(v+1,n)
        if need<0 or need>n-v-1: return
        for wins in combinations(candidates,need):
            wins=set(wins); b=[row[:] for row in a]
            for w in candidates: b[v][w]=int(w in wins); b[w][v]=1-b[v][w]
            if all(sum(b[w])<=target<=sum(b[w])+n-v-2 for w in candidates): visit(v+1,b)
    visit(0,[[0]*n for _ in range(n)])
    return found

def family(r, sigmas):
    a=[[0]*(3*r) for _ in range(3*r)]
    for i in range(r):
        for x in range(3): a[3*i+x][3*i+(x+1)%3]=1
    for (i,j),sigma in zip(combinations(range(r),2),sigmas):
        for x,y in product(range(3),repeat=2):
            u,v=3*i+x,3*j+y; a[u][v]=int(y==sigma[x]); a[v][u]=1-a[u][v]
    return a

def degree_maps(a):
    deg=[sum(row) for row in a]; cells=[tuple(v for v in range(len(a)) if deg[v]==d) for d in sorted(set(deg))]
    for ps in product(*(permutations(c) for c in cells)):
        p=list(range(len(a)))
        for c,q in zip(cells,ps):
            for v,w in zip(c,q): p[v]=w
        yield tuple(p)

def degree_group_and_key(a):
    c=code(a); group=[]; key=None
    for p in degree_maps(a):
        b=code(a,p)
        if b==c: group.append(p)
        if key is None or b<key: key=b
    return tuple(group),key

def signed_kernel(sigmas,r):
    edges=tuple(combinations(range(r),2)); signs=[]
    for sigma in sigmas:
        inv=sum(sigma[i]>sigma[j] for i,j in combinations(range(3),2))
        signs.append(-1 if inv%2 else 1)
    values=tuple(k for k in product(range(3),repeat=r) if all((k[j]-e*k[i])%3==0 for (i,j),e in zip(edges,signs)))
    return tuple(signs),values

def join(parts):
    total=sum(len(p) for p in parts); a=[[0]*total for _ in range(total)]; off=0
    for part in parts:
        for v in range(len(part)):
            for w in range(len(part)): a[off+v][off+w]=part[v][w]
            for w in range(off+len(part),total): a[off+v][w]=1
        off+=len(part)
    return a

def main():
    p=argparse.ArgumentParser(); p.add_argument('--receipt',required=True); args=p.parse_args()
    small={}; reps={}
    for n in range(1,7):
        labelled,classes=full_classes(n); reps[n]=classes
        for a,g in classes:
            check_refinement(a,g); require(excess(a,g)==0,'all classes through six zero excess')
        small[n]={'labels':labelled,'classes':len(classes)}
    require([small[n]['classes'] for n in range(1,7)]==[1,1,2,4,12,56],'small class totals')
    regs={}
    for n in (5,7):
        labelled,classes=full_classes(n,regular_codes(n))
        records=[]
        for a,g in classes:
            check_refinement(a,g)
            records.append({'code':code(a),'automorphism_order':len(g),'orbits':vertex_orbits(g,n),'error':excess(a,g),'strong':len(strong_components(a))==1})
        regs[n]={'labels':labelled,'classes':records}
    require(regs[5]['labels']==24 and len(regs[5]['classes'])==1,'unique regular five')
    require(regs[7]['labels']==2640 and len(regs[7]['classes'])==3,'regular seven census')
    require(sorted(c['automorphism_order'] for c in regs[7]['classes'])==[3,7,21],'regular seven group orders')
    require(sum(c['error'] for c in regs[7]['classes'])==14,'regular seven excess fourteen')
    require(all(c['strong'] for c in regs[7]['classes']),'regular seven all strong')
    perms3=tuple(permutations(range(3))); classes={}; signs_count=Counter()
    for sigmas in product(perms3,repeat=3):
        a=family(3,sigmas); g,key=degree_group_and_key(a); signs,ks=signed_kernel(sigmas,3)
        expected=[]
        for shifts in ks: expected.append(tuple(3*i+(x+shifts[i])%3 for i in range(3) for x in range(3)))
        require(set(g)==set(expected),'r3 literal group equals signed kernel')
        require(stable(a)==((0,1,2),(3,4,5),(6,7,8)),'r3 distinct degree cells')
        require(len(strong_components(a))==1,'r3 all strong')
        require(excess(a,g)==(18 if len(g)==1 else 0),'r3 error formula')
        record={'automorphism_order':len(g),'error':excess(a,g),'signs':signs}
        require(key not in classes or classes[key]==record,'r3 canonical class invariant')
        classes[key]=record; signs_count[signs]+=1
    dist=Counter(c['automorphism_order'] for c in classes.values())
    require(len(classes)==16 and dist=={1:4,3:12},'r3 sixteen classes split twelve four')
    require(sum(c['error'] for c in classes.values())==72,'r3 native excess seventy two')
    # Exact gauge-orbit enumeration on the entire r=4 affine family (6^6=46656).
    # No tournament isomorphism assumption is used in the orbit removal itself.
    r4_counts=Counter(); edges=tuple(combinations(range(4),2)); powers=tuple(product(range(3),repeat=4))
    for signs in product((1,-1),repeat=6):
        offsets=set(product(range(3),repeat=6)); class_count=0; kernel=None
        while offsets:
            t=min(offsets)
            orbit={tuple((t[h]+k[j]-signs[h]*k[i])%3 for h,(i,j) in enumerate(edges)) for k in powers}
            require(orbit<=offsets,'r4 gauge orbits disjoint')
            offsets.difference_update(orbit); class_count+=1
            size=81//len(orbit)
            require(kernel is None or kernel==size,'r4 kernel offset independent'); kernel=size
        r4_counts[(kernel,class_count)]+=1
    require(r4_counts=={(3,27):8,(1,9):56},'r4 exact gauge quotient')
    r4_total=sum(number*classes for (kernel,classes),number in r4_counts.items())
    require(r4_total==720,'r4 family classes seven twenty')
    # Literal degree-preserving maps for two order12 controls, no gauge pruning.
    control12=[]
    for sigmas in ((perms3[0],)*6,(perms3[0],)*5+(perms3[1],)):
        a=family(4,sigmas); g,key=degree_group_and_key(a); signs,ks=signed_kernel(sigmas,4)
        require(len(g)==len(ks),'r4 degree-only literal group control')
        control12.append({'code':key,'automorphism_order':len(g),'error':excess(a,g)})
    hostile7=next(decode(7,c['code']) for c in regs[7]['classes'] if c['error'])
    joins=[]
    # Degree blocks prevent factorial blowup here; SCC decomposition is checked independently.
    for sizes in ((1,7),(7,1),(1,1,7),(1,7,1),(7,1,1)):
        a=join([hostile7 if size==7 else [[0]] for size in sizes]); g,_=degree_group_and_key(a)
        require(excess(a,g)==14,'hostile7 sequence inheritance')
        require(sorted(map(len,stable(a)))==sorted(sizes),'hostile7 sequence refinement cells')
        joins.append({'sizes':sizes,'error':excess(a,g)})
    prime_families={}
    for prime in (3,5,7):
        remaining=set(permutations(range(prime))); distribution=Counter()
        while remaining:
            q=min(remaining)
            orbit={tuple((q[(x-alpha)%prime]+beta)%prime for x in range(prime))
                   for alpha,beta in product(range(prime),repeat=2)}
            require(orbit<=remaining,'prime double rotation orbit disjoint')
            remaining.difference_update(orbit)
            stabilizer=prime*prime//len(orbit)
            distribution[stabilizer]+=1
        ordinary=(factorial(prime-1)+(prime-1)**2)//prime
        require(sum(distribution.values())==ordinary,'prime family ordinary count')
        require(distribution[prime]==prime-1,'prime affine symmetric class count')
        contribution=2*prime*(prime-1)*distribution[1]
        prime_families[prime]={'classes':ordinary,'group_order_distribution':dict(sorted(distribution.items())),'strong_error_contribution':contribution}
    require(prime_families[5]['strong_error_contribution']==160,'strong ten lower bound')
    require(prime_families[7]['strong_error_contribution']==8568,'strong fourteen lower bound')
    require(160+2*72+6*14==388,'total ten lower bound convolution')
    report={'schema':'session50-independent-cell-holonomy-v1','small_census':small,'regular_census':regs,
            'r3_family':{'labelled_configurations':216,'classes':len(classes),'group_order_distribution':dict(sorted(dist.items())),'total_excess':72},
            'r4_family':{'affine_configurations':6**6,'classes':r4_total,'gauge_distribution':[{'automorphism_order':k,'classes_per_signing':c,'signings':v} for (k,c),v in sorted(r4_counts.items())],'literal_controls':control12},
            'two_prime_matching_families':prime_families,'SCC_hostile7_inheritance':joins,'derived':{'Estrong7':14,'Estrong8':0,'E8':28,'Estrong9_lower_bound':72,'E9_lower_bound':114,'Estrong10_lower_bound':160,'E10_lower_bound':388,'Estrong14_lower_bound':8568},'checks':dict(sorted(CHECKS.items()))}
    Path(args.receipt).write_text(json.dumps(report,indent=2,sort_keys=True)+'\n')
    print(json.dumps({'status':'PASS','checks':sum(CHECKS.values()),'small_classes':[small[n]['classes'] for n in range(1,7)],'r3_classes':len(classes),'r4_classes':r4_total},sort_keys=True))
if __name__=='__main__': main()
