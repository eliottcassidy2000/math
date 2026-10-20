#!/usr/bin/env python3
"""Complete odd non-prime-power subgroup census in degrees <=10.
The completeness reduction is proved in README.md; the tiny prime-degree
transitive classification is independently exhausted here. Standard library
only. Full isomorphism and global signs are never replaced by N(H)-equivalence.
All validation survives python -O.
"""
from itertools import permutations, combinations, product
from collections import Counter, defaultdict
from functools import lru_cache
import json
CHECKS=0

def need(test, message):
    global CHECKS
    CHECKS += 1
    if not test: raise RuntimeError(message)

def compose(p,q): return tuple(p[q[i]] for i in range(len(p)))
def inverse(p): return tuple(p.index(i) for i in range(len(p)))
def closure(gens,n):
    found={tuple(range(n))}; todo=list(found)
    while todo:
        x=todo.pop()
        for g in gens:
            y=compose(g,x)
            if y not in found: found.add(y); todo.append(y)
    return frozenset(found)
def cyc(n,cycles):
    p=list(range(n))
    for c in cycles:
        for a,b in zip(c,c[1:]+c[:1]): p[a]=b
    return tuple(p)
def pair_data(n,H):
    pairs=list(combinations(range(n),2)); index={e:i for i,e in enumerate(pairs)}
    remaining=set(pairs); orbits=[]; base=0
    while remaining:
        u,v=min(remaining); arcs={(g[u],g[v]) for g in H}; orb={tuple(sorted(e)) for e in arcs}
        need(len(arcs)==len(orb),'odd group reverses pair')
        orbits.append(sum(1<<index[e] for e in orb)); remaining-=orb
        base|=sum(1<<index[tuple(sorted((a,b)))] for a,b in arcs if a>b)
    return pairs,index,orbits,base

def masks(orbits,base=0):
    return [base ^ sum(o for i,o in enumerate(orbits) if s>>i&1) for s in range(1<<len(orbits))]

def normalizer(n,H,gens):
    left=set(range(n)); blocks=[]
    while left:
        u=min(left); b=sorted({g[u] for g in H}); blocks.append(b);left-=set(b)
    # In this proved list all nonsingleton orbit sizes are distinct.
    roots=[b[0] for b in blocks if len(b)==1]
    blocks=[b for b in blocks if len(b)>1]+([roots] if roots else [])
    need(len({len(b) for b in blocks if len(b)>1})==len([b for b in blocks if len(b)>1]) or n==10,'orbit partition convention')
    ans=[]
    for choices in product(*(list(permutations(b)) for b in blocks)):
        p=list(range(n))
        for b,c in zip(blocks,choices):
            for x,y in zip(b,c): p[x]=y
        p=tuple(p); ip=inverse(p)
        if all(compose(compose(p,g),ip) in H for g in gens): ans.append(p)
    return ans

def prep(p,pairs,index):
    images=[];rev=destrev=0
    for k,(u,v) in enumerate(pairs):
        a,b=p[u],p[v];j=index[tuple(sorted((a,b)))];images.append(1<<j)
        if a>b: rev|=1<<k;destrev|=1<<j
    return images,rev,destrev

def move(mask,r,tournament=False):
    out=0;x=mask
    while x:
        b=x&-x;out|=r[0][b.bit_length()-1];x-=b
    return out^(r[2] if tournament else 0)
def sign(mask,r): return -1 if (mask&r[1]).bit_count()%2 else 1

def matrix(mask,n,pairs,tournament=False):
    A=[[0]*n for _ in range(n)]
    for k,(u,v) in enumerate(pairs):
        bit=(mask>>k)&1
        if tournament: A[u][v]=1-bit; A[v][u]=bit
        else: A[u][v]=A[v][u]=bit
    return A

def isomorphisms(A,B,colorsA=None,colorsB=None):
    n=len(A)
    if len(B)!=n:return
    ca=colorsA[:] if colorsA is not None else [sum(r) for r in A]
    cb=colorsB[:] if colorsB is not None else [sum(r) for r in B]
    # Simultaneous stable equitable refinement, valid for directed and undirected.
    while True:
        palette=sorted(set(ca+cb)); da=[(ca[i],tuple(sum(A[i][j] for j in range(n) if ca[j]==c) for c in palette),tuple(sum(A[j][i] for j in range(n) if ca[j]==c) for c in palette)) for i in range(n)]
        db=[(cb[i],tuple(sum(B[i][j] for j in range(n) if cb[j]==c) for c in palette),tuple(sum(B[j][i] for j in range(n) if cb[j]==c) for c in palette)) for i in range(n)]
        keys={x:k for k,x in enumerate(sorted(set(da+db)))};na=[keys[x] for x in da];nb=[keys[x] for x in db]
        if Counter(na)!=Counter(nb):return
        stable=len(set(na))==len(set(ca)) and len(set(nb))==len(set(cb))
        ca,cb=na,nb
        if stable:break
    mapping={};used=set()
    def rec():
        if len(mapping)==n:
            yield tuple(mapping[i] for i in range(n));return
        best=None;choices=None
        for u in range(n):
            if u in mapping:continue
            opts=[v for v in range(n) if v not in used and ca[u]==cb[v] and all(A[u][x]==B[v][y] and A[x][u]==B[y][v] for x,y in mapping.items())]
            if not opts:return
            if choices is None or len(opts)<len(choices):best,choices=u,opts
        for v in choices:
            mapping[best]=v;used.add(v)
            yield from rec()
            del mapping[best];used.remove(v)
    yield from rec()

def iso(A,B): return next(isomorphisms(A,B),None) is not None

def full_graph_sign(A):
    n=len(A);remaining=set(range(n));classes=[]
    while remaining:
        u=min(remaining);C=tuple(v for v in sorted(remaining) if all(A[u][w]==A[v][w] for w in range(n) if w not in (u,v)))
        classes.append(C);remaining-=set(C)
    # A true-twin transposition is immediately odd.
    for C in classes:
        if len(C)>1 and A[C[0]][C[1]]:
            p=list(range(n));p[C[0]],p[C[1]]=C[1],C[0]
            return False,tuple(p),classes
    B=[[A[C[0]][D[0]] if i!=j else 0 for j,D in enumerate(classes)] for i,C in enumerate(classes)]
    colors=[len(C) for C in classes]
    for q in isomorphisms(B,B,colors,colors):
        p=list(range(n))
        for i,j in enumerate(q):
            for u,v in zip(classes[i],classes[j]):p[u]=v
        need(all(A[u][v]==A[p[u]][p[v]] for u in range(n) for v in range(n)),'twin quotient lift')
        if sum(A[u][v] for u,v in combinations(range(n),2) if p[u]>p[v])%2:
            return False,tuple(p),classes
    return True,None,classes

def relative_records(n,H,gens,tournament):
    pairs,index,orbits,base=pair_data(n,H);values=masks(orbits,base if tournament else 0)
    N=normalizer(n,H,gens);prepared=[prep(p,pairs,index) for p in N];seen=set();records=[]
    for x in sorted(values):
        if x in seen:continue
        copies={move(x,r,tournament) for r in prepared};seen|=copies
        stabilizers=[i for i,r in enumerate(prepared) if move(x,r,tournament)==x]
        relative_even=tournament or all(sign(x,prepared[i])==1 for i in stabilizers)
        if relative_even:records.append({'mask':x,'relative_stabilizer_order':len(stabilizers),'orbit_size':len(copies)})
    need(set(values)==seen,'normalizer orbits partition fixed cube')
    return records,pairs,orbits,N

def classify(records,n,pairs,tournament):
    reps=[];types=[]
    for r in records:
        A=matrix(r['mask'],n,pairs,tournament)
        j=next((i for i,B in enumerate(reps) if iso(A,B)),None)
        if j is None: j=len(reps);reps.append(A);types.append({'multiplicity':0,'representative_mask':r['mask']})
        types[j]['multiplicity']+=1
    return reps,types

def prime_degree_control(p):
    g=cyc(p,[list(range(p))]);H=closure([g],p);pairs,index,orbits,base=pair_data(p,H);out=[]
    for x in masks(orbits,base):
        A=matrix(x,p,pairs,True);auts=list(isomorphisms(A,A))
        # Independent brute-force all S_p permutations, not the refinement search.
        brute=[q for q in permutations(range(p)) if all(A[u][v]==A[q[u]][q[v]] for u,v in pairs)]
        need(set(auts)==set(brute),'independent prime-degree automorphism enumeration')
        affine={tuple((a*i+b)%p for i in range(p)) for a in range(1,p) for b in range(p)}
        need(set(auts)<=affine,'prime-degree automorphism group not affine')
        need(len(auts) in ([5] if p==5 else [7,21]),'unexpected prime-degree automorphism order')
        out.append({'mask':x,'automorphism_order':len(auts)})
    return {'degree':p,'all_cyclic_tournaments':len(out),'automorphism_order_histogram':dict(sorted(Counter(r['automorphism_order'] for r in out).items())),'records':out}

def cases():
    out=[]
    for n in range(7,11):
        t=cyc(n,[list(range(7))]);a=tuple([(2*i)%7 for i in range(7)]+list(range(7,n)))
        out.append((n,'F21_fixed', [t,a]))
        if n>=8:
            out.append((n,'C15_5plus3',[cyc(n,[list(range(5))]),cyc(n,[[5,6,7]])]))
        if n==10:
            b=cyc(n,[[7,8,9]])
            out.extend([(n,'C21_7plus3',[t,b]),(n,'F21_diagonal_7plus3',[t,compose(a,b)]),(n,'F21_times_C3',[t,a,b])])
    return out

def census_case(n,name,gens):
    H=closure(gens,n);T,pairs,orbits,N=relative_records(n,H,gens,True);E,_,_,N2=relative_records(n,H,gens,False)
    need(N==N2,'normalizer consistency');need(len(T)==len(E),'relative Walsh identity')
    tr,tt=classify(T,n,pairs,True);er,et=classify(E,n,pairs,False)
    even=[];escape=[]
    for A,ty in zip(er,et):
        good,witness,twins=full_graph_sign(A)
        ty['degrees']=sorted(map(sum,A));ty['twin_class_sizes']=sorted(map(len,twins))
        if good:even.append(ty)
        else:
            ty['odd_automorphism']=witness
            need(witness not in N,'escaping sign accidentally normalizes H')
            need(all(A[u][v]==A[witness[u]][witness[v]] for u in range(n) for v in range(n)),'odd witness not automorphism')
            need(sum(A[u][v] for u,v in combinations(range(n),2) if witness[u]>witness[v])%2==1,'odd witness has even sign')
            escape.append(ty)
    ft=sum(x['multiplicity']-1 for x in tt);fe=sum(x['multiplicity']-1 for x in even);z=sum(x['multiplicity'] for x in escape)
    need(len(tt)-len(even)==fe+z-ft,'defect balance')
    hist=lambda ts:dict(sorted(Counter(t['multiplicity'] for t in ts).items()))
    return {'n':n,'class':name,'order':len(H),'generators':[list(g) for g in gens],'normalizer_order':len(N),'pair_orbits':len(orbits),'fixed_objects_each':2**len(orbits),'relative_tournament_orbits':len(T),'relative_even_graph_orbits':len(E),'tau_T':len(tt),'tau_E':len(even),'difference':len(tt)-len(even),'tournament_multiplicity_histogram':hist(tt),'even_graph_multiplicity_histogram':hist(even),'escape_graph_multiplicity_histogram':hist(escape),'tournament_fusion_defect':ft,'even_graph_fusion_defect':fe,'sign_escape_relative_orbits':z,'sign_escape_global_types':len(escape),'tournament_types':tt,'even_graph_types':even,'sign_escape_types':escape}

def main():
    controls=[prime_degree_control(5),prime_degree_control(7)]
    rows=[census_case(n,name,gens) for n,name,gens in cases()]
    need(len(rows)==10,'complete universe has ten degree-class pairs')
    need([(r['n'],r['class']) for r in rows if r['difference']] == [(10,'C15_5plus3')],'unexpected equality failure')
    print(json.dumps({'status':'PASS','scope':'all odd non-prime-power subgroups H <= S_n, 1 <= n <= 10, up to S_n-conjugacy; completeness proof in README.md','checks':CHECKS,'prime_degree_controls':controls,'class_counts_by_degree':{str(n):sum(r['n']==n for r in rows) for n in range(1,11)},'rows':rows},indent=2,sort_keys=True))
if __name__=='__main__':main()
