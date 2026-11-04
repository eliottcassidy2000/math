#!/usr/bin/env python3
"""Independent 36-state gap observer and direct all-marking source audit.
No imports from the producer or other repository computation modules.
"""
from collections import deque, Counter
from fractions import Fraction
from functools import lru_cache
from itertools import product, permutations
import hashlib
import json
from pathlib import Path

HERE=Path(__file__).resolve().parent
SUFFIXES=('','0','1','2','00','10','20','30','01','11','21','31','000','100','200','010','210')
STRONG=('2','12','20','21','22','30','32','102','112','120','121','122','130','132','200','201','210')

def require(v,message):
    if not v:raise RuntimeError(message)

def graph(w):
    rows=[2,0]+[0]*len(w)
    for i,c in enumerate(w,2):
        for apex in (0,1):
            if int(c)&(1<<apex):rows[apex]|=1<<i
            else:rows[i]|=1<<apex
        for j in range(2,i):rows[j]|=1<<i
    return tuple(rows)

def hamilton(r):
    n=len(r)
    @lru_cache(None)
    def f(last,left):
        if not left:return 1
        return sum(f(v,left^(1<<v)) for v in range(n) if left>>v&1 and r[last]>>v&1)
    return sum(f(v,((1<<n)-1)^(1<<v)) for v in range(n))

def strong(r):
    n=len(r);reach=[v|1<<i for i,v in enumerate(r)]
    for k in range(n):
        for i in range(n):
            if reach[i]>>k&1:reach[i]|=reach[k]
    return all(v==(1<<n)-1 for v in reach)

def canon(r):
    n=len(r)
    cells={d:[i for i,v in enumerate(r) if v.bit_count()==d] for d in sorted(set(v.bit_count() for v in r))}
    best=None
    for blocks in product(*(permutations(v) for v in cells.values())):
        p=sum(blocks,());rep=tuple(sum(((r[u]>>v)&1)<<j for j,v in enumerate(p)) for u in p)
        if best is None or rep<best:best=rep
    return best

def frame(n,bits):
    r=[1<<(i+1) if i<n-1 else 0 for i in range(n)];k=0
    for i in range(n):
        for j in range(i+2,n):
            if bits>>k&1:r[i]|=1<<j
            else:r[j]|=1<<i
            k+=1
    return tuple(r)

def marks(r):
    n=len(r);out=set()
    for x in range(n):
        for y in range(n):
            if x==y or not(r[x]>>y&1):continue
            other=[i for i in range(n) if i not in (x,y)]
            original=tuple(other)
            other.sort(key=lambda i:sum(r[i]>>j&1 for j in original),reverse=True)
            if not all(r[u]>>v&1 for a,u in enumerate(other) for v in other[a+1:]):continue
            out.add(''.join(str(((r[x]>>v)&1)+2*((r[y]>>v)&1)) for v in other))
    return sorted(out)

def menu(r,targets,maximum_order):
    out={};cases=0;ww=marks(r)
    for w in ww:
        for length in range(maximum_order-len(r)+1):
            for cs in product('0123',repeat=length):
                s=''.join(cs);h=hamilton(graph(w+s));cases+=1
                if h in targets:out.setdefault(h,(w,s))
    return dict(words=ww,found=out,cases=cases)

# Full first-gap state remembers empty separately from good-last.
STATES=tuple(product(range(2),range(3),range(3),range(2)))
POS={v:i for i,v in enumerate(STATES)}
INIT=tuple(int(p==q==r==0) for o,p,q,r in STATES)
OUT=tuple(int(p!=2 and q!=2 and (o==0 or q==1)) for o,p,q,r in STATES)
ONE=tuple(int(o==0 and q==0 and r==0) for o,p,q,r in STATES)

def advance(v,c):
    ans=[0]*36
    for j,(o,p,q,r) in enumerate(STATES):
        val=v[j]
        if not val:continue
        a=(c>>(1 if o else 0))&1;b=(c>>(0 if o else 1))&1
        ans[POS[o,1 if not a else 2,q,r]]+=val
        if q or a:ans[POS[o,p,1 if not b else 2,r]]+=val
        if r or b:ans[POS[o,p,q,1]]+=val
    return tuple(ans)

def dot(a,b):return sum(x*y for x,y in zip(a,b))

def vector(w):
    v=INIT
    for c in w:v=advance(v,int(c))
    return v

def observed(v,s):
    for c in s:v=advance(v,int(c))
    return dot(OUT,v)

def independent_basis_insert(rows,v):
    w=list(map(Fraction,v))
    for pivot,row in rows:
        if w[pivot]:
            coeff=w[pivot];w=[a-coeff*b for a,b in zip(w,row)]
    if not any(w):return False
    p=next(i for i,x in enumerate(w) if x)
    rows.append((p,tuple(x/w[p] for x in w)))
    return True

def inverse_times(a,b):
    n=len(a);m=[list(map(Fraction,r))+[Fraction(t)] for r,t in zip(a,b)];det=Fraction(1)
    for j in range(n):
        k=next(i for i in range(j,n) if m[i][j])
        if k!=j:m[k],m[j]=m[j],m[k];det=-det
        pivot=m[j][j];det*=pivot
        for k in range(j,n+1):m[j][k]/=pivot
        for i in range(n):
            if i==j:continue
            factor=m[i][j]
            for k in range(j,n+1):m[i][k]-=factor*m[j][k]
    return tuple(r[-1] for r in m),det

def run():
    cases=0;word_classes={}
    for length in range(5):
        for chars in product('0123',repeat=length):
            w=''.join(chars);r=graph(w);require(observed(vector(w),'')==hamilton(r),('36 state graph count',w));cases+=1
            if strong(r):word_classes.setdefault(length+2,set()).add(canon(r))
    require({n:len(c) for n,c in word_classes.items()}=={3:1,4:1,5:6,6:34},'all small two-apex source class counts')
    reachable=[];basis=[];queue=deque([INIT])
    while queue:
        v=queue.popleft()
        if independent_basis_insert(basis,v):
            reachable.append(v)
            for c in range(4):queue.append(advance(v,c))
    vv=[vector(w) for w in STRONG]
    require(all(strong(graph(w)) for w in STRONG),'strong profile source')
    mat=[[observed(v,s) for s in SUFFIXES] for v in vv]
    # Literal graph re-evaluation of the complete17x17 certificate.
    require(mat==[[hamilton(graph(w+s)) for s in SUFFIXES] for w in STRONG],'literal17x17 profile')
    constant,det=inverse_times(mat,[1]*17)
    require(det==150994944,'strong determinant')
    require(all(dot(ONE,v)==1 for v in reachable),'independent full-gap constant')
    require(all(sum(c*observed(v,s) for c,s in zip(constant,SUFFIXES))==1 for v in reachable),'constant query representation')
    closure=0
    for suffix in SUFFIXES:
        for c in '0123':
            coefficients,_=inverse_times(mat,[observed(v,c+suffix) for v in vv])
            require(all(sum(k*observed(v,s) for k,s in zip(coefficients,SUFFIXES))==observed(v,c+suffix) for v in reachable),'36-state all-reachable closure')
            closure+=1
    sources=set();frames=0;all_strong={5:set(),6:set()}
    for n in (5,6):
        for bits in range(1<<((n-1)*(n-2)//2)):
            r=frame(n,bits);frames+=1
            if strong(r):
                all_strong[n].add(canon(r))
                if hamilton(r)==15:sources.add(canon(r))
    require(len(sources)==4,'four independent strong H15 cores')
    require({n:len(c) for n,c in all_strong.items()}=={5:6,6:35},'complete strong class controls')
    inventories=[];total_suffix=0
    for r in sorted(sources,key=lambda r:(len(r),r)):
        entry=menu(r,{47,49},8);total_suffix+=entry['cases']
        require(47 not in entry['found'],'every H15 marking misses47')
        if len(r)==5:require(49 not in entry['found'],'both H15 order5 markings miss49')
        inventories.append(dict(rows=r,order=len(r),**entry))
    require(sum(49 in e['found'] for e in inventories)==1,'exactly one H15 core reaches49')
    source9=canon(graph('112'));entry9=menu(source9,{31},7)
    require(31 not in entry9['found'],'H9 source all markings miss31')
    require(hamilton(graph('10102'))==31 and hamilton(source9)==9,'interior H9 to31')
    old=graph('112');new=graph('10102');embedding=(0,1,2,4,6)
    require(all(((old[i]>>j)&1)==((new[embedding[i]]>>embedding[j])&1) for i in range(5) for j in range(5)),'interior old arcs retained')
    weak_controls=[]
    for m in range(1,7):
        for c in '012':
            apex=0 if not(int(c)&1) else 1
            rr=graph('3'*m+c);vertices=[apex]+list(range(2,m+3))
            small=tuple(sum(1<<j for j,v in enumerate(vertices) if rr[u]>>v&1) for u in vertices)
            h=hamilton(small);require(h==2**m+1,'weak-source obstruction')
            weak_controls.append([m,c,h])
    positive={w:hamilton(graph(w)) for w in ('1202','1210','11122','10212')}
    require(list(positive.values())==[29,31,47,49],'positive targets')
    return dict(status='PASS',independent_full_gap_states=36,reachable_span=len(reachable),literal_word_checks=cases,
        strong_profile_determinant=str(det),closure_checks=closure,constant_coefficients=[str(c) for c in constant],
        independent_H15_frames=frames,word_source_classes={n:len(c) for n,c in word_classes.items()},
        complete_strong_classes={n:len(c) for n,c in all_strong.items()},H15_inventory=inventories,suffix_graph_counts=total_suffix,
        H9_all_marking31_hostile=entry9,interior_embedding=embedding,weak_unboundedness_controls=weak_controls,
        positive_targets=positive,
        producer_script_sha256=hashlib.sha256((HERE/'check.py').read_bytes()).hexdigest(),
        producer_receipt_sha256=hashlib.sha256((HERE/'check.json').read_bytes()).hexdigest(),
        independent_script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())

if __name__=='__main__':print(json.dumps(run(),indent=2,sort_keys=True))
