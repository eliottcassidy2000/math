#!/usr/bin/env python3
"""Positive two-apex transfer, strong-span certificate and complete suffix menus.

The gap decomposition and rank17 are inherited. This independent 24-state
implementation merges the empty/good first gap; no minimality is claimed.
"""
from collections import defaultdict, deque
from fractions import Fraction
from itertools import product, permutations
from math import prod
import hashlib
import json
from pathlib import Path

HERE=Path(__file__).resolve().parent
STATES=tuple(product(range(2),range(2),range(3),range(2)))
INDEX={s:i for i,s in enumerate(STATES)}
INITIAL=tuple(int(p==1 and q==0 and r==0) for o,p,q,r in STATES)
OUTPUT=tuple(int(p==1 and q<2 and (o==0 or q==1)) for o,p,q,r in STATES)
CONSTANT=tuple(int(o==0 and q==0 and r==0) for o,p,q,r in STATES)
SUFFIXES=('','0','1','2','00','10','20','30','01','11','21','31',
          '000','100','200','010','210')
STRONG_BASIS=('2','12','20','21','22','30','32','102','112','120',
              '121','122','130','132','200','201','210')

def need(ok,message):
    if not ok: raise RuntimeError(message)

def transitions(c):
    edges=[]
    for j,(o,p,q,r) in enumerate(STATES):
        a=(c>>(1 if o else 0))&1
        b=(c>>(0 if o else 1))&1
        edges.append((INDEX[o,1-a,q,r],j))
        if q or a: edges.append((INDEX[o,p,2 if b else 1,r],j))
        if r or b: edges.append((INDEX[o,p,q,1],j))
    return tuple(edges)

EDGES=tuple(transitions(c) for c in range(4))

def step(v,c):
    out=[0]*24
    for i,j in EDGES[c]: out[i]+=v[j]
    return tuple(out)

def adjoint(v,c):
    out=[0]*24
    for i,j in EDGES[c]: out[j]+=v[i]
    return tuple(out)

def dot(a,b): return sum(x*y for x,y in zip(a,b))

def state(w):
    v=INITIAL
    for c in w: v=step(v,int(c))
    return v

def suffix_row(w):
    row=OUTPUT
    for c in w[::-1]: row=adjoint(row,int(c))
    return row

ROWS=tuple(suffix_row(s) for s in SUFFIXES)
def profile(w):
    v=state(w)
    return tuple(dot(row,v) for row in ROWS)
def hamilton(w): return dot(OUTPUT,state(w))

def graph(w):
    n=len(w)+2;rows=[0]*n;rows[0]=2
    for i,c in enumerate(map(int,w),2):
        for a in range(2):
            if c>>a&1: rows[a]|=1<<i
            else: rows[i]|=1<<a
        for j in range(2,i): rows[j]|=1<<i
    return tuple(rows)

def direct_h(rows):
    n=len(rows);D=[[0]*n for _ in range(1<<n)]
    for v in range(n): D[1<<v][v]=1
    for mask in range(1,1<<n):
        for v in range(n):
            if not mask>>v&1: continue
            before=mask^(1<<v)
            if before: D[mask][v]=sum(D[before][u] for u in range(n)
                                      if before>>u&1 and rows[u]>>v&1)
    return sum(D[-1])

def strong_rows(rows):
    n=len(rows);full=(1<<n)-1
    for v in range(n):
        seen=1<<v
        while True:
            nxt=seen
            for u in range(n):
                if seen>>u&1: nxt|=rows[u]
            if nxt==seen: break
            seen=nxt
        if seen!=full: return False
    return True

def strong_word(w):
    if not w or w[0]=='0' or w[-1]=='3': return False
    seen=False
    for c in map(int,w):
        seen=seen or bool(c&2)
        if seen and not c&1: return True
    return False

def canon(rows):
    groups=defaultdict(list)
    for v,row in enumerate(rows): groups[row.bit_count()].append(v)
    candidates=[]
    for chunks in product(*(permutations(groups[d]) for d in sorted(groups))):
        p=sum(chunks,())
        candidates.append(tuple(sum(((rows[a]>>b)&1)<<j for j,b in enumerate(p)) for a in p))
    return min(candidates)

def floor_strong(n):
    need(n>=3,'strong floor domain')
    k,r=divmod(n,3)
    return (3*5**(k-1),5**k,9*5**(k-1))[r]

def max_strong(h):
    n=3
    while floor_strong(n+1)<=h: n+=1
    return n

def suffix_menu(w,targets):
    need(strong_word(w),'suffix menu requires a strong source')
    bound=max_strong(max(targets))-len(w)-2
    found={};visited=0
    def visit(v,s):
        nonlocal visited
        visited+=1;h=dot(OUTPUT,v)
        if h in targets: found.setdefault(h,s)
        if len(s)==bound or h>max(targets): return
        for c in range(4): visit(step(v,c),s+str(c))
    visit(state(w),'')
    return found,bound,visited

def insert_basis(basis,v):
    v=list(map(Fraction,v))
    for p,row in basis.items():
        if v[p]:
            a=v[p];v=[x-a*y for x,y in zip(v,row)]
    for p,x in enumerate(v):
        if x:
            basis[p]=tuple(y/x for y in v)
            return True
    return False

def solve(matrix,rhs):
    n=len(rhs);a=[list(map(Fraction,row))+[Fraction(b)] for row,b in zip(matrix,rhs)]
    determinant=Fraction(1)
    for j in range(n):
        k=next((i for i in range(j,n) if a[i][j]),None)
        need(k is not None,'singular certificate matrix')
        if k!=j: a[k],a[j]=a[j],a[k];determinant=-determinant
        pivot=a[j][j];determinant*=pivot;a[j]=[v/pivot for v in a[j]]
        for i in range(n):
            if i!=j and a[i][j]:
                q=a[i][j];a[i]=[x-q*y for x,y in zip(a[i],a[j])]
    return tuple(row[-1] for row in a),determinant

def algebra_certificate():
    # Close the physical reachable span rather than guessing a word-depth cutoff.
    basis={};vectors=[];words=[];queue=deque([(INITIAL,'')])
    while queue:
        v,w=queue.popleft()
        if not insert_basis(basis,v): continue
        vectors.append(v);words.append(w)
        for c in range(4): queue.append((step(v,c),w+str(c)))
    strong_matrix=[profile(w) for w in STRONG_BASIS]
    need(all(strong_word(w) for w in STRONG_BASIS),'strong spanning witnesses')
    constant_coeff,det=solve(strong_matrix,[1]*17)
    need(all(adjoint(CONSTANT,c)==CONSTANT for c in range(4)),'constant row invariant')
    need(dot(CONSTANT,INITIAL)==1,'constant initial value')
    def represent(row):
        coeff,_=solve(strong_matrix,[dot(row,state(w)) for w in STRONG_BASIS])
        difference=tuple(sum(coeff[j]*ROWS[j][i] for j in range(17))-row[i] for i in range(24))
        need(all(dot(difference,v)==0 for v in vectors),'all-reachable row identity')
        return coeff
    need(represent(CONSTANT)==constant_coeff,'constant representation')
    for row in ROWS:
        for c in range(4): represent(adjoint(row,c))
    return {'physical_reachable_dimension':len(vectors),'observable_dimension':17,
            'strong_basis':list(STRONG_BASIS),'strong_basis_determinant':str(det),
            'constant_coefficients':[str(x) for x in constant_coeff],
            'closure_identities':68,'reachable_basis_words':words}

def main():
    algebra=algebra_certificate();literal=0;groups_by_order={};allmenus={}
    for m in range(5):
        groups=defaultdict(list)
        for cs in product('0123',repeat=m):
            w=''.join(cs);rows=graph(w)
            need(hamilton(w)==direct_h(rows),'independent subset DP')
            need(strong_word(w)==strong_rows(rows),'strong language criterion')
            literal+=1
            if strong_word(w): groups[canon(rows)].append(w)
        if not m: continue
        entries=[]
        for rows,words in sorted(groups.items()):
            h=hamilton(words[0]);targets={3*h+2,3*h+4};found={};menus=[]
            for w in words:
                actual,bound,visits=suffix_menu(w,targets)
                for target,suffix in actual.items(): found.setdefault(target,{'word':w,'suffix':suffix})
                menus.append({'word':w,'bound':bound,'found':actual,'states_visited':visits})
            entry={'rows':rows,'h':h,'words':words,'found':found,
                   'missing':sorted(targets-set(found)),'menus':menus}
            entries.append(entry)
            for w in words: allmenus[w]=entry
        groups_by_order[m+2]=entries
    need(allmenus['21']['missing']==[],'remarking repairs order4 word21')
    need(19 not in suffix_menu('21',{17,19})[0],'fixed word21 misses19')
    need(17 not in suffix_menu('32',{17,19})[0],'fixed word32 misses17')
    h15=[entry for entries in groups_by_order.values() for entry in entries if entry['h']==15]
    need(len(h15)==4,'four distinct strong H15 classes')
    need(all(47 in entry['missing'] for entry in h15),'all H15 classes miss47 in every mark')
    # Independent full labelled universe supplies completeness of this H15 list.
    independent=set();labelled_checked=0
    for n in (5,6):
        pairs=[(i,j) for i in range(n) for j in range(i+1,n)]
        # Every class has the Hamilton frame 0->1->...->n-1.
        free=[(i,j) for i,j in pairs if j!=i+1]
        for bits in range(1<<len(free)):
            rows=[0]*n
            for i in range(n-1): rows[i]|=1<<(i+1)
            for k,(i,j) in enumerate(free):
                a,b=(i,j) if bits>>k&1 else (j,i);rows[a]|=1<<b
            labelled_checked+=1
            if strong_rows(rows) and direct_h(rows)==15: independent.add(canon(rows))
    need(independent=={tuple(e['rows']) for e in h15},'independent all-class H15 recovery')
    # Positive controls: the absent extensions are actual attainable targets.
    positive={'29':'1202','31':'1210','47':'11122','49':'10212'}
    for target,w in positive.items(): need(direct_h(graph(w))==int(target),'positive target control')
    # One interior-insertion repair, with all old vertices and apices retained.
    need(hamilton('112')==9 and direct_h(graph('10102'))==31,'interior insertion repair')
    need(31 in allmenus['112']['missing'],'all-marking obstruction at the interior repair source')
    # The proposed repeating repair has a different slope; certify all lengths.
    family_basis={};family_vectors=[];v=state('1010')
    while insert_basis(family_basis,v):
        family_vectors.append(v);v=step(v,1)
    row=suffix_row('2')
    for v in family_vectors:
        need(dot(row,step(step(v,1),1))-3*dot(row,step(v,1))+2*dot(row,v)==0,
             'closed (shift-1)(shift-2) family identity')
    need(hamilton('10102')==31 and hamilton('101012')==53,'family initial values')
    out={'status':'PASS','scope':'All-length suffix menus from every two-apex marking of each listed strong source; full spectrum remains open.',
         'algebra':algebra,'literal_words':literal,'independent_H15_frames':labelled_checked,
         'class_summary':{n:{'classes':len(entries),'both_branches':sum(not e['missing'] for e in entries)} for n,entries in groups_by_order.items()},
         'classes':groups_by_order,'positive_target_words':positive,
         'interior_repair':{'source':'112','target':'10102','old_chain_positions':[0,2,4],'H':[9,31]},
         'repeated_interior_family':{'formula':'11*2^(m-2)+9, m>=3','Krylov_dimension':len(family_vectors),'initial':[31,53]},
         'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    print(json.dumps(out,indent=2,sort_keys=True))

if __name__=='__main__': main()
