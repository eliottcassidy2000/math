#!/usr/bin/env python3
"""Independent Session145 bridge controls. Imports no discovery implementations."""
from collections import Counter
from itertools import combinations, permutations, product
import json, math

def check(ok, message):
    if not ok: raise RuntimeError(message)

def matchings(left):
    if not left:
        yield (); return
    a = left[0]
    for b in left[1:]:
        for tail in matchings(tuple(v for v in left if v != a and v != b)):
            yield ((a,b),) + tail

def gf2rank(rows):
    pivots={}
    for x in rows:
        while x:
            p=x.bit_length()-1
            if p not in pivots:
                pivots[p]=x;break
            x ^= pivots[p]
    return len(pivots)

def components(P,Q):
    # Rank of the unoriented incidence matrix over F2 equals |V|-components.
    return 2*len(P)-gf2rank([(1<<a)|(1<<b) for a,b in P+Q])

def perm_cycles(p):
    unseen=set(range(len(p))); k=0
    while unseen:
        x=next(iter(unseen)); k+=1
        while x in unseen:
            unseen.remove(x); x=p[x]
    return k

matching_rows=[]
for r in range(0,8):
    P=tuple((2*i,2*i+1) for i in range(r))
    M=list(matchings(tuple(range(2*r))))
    hist=Counter(components(P,Q) for Q in M)
    perm_hist=Counter(perm_cycles(p) for p in permutations(range(r)))
    expected={k: v*2**(r-k) for k,v in perm_hist.items()}
    check(dict(hist)==expected,'matching incidence rank / permutation-cycle coefficients')
    check(sum(hist.values())==math.prod(range(1,2*r,2)), 'matching total')
    if r:
        check(hist[1]==2**(r-1)*math.factorial(r-1),'one-component coefficient')
    unique=None
    if r<=4:
        rows=[tuple(components(A,B) for B in M) for A in M]
        check(len(set(rows))==len(M),'all contextual rows distinct')
        check(all((rows[i][j]==r)==(i==j) for i in range(len(M)) for j in range(len(M))), 'diagonal separator')
        unique=len(set(rows))
    matching_rows.append({'r':r,'coefficients_constant_first':[hist[k] for k in range(r+1)],'matchings':len(M),'contextual_rows':unique})

# Check the stronger actual union-of-two-linear-forests condition on SIX vertices.
n=6; E=list(combinations(range(n),2)); full=(1<<n)-1

def adj(mask):
    A=[0]*n
    for i,(u,v) in enumerate(E):
        if mask>>i&1:
            A[u]|=1<<v; A[v]|=1<<u
    return A

def linear_forest(A):
    if any(x.bit_count()>2 for x in A): return False
    todo=full
    while todo:
        reach=todo&-todo; frontier=reach
        while frontier:
            b=frontier&-frontier; frontier-=b
            new=A[b.bit_length()-1]&~reach
            reach|=new;frontier|=new
        vertices=[i for i in range(n) if reach>>i&1]
        if sum(A[i].bit_count() for i in vertices)//2 != len(vertices)-1: return False
        todo &= ~reach
    return True

LF=[mask for mask in range(1<<len(E)) if linear_forest(adj(mask))]
unions={a|b for a in LF for b in LF}
for mask in unions:
    A=adj(mask);remaining=full;order=[]
    while remaining:
        eligible=[i for i in range(n) if remaining>>i&1 and (A[i]&remaining).bit_count()<=3]
        check(bool(eligible),'two linear forests must be 3-degenerate')
        v=eligible[0];order.append(v);remaining &= ~(1<<v)
    colors={}
    for v in reversed(order):
        available=set(range(4))-{colors[w] for w in colors if A[v]>>w&1}
        check(bool(available),'4-color reverse greedy')
        colors[v]=min(available)

# Fox coloring: solve full cyclic crossing constraint systems, rather than powers.
def rank_mod3(A):
    A=[row[:] for row in A];k=0
    for j in range(len(A[0])):
        pivot=next((i for i in range(k,len(A)) if A[i][j]%3),None)
        if pivot is None: continue
        A[k],A[pivot]=A[pivot],A[k]
        s=pow(A[k][j]%3,-1,3);A[k]=[(s*x)%3 for x in A[k]]
        for i in range(len(A)):
            if i!=k:
                s=A[i][j];A[i]=[(a-s*b)%3 for a,b in zip(A[i],A[k])]
        k+=1
    return k
fox={}
for m in range(1,13):
    A=[]
    for t in range(m):
        nxt=(t+1)%m
        row=[0]*(2*m)
        for idx,value in [(2*nxt,1),(2*t,-2),(2*t+1,1)]: row[idx]=(row[idx]+value)%3
        A.append(row)
        row=[0]*(2*m)
        for idx,value in [(2*nxt+1,1),(2*t,-1)]:row[idx]=(row[idx]+value)%3
        A.append(row)
    fox[m]=3**(2*m-rank_mod3(A))
    check(fox[m]==(9 if m%3==0 else 3),'Fox full diagram constraints')
def cross(x,i):
    y=list(x);y[i],y[i+1]=(2*x[i]-x[i+1])%3,x[i];return tuple(y)
for x in product(range(3),repeat=3):
    check(cross(cross(cross(x,0),1),0)==cross(cross(cross(x,1),0),1),'braid relation all 27 colors')

# Enumerate actual words in the unique square-affine chart start=(0,1),
# then expand all 55 group elements and count every endpoint/arc/cycle.
q=11;QR={1,3,4,5,9}; canonical=[]
for tail in permutations(range(2,q)):
    word=(0,1)+tail
    if all((b-a)%q in QR for a,b in zip(word,word[1:])):canonical.append(word)
check(len(canonical)==1729,'square-affine canonical orbit count')
words=set();endpoint=[[0]*q for _ in range(q)];arcs=Counter();cycles=set();closed_words=0
for word in canonical:
    for scale in QR:
        for shift in range(q):
            w=tuple((scale*x+shift)%q for x in word)
            check(w not in words,'affine action free and canonical chart unique')
            words.add(w);endpoint[w[0]][w[-1]]+=1;arcs.update(zip(w,w[1:]))
            if (w[0]-w[-1])%q in QR:
                closed_words+=1;i=w.index(0);cycles.add(w[i:]+w[:i])
check(len(words)==95095,'expanded word total')
check(len(arcs)==55 and set(arcs.values())=={17290},'all arc incidences')
check(len(cycles)==5505 and closed_words==11*len(cycles),'rotation identified cycles')
for a in range(q):
    for b in range(q):
        check(endpoint[a][b]==(0 if a==b else 628 if (b-a)%q in QR else 1101),'full endpoint matrix')
# Scope-amendment control: add m disjoint three-edge circles, retaining every vertex.
closed_component_cases=0
for r in range(1,5):
    P=tuple((2*i,2*i+1) for i in range(r))
    for Q in matchings(tuple(range(2*r))):
        rows=[(1<<a)|(1<<b) for a,b in P+Q]
        bare=2*r-gf2rank(rows)
        for m in range(4):
            augmented=list(rows)
            for j in range(m):
                a=2*r+3*j
                augmented.extend([(1<<a)|(1<<(a+1)),(1<<(a+1))|(1<<(a+2)),(1<<(a+2))|(1<<a)])
            observed=2*r+3*m-gf2rank(augmented)
            check(observed==bare+m,'pre-existing closed components add to matching loops')
            if m:check(observed!=bare,'same endpoint matching cannot see additional circles')
            closed_component_cases+=1

print(json.dumps({'status':'PASS','closed_component_cases':closed_component_cases,'matchings':matching_rows,'support_graphs':{'n':n,'linear_forests':len(LF),'ordered_pairs':len(LF)**2,'distinct_unions':len(unions),'verified_degeneracy_bound':3,'verified_chromatic_bound':4},'fox_full_diagram_colorings':fox,'paley11':{'affine_word_orbits':len(canonical),'H':len(words),'endpoint_row_0':endpoint[0],'all_55_arc_counts':17290,'closed_rooted_words':closed_words,'direct_cycles_mod_rotation':len(cycles)}},indent=2,sort_keys=True))
