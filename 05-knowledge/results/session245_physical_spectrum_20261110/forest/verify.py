"""Exact leaf-run arithmetic compiler and independent controls. Python 3 stdlib."""
from itertools import permutations, product
from collections import defaultdict
from fractions import Fraction
import json

def compiler(k):
    if not isinstance(k,int) or not 1<=k<=4: raise ValueError('k must be 1..4')
    states=list(product(*([range(3)]*k+[(0,2)]))); index={s:i for i,s in enumerate(states)}
    tables=[]
    for p in permutations(range(k)):
        letters=[]; operators=[]
        for a in range(-1,k):
            bins=[]
            for g in range(k+1):
                op=[]
                for s in states:
                    if s[g]==0 and g>0 and a==p[g-1]:op.append(-1);continue
                    t=list(s);t[g]=2 if g==k or a==p[g] else 1
                    op.append(index[tuple(t)])
                bins.append(op)
            operators.append(bins)
            letters.append([[op[i] for op in bins if op[i]>=0] for i in range(len(states))])
        accept=[all(s[g]!=1 for g in range(k)) and all(s[g]!=0 or p[g-1]<p[g] for g in range(1,k)) for s in states]
        tables.append((letters,accept,operators))
    return states,tables

def advance(v,t):
    out=[0]*len(v)
    for i,x in enumerate(v):
        if x:
            for j in t[i]:out[j]+=x
    return out

def initial(comp):return tuple(tuple([1]+[0]*(len(comp[0])-1)) for _ in comp[1])
def step(v,a,comp):return tuple(tuple(advance(x,table[0][a+1])) for x,table in zip(v,comp[1]))
def output(v,comp):return sum(sum(x for x,ok in zip(vec,tab[1]) if ok) for vec,tab in zip(v,comp[1]))
def count(word,k=2,comp=None):
    comp=comp or compiler(k);v=initial(comp)
    for a in word:
        if not -1<=a<k:raise ValueError('bad leaf symbol')
        v=step(v,a,comp)
    return output(v,comp)

def graph(word,k):
    n=k+len(word);out=[sum(1<<j for j in range(i+1,n)) for i in range(n)]
    for v,a in enumerate(word,k):
        if a>=0:out[a]^=1<<v;out[v]|=1<<a
    return out

def hp(out):
    n=len(out);d=[[0]*n for _ in range(1<<n)]
    for v in range(n):d[1<<v][v]=1
    for s in range(1,1<<n):
        for v in range(n):
            if s>>v&1:
                prev=s^(1<<v)
                if prev:d[s][v]=sum(d[prev][u] for u in range(n) if prev>>u&1 and out[u]>>v&1)
    return sum(d[-1])

def require(x,msg):
    if not x:raise RuntimeError(msg)

def coeffs(k,pre,post,a):
    comp=compiler(k);vals=[count(pre+(a,)*m+post,k,comp) for m in range(1,k+2)]
    aug=[[Fraction(j**m) for j in range(1,k+2)]+[Fraction(vals[m-1])] for m in range(1,k+2)]
    for i in range(k+1):
        pivot=aug[i][i];aug[i]=[x/pivot for x in aug[i]]
        for r in range(k+1):
            if r!=i:
                q=aug[r][i];aug[r]=[x-q*y for x,y in zip(aug[r],aug[i])]
    return [row[-1] for row in aug]

def main():
    checked={}; projector_checks=0
    for k,maxleaf in [(1,6),(2,6),(3,4)]:
        comp=compiler(k);ncheck=0
        for m in range(maxleaf+1):
            for w in product(range(-1,k),repeat=m):
                require(count(w,k,comp)==hp(graph(w,k)),('DP mismatch',k,w));ncheck+=1
        checked[str(k)]=ncheck
        for _,_,ops in comp[1]:
            for bins in ops:
                for i in range(len(comp[0])):
                    for x in bins:
                        xi=x[i];require(xi<0 or x[xi]==xi,'idempotence');projector_checks+=1
                        for y in bins:
                            yi=y[i];xy=-1 if yi<0 else x[yi];yx=-1 if xi<0 else y[xi]
                            require(xy==yx,'commutation');projector_checks+=1
    comp=compiler(2)
    for a in range(10):
        require(count((-1,)*a+(0,1),2,comp)==2*3**a+2**(a+2)-1,'two-star formula')
    families={}
    from math import factorial
    for k in range(1,5):
        c=coeffs(k,(),tuple(range(k)),-1);cp=compiler(k)
        require(c[-1]==factorial(k),'leading coefficient')
        for a in range(4):
            w=(-1,)*a+tuple(range(k))
            require(count(w,k,cp)==hp(graph(w,k)),'independent family DP')
        for m in range(1,9):require(count((-1,)*m+tuple(range(k)),k,cp)==sum(c[j-1]*j**m for j in range(1,k+2)),'run formula')
        families[str(k)]={'positive_run_coefficients':[str(x) for x in c],'zero_run':count(tuple(range(k)),k,cp)}
    # Hostile controls: a scalar does not determine its continuation menu.
    bycount={};hostile=None
    for m in range(2,6):
        bycount={}
        for w in product(range(-1,2),repeat=m):
            h=count(w,2,comp);menu=tuple(count(w+(a,),2,comp) for a in range(-1,2))
            if h>1 and h in bycount and bycount[h][1]!=menu:
                hostile={'H':h,'words':[bycount[h][0],w],'ordered_child_counts':[bycount[h][1],menu]};break
            bycount[h]=(w,menu)
        if hostile:break
    require(hostile is not None,'missing hostile')
    witness173=(-1,0,0,0,-1,1)
    require(count(witness173,2,comp)==173 and hp(graph(witness173,2))==173,'173 witness')
    # Complete word census by length (fixed initial roots; not full forests).
    images=[set() for _ in range(11)];totals=[0]*11
    def visit(v,depth):
        images[depth].add(output(v,comp));totals[depth]+=1
        if depth<10:
            for a in range(-1,2):visit(step(v,a,comp),depth+1)
    visit(initial(comp),0)
    census=[]
    for m,im in enumerate(images):
        require(totals[m]==3**m,'word count')
        miss=next(h for h in range(1,max(im)+100,2) if h not in im and h not in(7,21))
        census.append({'leaves':m,'orders':m+2,'words':totals[m],'distinct':len(im),'maximum':max(im),'first_allowed_miss':miss,'contains_173':173 in im})
    print(json.dumps({'status':'PROVED compiler and fixed-run spectral identity; FINITE-EXACT controls; all-target coverage OPEN','independent_DP_words':checked,'projector_checks':projector_checks,'two_star_family':'H(a)=2*3^a+2^(a+2)-1, a>=0','k_star_families':families,'scalar_hostile':hostile,'witness173':{'word':witness173,'outgoing_masks':graph(witness173,2)},'two_root_census':census},indent=2,sort_keys=True))
if __name__=='__main__':main()
