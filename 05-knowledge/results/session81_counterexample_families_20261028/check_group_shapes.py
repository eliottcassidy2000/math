#!/usr/bin/env python3
"""Standard-library checks; literal permutations, no n! search or asserts."""
from itertools import combinations
from collections import deque
import json

def require(x, msg):
    if not x: raise RuntimeError(msg)

def pair(h):
    d,e=(0,),(1,)
    for _ in range(h): d,e=d+d+e,e+e+d
    return d,e

def witness(h,r):
    d,_=pair(r-1)
    s=d*3
    for j in range(r,h):
        d,e=pair(j); s=s+d+e
    return s

def typed(word):
    if len(word)==1: return str(word[0]),0
    a=len(word)//3
    x=[typed(word[i*a:(i+1)*a]) for i in range(3)]
    ts=[t for t,_ in x]
    rots=[tuple(ts[i:]+ts[:i]) for i in range(3)]
    return '('+','.join(min(rots))+')',sum(e for _,e in x)+(len(set(ts))==1)

def arc(word,u,v):
    bu,qu=divmod(u,4); bv,qv=divmod(v,4)
    if bu==bv:
        forward=(qu==0) if 0 in (qu,qv) else ((qv-qu)%3==1)
        return forward ^ bool(word[bu])
    size=len(word)
    while True:
        size//=3
        x,y=bu//size,bv//size
        if x!=y: return (y-x)%3==1
        bu%=size; bv%=size

def triples(word):
    return {t for t in combinations(range(4*len(word)),3)
            if arc(word,t[0],t[1])==arc(word,t[1],t[2])==arc(word,t[2],t[0])}

def compose(p,q): return bytes(p[x] for x in q)
def inverse(p):
    q=[0]*len(p)
    for i,x in enumerate(p): q[x]=i
    return bytes(q)

def generators(h,r):
    L=3**h; n=4*L; a=3**(r-1)
    gs=[]
    for k in range(L):
        p=list(range(n))
        for q in (1,2,3): p[4*k+q]=4*k+q%3+1
        gs.append(bytes(p))
    p=list(range(n))
    for i in range(3):
        for k in range(a):
            for q in range(4): p[4*(i*a+k)+q]=4*(((i+1)%3)*a+k)+q
    gs.append(bytes(p))
    return gs

def closure(gs,n):
    one=bytes(range(n)); seen={one}; todo=deque([one])
    while todo:
        p=todo.popleft()
        for g in gs:
            q=compose(p,g)
            if q not in seen: seen.add(q); todo.append(q)
    return seen

def enum_case(h,r):
    L=3**h; a=3**(r-1); n=4*L
    w=witness(h,r); gs=generators(h,r); one=bytes(range(n))
    for g in gs:
        require(all(arc(w,u,v)==arc(w,g[u],g[v]) for u,v in combinations(range(n),2)), 'generator not automorphism')
    G=closure(gs,n)
    require(len(G)==3**(L+1),'group order')
    center=[p for p in G if all(compose(p,g)==compose(g,p) for g in gs)]
    require(len(center)==3**(L-2*a),'center size')
    # All generator commutators; closure plus normality certifies derived group.
    cs=[compose(compose(compose(inverse(g),inverse(k)),g),k) for g in gs for k in gs]
    D=closure(list(set(cs)-{one}),n)
    require(len(D)==3**(2*a),'derived size')
    require(all(compose(compose(inverse(g),d),g) in D for g in gs for d in D),'derived normality')
    require(all(compose(compose(d,d),d)==one for d in D),'derived exponent')
    c3s={compose(compose(compose(inverse(d),inverse(g)),d),g) for d in D for g in gs}
    C3=closure(list(c3s-{one}),n)
    require(len(C3)==3**a and C3.issubset(set(center)), 'class three and gamma3')
    orders={}
    for p in G:
        p3=compose(compose(p,p),p)
        order=1 if p==one else 3 if p3==one else 9
        require(compose(compose(p3,p3),p3)==one,'exponent > 9')
        orders[order]=orders.get(order,0)+1
    require(9 in orders,'exponent below 9')
    require(orders.get(1,0)+orders.get(3,0)==3**L+2*3**(L-a),'cube-kernel count')
    require(all(compose(compose(p,p),p)==one for p in center),'center exponent')
    return dict(h=h,r=r,n=n,order=len(G),center=len(center),derived=len(D),gamma3=len(C3),nilpotency_class=3,element_orders=orders)

def main():
    checked=[]
    for h in range(1,8):
        types=[]
        for r in range(1,h+1):
            w=witness(h,r); t,e=typed(w)
            require(len(w)==3**h and e==1,'unique active upper node')
            types.append(t)
        require(len(set(types))==h,'witness isotypes distinct')
        for j in range(h+1):
            d,e=pair(j)
            require(typed(d)[0]!=typed(e)[0] and typed(d)[1]==typed(e)[1]==0,'minimal pair')
        checked.append(h)
    for h in (1,2,3):
        baseline=triples(witness(h,1))
        for r in range(1,h+1): require(triples(witness(h,r))==baseline,'literal C3 fiber')
    result={'unique_node_and_type_checks_heights':checked,'literal_C3_checks_heights':[1,2,3],
            'full_permutation_enumerations':[enum_case(h,r) for h in (1,2) for r in range(1,h+1)]}
    print(json.dumps(result,sort_keys=True,indent=2))
if __name__=='__main__': main()
