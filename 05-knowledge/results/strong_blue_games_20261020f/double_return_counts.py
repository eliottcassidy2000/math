"""Independent path/subset controls for the two-return-arc family."""
from itertools import combinations, permutations, product
from math import gcd
from pathlib import Path
import hashlib
import json
from boundary_transport_probe import paths, G

HERE=Path(__file__).resolve().parent


def family(a,b):
    n=a+b+1
    t=[sum(1<<v for v in range(u+1,n)) for u in range(n)]
    for u,v in ((0,a),(a,a+b)):
        t[u]^=1<<v;t[v]^=1<<u
    return tuple(t)


def subsets(xs):
    for r in range(len(xs)+1):
        for x in combinations(xs,r):yield frozenset(x)


def construction(a,b):
    A=frozenset(range(a));B=frozenset(range(a+1,a+b+1))
    Is=list(subsets(tuple(range(1,a))))+[A]
    Js=[frozenset()]+[s|{a+b} for s in subsets(tuple(range(a+1,a+b)))]
    return {tuple(sorted(I))+tuple(sorted(J))+(a,)+tuple(sorted(A-I))+tuple(sorted(B-J))
            for I,J in product(Is,Js)}


def automorphisms(t):
    classes={}
    for v,row in enumerate(t):classes.setdefault(row.bit_count(),[]).append(v)
    cs=list(classes.values());answer=0
    for images in product(*(permutations(c) for c in cs)):
        p={v:w for c,im in zip(cs,images) for v,w in zip(c,im)}
        if all(bool(t[u]>>v&1)==bool(t[p[u]]>>p[v]&1) for u in range(len(t)) for v in range(len(t))):answer+=1
    return answer


def main():
    products=[];balanced=[]
    for a,b in product(range(2,6),repeat=2):
        t=family(a,b);ps=set(paths(t));expected=construction(a,b)
        G.need(ps==expected,'path/subset bijection')
        G.need(len(ps)==(2**(a-1)+1)*(2**(b-1)+1),'product count')
        products.append(dict(a=a,b=b,H=len(ps)))
    for k in range(2,8):
        t=family(k,k);ps=set(paths(t));q=2**(k-1)+1
        blue=[p for p in ps if G.blue(t,p)]
        disjoint=[p for p in ps if all(abs(u-v)!=1 for u,v in zip(p,p[1:]))]
        A=range(k);B=range(k+1,2*k+1)
        Q=tuple(v for v in A if v%2)+tuple(v for v in B if not v%2)+(k,)+tuple(v for v in A if not v%2)+tuple(v for v in B if v%2)
        G.need(ps==construction(k,k),'balanced bijection')
        G.need(len(ps)==q*q and len(blue)==q,'square and diagonal count')
        G.need(disjoint==[Q] and Q in blue,'unique all-blue square hostile')
        if k<=4:G.need(automorphisms(t)==1,'score-refined automorphism enumeration')
        balanced.append(dict(k=k,n=2*k+1,H=len(ps),b=len(blue),deficit=len(ps)-len(blue),unique_arc_disjoint_path=Q))
    arithmetic=[]
    for k in range(2,101):
        q=2**(k-1)+1;H=q*q;L=(H-1)//2;C=(H+1)//2
        G.need(q*q+L*L==C*C and gcd(q,L)==1 and C-L==1,'primitive difference-one triple')
        if k<=6:arithmetic.append(dict(k=k,triple=[q,L,C]))
    out=dict(status='FINITE-EXACT controls; infinite scope proved in double_return_counts.md',
             source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
             product_universe='all a,b=2..5, direct DFS versus subset words',products=products,
             balanced_universe='all k=2..7; all directed HPs and positional blue predicates',balanced=balanced,
             arithmetic_control='k=2..100',first_triples=arithmetic)
    (HERE/'double_return_counts.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(dict(products=len(products),balanced=balanced,first_triples=arithmetic),indent=2))


if __name__=='__main__':main()
