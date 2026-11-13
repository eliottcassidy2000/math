"""Literal marked cover controls for the slot phase note."""
from itertools import combinations, permutations, product
from fractions import Fraction

def ck(test, receipt):
    if not test:
        raise AssertionError(receipt)

def arcs(blocks):
    ans=set()
    for i,A in enumerate(blocks):
        for B in blocks[i+1:]:
            ans.update(product(A,B))
        if len(A)==3:
            ans.update(zip(A,A[1:]+A[:1]))
    return ans

def covers(n,k,a):
    out=set()
    for v in permutations(range(n)):
        for c in combinations(range(1,n),k-1):
            t=(0,)+c+(n,)
            p=tuple(v[i:j] for i,j in zip(t,t[1:]))
            if all(e in a for path in p for e in zip(path,path[1:])):
                out.add(tuple(sorted(p)))
    return out

def phase_count(g,D,q,geometric):
    todo=set(D)
    f=1
    while todo:
        p=todo.pop()
        h=tuple(g[x] for x in p)
        L=1
        while h!=p:
            todo.remove(h)
            h=tuple(g[x] for x in h)
            L+=1
        f*=((q**L if geometric else q)-1)**(len(p)-1)
    return f

def stable(g,D):
    return {tuple(g[x] for x in p) for p in D}==set(D)

def run():
    A=((0,),(1,),(2,3,4),(5,))
    B=((2,3,4),(0,),(1,),(5,))
    j={v:i for i,C in enumerate(B) for v in C}
    for k in range(1,7):
        s=covers(6,k,arcs(A))
        t=covers(6,k,arcs(B))
        def reorder(p):
            return tuple(sorted(p,key=j.__getitem__))
        u={tuple(sorted(reorder(p) for p in D)) for D in s}
        ck(u==t,(k,len(s),len(t)))
        ck(all(len({frozenset(p) for p in D})==k for D in s), k)
    print("all six marked occupancy layers agree under SCC reorder")
    D=covers(6,3,arcs(((0,1,2),(3,4,5))))
    G=[tuple((x+a)%3 if x<3 else 3+(x+b-3)%3 for x in range(6))
       for a,b in product(range(3),repeat=2)]
    for q in (2,3,4):
        split,geom=[Fraction(sum(phase_count(g,F,q,tag) for g in G for F in D if stable(g,F)),len(G))
                    for tag in (False,True)]
        ck(geom>split and split.denominator==geom.denominator==1,(q,split,geom))
        print("two triangle twist, q=%s: split=%s geometric=%s" % (q,split,geom))

if __name__=="__main__":
    run()
