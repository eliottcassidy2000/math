"""Exact joint line-event compiler retaining full bipartite incidence.

All ordered pairs of non-axis collinear grid triples for n=4,...,8.
All even-cycle skeleton types at each size; exact shore-label expectations.
No simulation, independence premise, or asymptotic inference.
"""
from collections import Counter, defaultdict
from fractions import Fraction as F
from itertools import combinations
from math import factorial, gcd


def need(ok,msg):
    if not ok:
        raise RuntimeError(msg)


def falling(n,k):
    return factorial(n)//factorial(n-k) if k<=n else 0


def signature(points):
    adj=defaultdict(set)
    for x,y in points:
        adj[(0,x)].add((1,y)); adj[(1,y)].add((0,x))
    if any(len(v)>2 for v in adj.values()):
        return None
    seen=set(); result=[]
    for v in adj:
        if v in seen:
            continue
        stack=[v]; component=set()
        while stack:
            u=stack.pop()
            if u not in component:
                component.add(u); stack.extend(adj[u])
        seen.update(component)
        nr=sum(u[0]==0 for u in component)
        nc=len(component)-nr
        edges=sum(len(adj[u]) for u in component)//2
        cyc=all(len(adj[u])==2 for u in component)
        result.append((edges,nr,nc,cyc))
    return tuple(sorted(result))


def grid_copies(n,sig):
    nr=sum(c[1] for c in sig); nc=sum(c[2] for c in sig)
    aut=1
    for (edges,r,c,cyc),multiplicity in Counter(sig).items():
        one=edges if cyc else (2 if edges%2==0 else 1)
        aut*=one**multiplicity*factorial(multiplicity)
    numerator=falling(n,nr)*falling(n,nc)
    need(numerator%aut==0,'shore automorphism denominator')
    return numerator//aut


def collinear_triples(n):
    # Build maximal nonaxis lines from all grid point pairs.
    lines=defaultdict(set)
    points=list((x,y) for x in range(n) for y in range(n))
    for p,q in combinations(points,2):
        dx,dy=q[0]-p[0],q[1]-p[1]
        if not dx or not dy:
            continue
        d=gcd(dx,dy); a,b=dy//d,-dx//d
        if a<0:
            a,b=-a,-b
        lines[a,b,a*p[0]+b*p[1]].update((p,q))
    return [frozenset(t) for line in lines.values() for t in combinations(sorted(line),3)]


def cycle_skeleton(parts):
    points=[]; start=0
    for r in parts:
        need(r>=2,'simple even cycle')
        points.extend((start+i,start+j) for i in range(r) for j in (i,(i+1)%r))
        start+=r
    return points


def partitions(n,lo=2):
    if n==0:
        yield ()
    for r in range(lo,n+1):
        for tail in partitions(n-r,r):
            yield (r,)+tail


def skeleton_profiles(points):
    out=Counter()
    for k in range(3,7):
        out.update(signature(s) for s in combinations(points,k))
    need(None not in out,'degree two skeleton')
    return out


def main():
    print('Exact full joint-incidence compiler: n=4,...,8, every even-cycle skeleton type')
    print('Cycle parts denote numbers of row vertices, half the graph-cycle lengths.')
    for n in range(4,9):
        triples=collinear_triples(n)
        events=Counter(signature(a|b) for a in triples for b in triples)
        need(None not in events,'union of two rook triples has degree at most two')
        need(sum(events.values())==len(triples)**2,'ordered-pair universe')
        denominator={sig:grid_copies(n,sig) for sig in events}
        need(all(v>0 for v in denominator.values()),'all event types embed in the grid')
        mean=F(len(triples)*4*(2*n-5),n*(n-1)**2*(n-2))
        data={}
        counts={}
        for parts in partitions(n):
            profiles=skeleton_profiles(cycle_skeleton(parts))
            ex2=sum((F(mult*profiles[sig],denominator[sig]) for sig,mult in events.items()),F(0))
            data[parts]=ex2-mean**2
            counts[parts]=profiles
        print('n=',n,'grid_triples=',len(triples),'ordered_pairs=',len(triples)**2,
              'incidence_types=',len(events),'mean=',str(mean))
        for parts,var in sorted(data.items()):
            print('  parts=',parts,'Var[X3]=',str(var))
        if n==4:
            need(data[(2,2)]==F(56,9) and data[(4,)]==F(25,18),'independent complete n4 variance')
        if n==6:
            need(data[(2,2,2)]==F(285968,16875),'n6 disconnected full census control')
            need(data[(2,4)]==F(1313533,101250),'n6 four-cycle full census control')
            need(data[(3,3)]==F(245599,22500),'n6 six-cycle full census control')
            need(data[(6,)]==F(247199,22500),'n6 long-cycle full census control')
        if n>=6:
            base=data[(n,)]
            a=data[tuple(sorted((2,n-2)))]-base
            b=(data[tuple(sorted((3,n-3)))]-base)/(2 if n==6 else 1)
            for parts,var in data.items():
                need(var==base+parts.count(2)*a+parts.count(3)*b,'affine short-cycle variance')
            print('  affine variance coefficients (base,c4,c6)=',tuple(map(str,(base,a,b))))
            # Stronger finite control: all relevant edge-subgraph profiles obey
            # the same affine law, before the geometry-weighted sum.
            k4=tuple(sorted((2,n-2))); k6=tuple(sorted((3,n-3)))
            for sig in events:
                v0=counts[(n,)][sig]
                v4=counts[k4][sig]-v0
                v6=F(counts[k6][sig]-v0,2 if n==6 else 1)
                for parts,profiles in counts.items():
                    need(profiles[sig]==v0+parts.count(2)*v4+parts.count(3)*v6,
                         'affine short-cycle edge-profile law')
    print('PASS. Compiler formulas are exact; the universal short-cycle reduction needs its separate proof.')


if __name__=='__main__':
    main()
