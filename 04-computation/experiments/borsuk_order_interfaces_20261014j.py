"""Exact finite Borsuk direction/order interfaces; no unrestricted Borsuk oracle."""
from dataclasses import dataclass, replace
from fractions import Fraction as F
from itertools import combinations, permutations, product
from functools import lru_cache
from math import factorial
import json


def integer(x, lower=0):
    if type(x) is not int or x < lower:
        raise ValueError('exact integer outside domain')
    return x


def rational(x):
    if type(x) not in (int, F):
        raise ValueError('exact rational required')
    return F(x)


def vector(x, dim=None):
    if type(x) is not tuple or (dim is not None and len(x) != dim):
        raise ValueError('exact vector shape')
    return tuple(rational(a) for a in x)


def dot(x, y):
    return sum((a*b for a,b in zip(x,y)),F(0))


def sub(x,y):
    return tuple(a-b for a,b in zip(x,y))


def solve_unique(rows, rhs):
    """Exact possibly rectangular system; return None unless unique consistent."""
    if not rows:
        return () if not rhs else None
    n=len(rows[0]); a=[list(map(F,row))+[F(b)] for row,b in zip(rows,rhs)]
    pivot=[]; k=0
    for j in range(n):
        z=next((i for i in range(k,len(a)) if a[i][j]),None)
        if z is None: continue
        a[k],a[z]=a[z],a[k]
        c=a[k][j];a[k]=[x/c for x in a[k]]
        for i in range(len(a)):
            if i!=k and a[i][j]:
                c=a[i][j];a[i]=[x-c*y for x,y in zip(a[i],a[k])]
        pivot.append(j);k+=1
    if any(not any(row[:n]) and row[n] for row in a) or k!=n:
        return None
    out=[F(0)]*n
    for i,j in enumerate(pivot):out[j]=a[i][n]
    return tuple(out)


def positive_definite(g):
    """Exact symmetric LDL pivots, no floating eigenvalues."""
    d=len(g)
    if any(g[i][j]!=g[j][i] for i in range(d) for j in range(d)):return False
    a=[list(row) for row in g]
    for i in range(d):
        if a[i][i]<=0:return False
        for j in range(i+1,d):
            for k in range(i+1,d):
                a[j][k]-=a[j][i]*a[i][k]/a[i][i]
    return True


@dataclass(frozen=True)
class Scene:
    points: tuple
    metric: tuple
    threshold2: F


def scene(points, metric, threshold2):
    if type(metric) is not tuple or not 1<=len(metric)<=4:
        raise ValueError('dimension one through four')
    d=len(metric);g=tuple(vector(row,d) for row in metric)
    if not positive_definite(g):raise ValueError('positive definite metric required')
    if type(points) is not tuple or not 2<=len(points)<=8:
        raise ValueError('two through eight marked points')
    p=tuple(vector(x,d) for x in points)
    if len(set(p))!=len(p):raise ValueError('distinct labelled points required')
    t=rational(threshold2)
    if t<=0:raise ValueError('positive threshold required')
    return Scene(p,g,t)


def audit_scene(s):
    if type(s) is not Scene or s!=scene(s.points,s.metric,s.threshold2):
        raise ValueError('typed scene required')
    return s


def distance2(s,i,j):
    v=sub(s.points[i],s.points[j])
    return sum(v[a]*s.metric[a][b]*v[b] for a in range(len(v)) for b in range(len(v)))


def edges(s):
    audit_scene(s)
    return tuple((i,j) for i,j in combinations(range(len(s.points)),2)
                 if distance2(s,i,j)>=s.threshold2)


def diameter2(s):
    audit_scene(s)
    return max(distance2(s,i,j) for i,j in combinations(range(len(s.points)),2))


def marked_differences(s, selected):
    audit_scene(s);integer(selected)
    if selected>>len(s.points):raise ValueError('subset outside source')
    out=[]
    for i,j in edges(s):
        if selected>>i&1:out.append(((i,j),sub(s.points[i],s.points[j])))
        if selected>>j&1:out.append(((j,i),sub(s.points[j],s.points[i])))
    return tuple(out)


@dataclass(frozen=True)
class Separation:
    source: Scene
    selected: int
    kind: str
    direction: tuple
    labels: tuple
    weights: tuple


def separate(s, selected):
    """Direction with every gap >=1, or <=d+1 labelled convex obstruction."""
    marked=marked_differences(s,selected);d=len(s.metric)
    if not marked:return Separation(s,selected,'direction',(F(0),)*d,(),())
    ds=tuple(v for _,v in marked)
    # Projection of 0 onto {w: d_i.w>=1}; active normal basis has <=d members.
    for k in range(1,min(d,len(ds))+1):
        for ids in combinations(range(len(ds)),k):
            vs=tuple(ds[i] for i in ids)
            lam=solve_unique(tuple(tuple(dot(a,b) for b in vs) for a in vs),(1,)*k)
            if lam is None or min(lam)<0:continue
            w=tuple(sum(lam[i]*vs[i][j] for i in range(k)) for j in range(d))
            if all(dot(w,v)>=1 for v in ds):
                return Separation(s,selected,'direction',w,(),())
    # An affinely independent convex dependence has <=d+1 vectors.
    for k in range(1,min(d+1,len(ds))+1):
        for ids in combinations(range(len(ds)),k):
            vs=tuple(ds[i] for i in ids)
            rows=tuple(tuple(v[j] for v in vs) for j in range(d))+((1,)*k,)
            lam=solve_unique(rows,(0,)*d+(1,))
            if lam is not None and min(lam)>0:
                return Separation(s,selected,'blocker',(),tuple(marked[i][0] for i in ids),lam)
    raise ArithmeticError('strict separation alternative failed')


def audit_separation(packet,s,selected):
    audit_scene(s);marked=dict(marked_differences(s,selected));d=len(s.metric)
    if type(packet) is not Separation:raise ValueError('typed separation required')
    audit_scene(packet.source);integer(packet.selected)
    if packet.source!=s or packet.selected!=selected:raise ValueError('wrong marked source')
    if type(packet.kind) is not str:raise ValueError('exact kind required')
    if packet.kind=='direction':
        w=vector(packet.direction,d)
        if packet.labels!=() or packet.weights!=() or any(dot(w,v)<1 for v in marked.values()):
            raise ValueError('false direction certificate')
    elif packet.kind=='blocker':
        if packet.direction!=() or type(packet.labels) is not tuple:
            raise ValueError('false blocker shape')
        for pair in packet.labels:
            if type(pair) is not tuple or len(pair)!=2 or any(type(x) is not int for x in pair):
                raise ValueError('exact oriented source pairs required')
        lam=vector(packet.weights,len(packet.labels))
        if not 1<=len(lam)<=d+1 or len(set(packet.labels))!=len(lam) or min(lam)<=0 or sum(lam)!=1:
            raise ValueError('convex support guard')
        if any(pair not in marked for pair in packet.labels):raise ValueError('unmarked blocker')
        if any(sum(lam[i]*marked[pair][j] for i,pair in enumerate(packet.labels)) for j in range(d)):
            raise ValueError('nonzero blocker sum')
    else:raise ValueError('unknown separation kind')
    return packet


def selected_by(s,w):
    audit_scene(s);w=vector(w,len(s.metric));bad=set()
    for i,j in edges(s):
        z=dot(w,sub(s.points[i],s.points[j]))
        if z<=0:bad.add(i)
        if z>=0:bad.add(j)
    return sum(1<<i for i in range(len(s.points)) if i not in bad)


@dataclass(frozen=True)
class Cover:
    source: Scene
    directions: tuple
    colours: tuple


def compile_cover(s,directions,k=5):
    audit_scene(s);integer(k,1)
    if type(directions) is not tuple or not 1<=len(directions)<=k:raise ValueError('direction budget')
    ws=tuple(vector(w,len(s.metric)) for w in directions)
    masks=tuple(selected_by(s,w) for w in ws)
    if any(not any(m>>i&1 for m in masks) for i in range(len(s.points))):raise ValueError('uncovered vertex')
    c=tuple(next(j for j,m in enumerate(masks) if m>>i&1) for i in range(len(s.points)))
    if any(c[i]==c[j] for i,j in edges(s)):raise ArithmeticError('selection independence failed')
    return Cover(s,ws,c)


def audit_cover(packet,s,k=5):
    if type(packet) is not Cover or type(packet.colours) is not tuple or any(type(x) is not int for x in packet.colours):
        raise ValueError('exact cover receipt required')
    audit_scene(packet.source)
    if packet.source!=s or packet!=compile_cover(s,packet.directions,k):raise ValueError('source/cover mismatch')
    return packet


def colouring(n,es,k):
    integer(n);integer(k,1);adj=[set() for _ in range(n)]
    for a,b in es:adj[a].add(b);adj[b].add(a)
    c=[-1]*n
    def rec(left):
        if not left:return tuple(c)
        v=max(left,key=lambda x:(len({c[y] for y in adj[x] if c[y]>=0}),len(adj[x]),-x))
        used={c[y] for y in adj[v] if c[y]>=0}
        for a in range(min(k,max(c,default=-1)+2)):
            if a not in used:
                c[v]=a;r=rec(left-{v})
                if r is not None:return r
        c[v]=-1
        return None
    return rec(set(range(n)))


@dataclass(frozen=True)
class Partition:
    source: Scene
    colours: tuple


def diameter_partition(s):
    """Total <=8-point R4 theorem, on the exact diameter graph only."""
    audit_scene(s)
    if s.threshold2!=diameter2(s):raise ValueError('exact diameter threshold required')
    c=colouring(len(s.points),edges(s),5)
    if c is None:raise ArithmeticError('proved eight-point theorem failed')
    return Partition(s,c)


def audit_partition(packet,s):
    audit_scene(s)
    if type(packet) is not Partition or type(packet.colours) is not tuple:
        raise ValueError('typed partition required')
    audit_scene(packet.source)
    if packet.source!=s or s.threshold2!=diameter2(s) or len(packet.colours)!=len(s.points):
        raise ValueError('wrong diameter source')
    for c in packet.colours:
        integer(c)
        if c>=5:raise ValueError('five-colour budget')
    if any(packet.colours[i]==packet.colours[j] for i,j in edges(s)):
        raise ValueError('diameter pair in one part')
    return packet


def profiles(n):
    ts=sorted([(1,'I')]+[(i,'P') for i in range(2,n+1)]+[(i,'C') for i in range(3,n+1)])
    def rec(left,start,chosen):
        if not left:yield chosen;return
        for j in range(start,len(ts)):
            size,kind=ts[j]
            if size<=left:yield from rec(left-size,j,chosen+(ts[j],))
    return tuple(rec(n,0,()))


def profile_graph(p):
    off=0;es=[]
    for size,kind in p:
        es.extend((off+i,off+i+1) for i in range(size-1))
        if kind=='C':es.append((off,off+size-1))
        off+=size
    return off,tuple(sorted(es))


def complement(n,es):
    return tuple(x for x in combinations(range(n),2) if x not in set(es))


def near_clique(n,es):
    e=set(es)
    return next((a for a in combinations(range(n),6) if sum((i,j) in e for i,j in combinations(a,2))>=14),None)


def linear_extensions(n,relations):
    return tuple(p for p in permutations(range(n)) if all(p.index(a)<p.index(b) for a,b in relations))


def order_direction(points,p):
    d=len(points[0]);ds=tuple(sub(points[p[i+1]],points[p[i]]) for i in range(len(p)-1))
    # Reuse strict solver without changing metric/edge semantics: explicit active-set solver.
    for k in range(1,min(d,len(ds))+1):
        for ids in combinations(range(len(ds)),k):
            vs=tuple(ds[i] for i in ids)
            lam=solve_unique(tuple(tuple(dot(a,b) for b in vs) for a in vs),(1,)*k)
            if lam is None or min(lam)<0:continue
            w=tuple(sum(lam[i]*vs[i][j] for i in range(k)) for j in range(d))
            if all(dot(w,v)>=1 for v in ds):return w
    return None


def main():
    checks=0
    def check(z):
        nonlocal checks
        checks+=1
        if not z:raise AssertionError('exact control failed')
    def rejects(f):
        try:f()
        except ValueError:check(True)
        else:check(False)
    eye=lambda d:tuple(tuple(F(i==j) for j in range(d)) for i in range(d))
    simplex=scene(((0,0,0,0),(1,0,0,0),(0,1,0,0),(0,0,1,0),(0,0,0,1)),
                  tuple(tuple(2 if i==j else 1 for j in range(4)) for i in range(4)),2)
    u=((2,-1,-1),(-1,2,-1),(-1,-1,2))
    hexagon=scene(u+tuple(tuple(-a for a in p) for p in u),eye(3),24)
    cross=scene(tuple(tuple(s if j==i else 0 for j in range(4)) for i in range(4) for s in (-1,1)),eye(4),4)
    examples=[]
    for name,s in (('regular_simplex5',simplex),('alternating_hexagon',hexagon),('cross_polytope8',cross)):
        good=[];blockers=0
        for mask in range(1<<len(s.points)):
            r=separate(s,mask);check(audit_separation(r,s,mask)==r)
            if r.kind=='direction':
                check(selected_by(s,r.direction)&mask==mask);good.append(r)
            else:blockers+=1
        masks={selected_by(s,r.direction):r.direction for r in good}
        best=None
        for k in range(1,6):
            for chosen in combinations(masks,k):
                union=0
                for a in chosen:union|=a
                if union==(1<<len(s.points))-1:
                    best=compile_cover(s,tuple(masks[a] for a in chosen));break
            if best:break
        check(best is not None and audit_cover(best,s)==best)
        check(audit_partition(diameter_partition(s),s).source==s)
        check(all(distance2(s,i,j)<diameter2(s) for i,j in combinations(range(len(s.points)),2)
                  if best.colours[i]==best.colours[j]))
        examples.append(dict(name=name,vertices=len(s.points),edges=len(edges(s)),
                             selectable_subsets=len(good),blockers=blockers,
                             minimum_direction_cover=len(best.directions),colours=best.colours,
                             directions=best.directions))
    # Independent-set != one-direction class; exact three-vector zero dependence.
    h=separate(hexagon,7)
    check(h.kind=='blocker' and len(h.weights)==3 and h.weights==(F(1,3),)*3)
    check(not any((7>>i&1) and (7>>j&1) for i,j in edges(hexagon)))
    # All120 regular-simplex orders really arise as projections.
    check(len(linear_extensions(5,()))==120)
    for p in permutations(range(5)):
        w=order_direction(simplex.points,p);check(w is not None)
        check(all(dot(w,sub(simplex.points[p[j+1]],simplex.points[p[j]]))>=1 for j in range(4)))
    # Three-point poset splits 1/3--2/3 but its larger branch is unrealizable.
    line=((F(0),),(F(1),),(F(2),));ext=linear_extensions(3,((0,1),))
    small=tuple(p for p in ext if p.index(1)<p.index(2));large=tuple(p for p in ext if p not in small)
    check(len(ext)==3 and len(small)==1 and len(large)==2)
    check(all(order_direction(line,p) is None for p in large))
    check(order_direction(line,small[0]) is not None)
    # Exact all-component profiles for complement maximum degree2.
    profile_rows=[]
    for n in range(1,9):
        candidates=[]
        for p in profiles(n):
            _,h_edges=profile_graph(p);g=complement(n,h_edges)
            saving=sum(s-1 if k=='C' and s==3 else s//2 for s,k in p)
            chi=n-saving
            check(colouring(n,g,chi) is not None)
            if chi>1:check(colouring(n,g,chi-1) is None)
            alpha=sum(1 if k=='I' else ((s+1)//2 if k=='P' else s//2) for s,k in p)
            clique=max((len(a) for j in range(1,n+1) for a in combinations(range(n),j)
                        if all((i,k) in g for i,k in combinations(a,2))),default=0)
            check(alpha==clique)
            if chi>=6 and alpha<=5:
                candidates.append(p);check(near_clique(n,g) is not None)
        profile_rows.append(dict(n=n,profiles=len(profiles(n)),non5_K6free=candidates))
    check(profile_rows[-1]['non5_K6free']==[((1,'I'),(1,'I'),(1,'I'),(5,'C'))])
    # Bipyramid exact metric: four tetrahedron anchors, two opposite tips.
    tetra=((1,1,1,0),(1,-1,-1,0),(-1,1,-1,0),(-1,-1,1,0))
    bipyramid=scene(tetra+((0,0,0,1),(0,0,0,-1)),
                    ((1,0,0,0),(0,1,0,0),(0,0,1,0),(0,0,0,5)),8)
    check(all(distance2(bipyramid,i,j)==8 for i,j in combinations(range(6),2) if (i,j)!=(4,5)))
    check(distance2(bipyramid,4,5)==20==F(5,2)*8)
    check(len(edges(bipyramid))==15 and colouring(6,edges(bipyramid),5) is None)
    rejects(lambda:diameter_partition(bipyramid))
    # It's a unit/threshold graph carrier, but its diameter graph has only the tip edge.
    actual=scene(bipyramid.points,bipyramid.metric,20)
    check(edges(actual)==((4,5),))
    check(audit_partition(diameter_partition(actual),actual).source==actual)
    rejects(lambda:audit_partition(replace(diameter_partition(actual),colours=(True,)*6),actual))
    rejects(lambda:audit_partition(diameter_partition(actual),bipyramid))
    rejects(lambda:audit_separation(replace(h,selected=True),hexagon,7))
    rejects(lambda:audit_separation(replace(h,weights=(F(1,2),)*3),hexagon,7))
    rejects(lambda:audit_separation(replace(h,labels=((False,3),)+h.labels[1:]),hexagon,7))
    rejects(lambda:audit_separation(h,cross,7))
    rejects(lambda:scene(((0,),(0,)),((1,),),1))
    rejects(lambda:scene(((0,),(1,)),((0,),),1))
    rejects(lambda:scene(((0,),(1,)),((1,),),True))
    rejects(lambda:compile_cover(simplex,((1,1,1,1),)))
    receipt=compile_cover(hexagon,((1,2,4),(-1,-2,-4)))
    rejects(lambda:audit_cover(replace(receipt,colours=(True,)+receipt.colours[1:]),hexagon))
    print(json.dumps(dict(checks=checks,status='PROVED finite interfaces and <=8-point R4 case; Borsuk4 OPEN',
        examples=examples,hexagon_blocker=dict(labels=h.labels,weights=h.weights),
        order_split=dict(extensions=ext,small=small,geometrically_empty_large=large),
        degree_two_complement_profiles=profile_rows,
        bipyramid=dict(shared_distance2=8,tip_distance2=20,ratio=F(5,2)),
        limitations='No unrestricted five-direction cover, no dimension-four asymptotic substitution, no balanced-pair conjecture proof'),
        indent=2,default=str))


if __name__=='__main__':main()
