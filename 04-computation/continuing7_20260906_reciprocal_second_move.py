"""Exact integer-geometry second-move obstruction. No repository imports."""
from collections import defaultdict
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
from math import gcd
from pathlib import Path
import json
import sys
sys.stdout.reconfigure(newline='\n')

GATES = 0
def check(ok, label):
    global GATES
    GATES += 1
    if not ok:
        raise RuntimeError(label)

def base(p):
    return [1, 0] + [pow(x, -1, p) for x in range(2, p)]

def swap(g, u, v):
    return [v if y == u else u if y == v else y for y in g]

def det(a, b, c):
    return (b[0]-a[0])*(c[1]-a[1])-(c[0]-a[0])*(b[1]-a[1])

def literal(g):
    n = len(g)
    return {(a,b,c) for a,b,c in combinations(range(n),3)
            if (b-a)*(g[c]-g[a]) == (c-a)*(g[b]-g[a])}

def anchor_decoder(g, p):
    """Partition by number of off-conic points; native unoriented directions."""
    on = [(x,y) for x,y in enumerate(g) if x*y % p == 1]
    off = [(x,y) for x,y in enumerate(g) if x*y % p != 1]
    check(len(off) <= 4, 'at most four changed/off-conic anchors')
    answer = set()
    for a in off:
        groups = defaultdict(list)
        for b in on:
            dx,dy = b[0]-a[0], b[1]-a[1]
            d = gcd(dx,dy)
            dx,dy = dx//d,dy//d
            if dx < 0 or (dx == 0 and dy < 0):
                dx,dy = -dx,-dy
            groups[dx,dy].append(b)
        for points in groups.values():
            for b,c in combinations(points,2):
                answer.add(tuple(sorted((a[0],b[0],c[0]))))
    for a,b in combinations(off,2):
        for c in on:
            if det(a,b,c) == 0:
                answer.add(tuple(sorted((a[0],b[0],c[0]))))
    for a,b,c in combinations(off,3):
        if det(a,b,c) == 0:
            answer.add(tuple(sorted((a[0],b[0],c[0]))))
    return answer

def packets(g):
    p = len(g)
    triples = literal(g)
    pt = sorted(t for t in triples if 0 in t)
    check(len(triples) == 2*len(pt), 'transpose twin count')
    check(all((0 in t) ^ (1 in t) for t in triples), 'all old triples one anchor')
    out = []
    for t in pt:
        xs = set(t)-{0}
        ys = {g[x] for x in xs}
        check(not xs & ys, 'one packet shores disjoint')
        check(p-1 not in xs | ys, 'fixed point absent')
        out.append((xs,ys))
    for (x,y),(xx,yy) in combinations(out,2):
        check(not x & xx and not y & yy, 'same-anchor endpoint disjointness')
        check(len(x & yy) <= 1 and len(y & xx) <= 1, 'cross intersections at most one')
        check({g[z] for z in x & yy} == y & xx, 'inverse cross matching')
    return triples,out

def finite(p, all_moves):
    g = base(p)
    old, pack = packets(g)
    supports = [{g[x] for x in t} for t in old]
    rows, hits, safe = [], [], []
    for u,v in combinations(range(p),2):
        hit = all({u,v} & s for s in supports)
        if all_moves or hit:
            h = swap(g,u,v)
            direct = literal(h)
            decoded = anchor_decoder(h,p)
            check(direct == decoded, 'complete native decoder equals literal triangles')
            if not hit:
                check(bool(direct), 'retained old triangle persists')
            rows.append([u,v,len(direct)])
            if hit:
                hits.append([u,v,[list(t) for t in sorted(direct)]])
            if not direct:
                safe.append([u,v])
        else:
            check(any(not ({u,v} & s) for s in supports), 'unhit certificate')
    if p == 37:
        check(len(rows) == 666 and len(hits) == 9 and not safe, 'p37 exhaustive obstruction')
    if p == 43:
        check(len(rows) == 903 and len(hits) == 9 and safe == [[0,13],[1,10]], 'p43 exact positive moves')
    if p in (113,197):
        check(len(pack) == {113:2,197:3}[p] and not safe, 'multiple-packet obstruction')
        check(all(row[:2] == [0,1] for row in hits), 'only undo hits all old triples in control')
    # A move involving the retained fixed point cannot cover an old twin pair.
    check(all(not all({u,p-1} & s for s in supports) for u in range(p-1)), 'fixed-point move cannot repair')
    print('p',p,'old triples',len(old),'native packets',len(pack),'covers',len(hits),'safe',safe)
    return {'p':p,'old':sorted(old),'packets':[[sorted(x),sorted(y)] for x,y in pack],
            'all_moves_literal':all_moves,'move_counts':rows,'covers':hits,'safe':safe,
            'safe_permutations':[swap(g,*uv) for uv in safe]}

# Tiny exact polynomials, ascending coefficients, used for all r >= 0.
def trim(a):
    a = list(map(F,a))
    while len(a)>1 and not a[-1]: a.pop()
    return tuple(a)
def add(a,b):
    return trim([(a[i] if i<len(a) else 0)+(b[i] if i<len(b) else 0) for i in range(max(len(a),len(b)))])
def neg(a): return tuple(-x for x in a)
def sub(a,b): return add(a,neg(b))
def mul(a,b):
    c = [F(0)]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b): c[i+j] += x*y
    return trim(c)
def scale(a,k): return trim([k*x for x in a])
def divide(a,b):
    a = list(trim(a)); b = trim(b)
    q = [F(0)]*max(1,len(a)-len(b)+1)
    while len(a)>=len(b) and any(a):
        i = len(a)-len(b); t = a[-1]/b[-1]; q[i] = t
        for j,x in enumerate(b): a[i+j] -= t*x
        a = list(trim(a))
    return trim(q),trim(a)
def positive(a):
    a = trim(a)
    return len(a)<=2 and a[0]>0 and (len(a)==1 or a[1]>=0)
def unequal(a,b):
    d = sub(a,b)
    if len(d)==1: return bool(d[0])
    r = -d[0]/d[1]
    return r < 0 or r.denominator != 1
def pdet(a,b,c):
    return sub(mul(sub(b[0],a[0]),sub(c[1],a[1])),mul(sub(c[0],a[0]),sub(b[1],a[1])))

def uniform_family():
    p = (F(37),F(360)); one = (F(1),); zero = (F(0),)
    X1,Y1,X2,Y2 = (24,216),(17,160),(30,270),(21,200)
    # Each pair lists native x and y as affine polynomials in r.
    R,S = (X1,Y1),(X2,Y2)
    P,Q = (zero,one),(one,zero)
    order = [zero,one,Y1,Y2,X1,X2,sub(p,one)]
    check(all(positive(sub(b,a)) for a,b in zip(order,order[1:])), 'all-r strict template order')
    check(pdet(P,R,S) == (0,), 'all-r old P triangle')
    check(pdet(Q,(Y1,X1),(Y2,X2)) == (0,), 'all-r old Q triangle')
    witnesses = [
        ((one,X1), ((Y1,one),((19,180),(2,)),((33,320),(9,)))),
        ((one,X2), ((Y2,one),((25,240),(3,)),((31,300),(6,)))),
        ((Y1,X2), ((X1,X2),((23,216),(29,270)),((8,90),(14,144)))),
    ]
    cert = []
    def reciprocal(point):
        quotient,remainder = divide(sub(mul(*point),one),p)
        check(remainder == (0,) and all(q.denominator == 1 for q in quotient), 'all-r native reciprocal product identity')
        for z in point:
            check(positive(sub(z,one)) and positive(sub(p,z)), 'all-r native untouched bounds')
        return quotient
    reciprocal(R); reciprocal(S)
    for labels, points in witnesses:
        u,v = labels
        check(pdet(*points) == (0,), 'all-r created actual triangle')
        check(all(unequal(a[0],b[0]) for a,b in combinations(points,2)), 'all-r distinct triangle rows')
        check(points[0][1] == v or points[0][1] == u, 'moved output recorded')
        old_y = u if points[0][1] == v else v
        reciprocal((points[0][0],old_y))
        qs=[]
        for point in points[1:]:
            qs.append(reciprocal(point))
            check(unequal(point[1],u) and unequal(point[1],v), 'all-r other witness points are unchanged')
        # Corresponding transpose operation is the conjugate g tau g.
        old_inverse = {tuple(one):zero, tuple(X1):Y1, tuple(X2):Y2, tuple(Y1):X1, tuple(Y2):X2}
        trans_labels = (old_inverse[tuple(u)],old_inverse[tuple(v)])
        trans_points = tuple((y,x) for x,y in points)
        check(pdet(*trans_points) == (0,), 'transpose creates companion actual triangle')
        cert.append({'labels':labels,'points':points,'product_quotients':qs,
                     'transpose_labels':trans_labels,'transpose_points':trans_points})
    covers = {frozenset((tuple(u),tuple(v))) for u in (one,Y1,Y2) for v in (zero,X1,X2)}
    covered = {frozenset((tuple(a),tuple(b))) for row in cert for a,b in (row['labels'],row['transpose_labels'])}
    for u,v in ((one,zero),(Y1,X1),(Y2,X2)):
        covered.add(frozenset((tuple(u),tuple(v))))
        check(pdet((u,u),(v,v),(sub(p,one),sub(p,one))) == (0,), 'three diagonal cover witnesses')
    check(covers == covered and len(covered)==9, 'all nine necessary covers obstructed uniformly')
    check(gcd(37,360)==1, 'Dirichlet coprime progression')
    print('all-r progression p=37+360r: six nondiagonal covers blocked by three exact witnesses and transposes; three diagonal covers blocked')
    return cert

def main():
    data = {'finite':[finite(p,p in (37,43)) for p in (37,43,113,197)],
            'uniform_family':uniform_family()}
    def default(v):
        if isinstance(v,F): return str(v)
        raise TypeError(type(v))
    blob = (json.dumps(data,sort_keys=True,indent=2,default=default)+'\n').encode()
    here=Path(__file__).resolve().parent
    dest=here.parent/'05-knowledge/results' if here.name=='04-computation' else here
    path=dest/(Path(__file__).stem+'_certificate.json')
    path.write_bytes(blob)
    print('certificate SHA256',sha256(blob).hexdigest())
    print('PASS',GATES,'always-active exact gates')
    print('Scope: prescribed first swap and one ordinary second output transposition; actual p-point integer board only.')

if __name__ == '__main__': main()
