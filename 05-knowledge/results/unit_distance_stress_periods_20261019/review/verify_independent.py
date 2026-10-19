#!/usr/bin/env python3
"""Independent exact controls. Standard library only; no producer imports/files.

Run from any directory. Writes only independent-results.json beside this file.
All ranks use modular integer elimination, never floating point.
"""
import json
from itertools import product
from pathlib import Path


def rank(rows, p):
    basis = {}
    for row in rows:
        row = {i: a % p for i, a in enumerate(row) if a % p}
        while row:
            j = min(row)
            if j not in basis:
                inv = pow(row[j], -1, p)
                basis[j] = {i: a * inv % p for i, a in row.items()}
                break
            c = row[j]
            for i, a in basis[j].items():
                value = (row.get(i, 0) - c * a) % p
                if value:
                    row[i] = value
                else:
                    row.pop(i, None)
    return len(basis)


def vertices(p):
    return list(product(range(p), repeat=2))


def add(v, d, p):
    return ((v[0] + d[0]) % p, (v[1] + d[1]) % p)


def edge(a, b):
    return tuple(sorted((a, b)))


def directions(p, norm):
    return [d for d in vertices(p) if norm(d) % p == 1
            and d < ((-d[0]) % p, (-d[1]) % p)]


def lines(p, d):
    remaining = set(vertices(p))
    while remaining:
        v = min(remaining)
        line = [add(v, (t * d[0], t * d[1]), p) for t in range(p)]
        remaining.difference_update(line)
        yield [edge(line[t], line[(t + 1) % p]) for t in range(p)]


def reduced_rows(p, ds):
    for a, b in ds:
        for x, y in vertices(p):
            row = [0] * (2 * p)
            row[(y + b) % p] += a
            row[y] -= a
            row[p + (x + a) % p] += b
            row[p + x] -= b
            yield row


def check_rigidity(p):
    D = next(d for d in range(2, p) if pow(d, (p - 1) // 2, p) == p - 1)
    norm = lambda z: z[0] ** 2 - D * z[1] ** 2
    ds = directions(p, norm)
    assert len(ds) == (p + 1) // 2
    u = (1, 0)
    v = next(d for d in ds if d[1])
    # Coordinates relative to independent unit vectors u,v.
    cs = [((d[0] - d[1] * pow(v[1], -1, p) * v[0]) % p,
           d[1] * pow(v[1], -1, p) % p) for d in ds]
    offset_rank = rank(reduced_rows(p, cs), p)
    expected = 0 if p == 3 else 2 * p - 3
    assert offset_rank == expected
    full_rank = 2 * p * (p - 1) + offset_rank
    assert full_rank == (12 if p == 3 else 2 * p * p - 3)
    # Separate direct full-matrix checks, without using the reduction.
    direct = None
    if p <= 7:
        V = vertices(p)
        ix = {v: i for i, v in enumerate(V)}
        rows = []
        for d in ds:
            for a in V:
                b = add(a, d, p)
                row = [0] * (2 * len(V))
                gx, gy = 2 * d[0], -2 * D * d[1]
                row[2 * ix[a]] -= gx
                row[2 * ix[a] + 1] -= gy
                row[2 * ix[b]] += gx
                row[2 * ix[b] + 1] += gy
                rows.append(row)
        direct = rank(rows, p)
        assert direct == full_rank
    E = p * p * len(ds)
    return dict(p=p, edges=E, offset_rank=offset_rank, full_rank=full_rank,
                direct_full_matrix_rank=direct, stress_dimension=E-full_rank,
                extra_after_cycles=E-full_rank-p*len(ds))


def check_fixed_error(p):
    D = 2
    V = vertices(p)
    norm = lambda z: z[0] ** 2 - D * z[1] ** 2
    # Frozen q5 positive control copied explicitly from the earlier checkpoint.
    hs = [(3,1),(0,4),(2,3),(4,2),(1,1),(2,2),(0,3),(2,2),(4,1),
          (1,4),(2,0),(0,2),(1,3),(4,4),(1,3),(2,4),(0,0),(1,1),
          (3,0),(1,1),(2,2),(4,1),(1,0),(3,4),(1,0)]
    X = ({(x,y):(x,7*y) for x,y in V} if p == 3 else
         {v:(v[0]+p*h[0],v[1]+p*h[1]) for v,h in zip(V,hs)})
    errors = {}
    for d in directions(p, norm):
        for a in V:
            b = add(a, d, p)
            delta = (X[b][0]-X[a][0], X[b][1]-X[a][1])
            residue_delta = (b[0]-a[0], b[1]-a[1])
            assert (norm(delta)-1) % p == 0
            z = ((norm(delta)-1)//p) % p
            assert z == ((norm((-delta[0],-delta[1]))-1)//p) % p
            h = ((X[b][0]-b[0]-X[a][0]+a[0])//p,
                 (X[b][1]-b[1]-X[a][1]+a[1])//p)
            c = ((norm(residue_delta)-1)//p) % p
            R = (2*d[0]*h[0]-2*D*d[1]*h[1]) % p
            assert z == (c+R) % p
            errors[edge(a,b)] = z
    C = [C for d in directions(p,norm) for C in lines(p,d)]
    for cycle in C:
        assert sum(errors[e] for e in cycle) % p == (-2) % p
        assert [errors[e] for e in cycle if errors[e]] == [(-2) % p]
    assert len(errors) == p*len(C)
    return dict(p=p, edges=len(errors), failures=sum(bool(x) for x in errors.values()),
                cycles=len(C), fixed_error=(-2)%p, affine_and_orientation_checked=True)


def triangular(p, full=False):
    V = vertices(p)
    Q = lambda z: z[0]**2+z[0]*z[1]+z[1]**2
    mul = lambda a,b: ((a[0]*b[0]-a[1]*b[1])%p,
                       (a[0]*b[1]+a[1]*b[0]+a[1]*b[1])%p)
    hexagon = [(1,0),(0,1),(p-1,1),(p-1,0),(0,p-1),(1,p-1)]
    base = set()
    for x,y in V:
        for d,choose in [((1,0),x==0),((0,1),y==0),((-1,1),x==0)]:
            if choose:
                base.add(edge((x,y),add((x,y),d,p)))
    reps = [(1,0)]
    if full:
        remaining = {d for d in V if Q(d)%p == 1}
        reps = []
        while remaining:
            u = min(remaining)
            reps.append(u)
            remaining.difference_update(mul(u,d) for d in hexagon)
    selected = {edge(mul(u,a),mul(u,b)) for u in reps for a,b in base}
    families = [[mul(u,d) for d in hexagon] for u in reps]
    ds = [d for H in families for d in H[:3]]
    all_edges = {edge(v,add(v,d,p)) for d in ds for v in V}
    C = [C for d in ds for C in lines(p,d)]
    assert len(all_edges) == 3*p*p*len(reps)
    assert len(selected) == len(C) == 3*p*len(reps)
    assert all(len(set(cycle)&selected)==1 for cycle in C)
    def carry(e):
        a,b=e
        delta=(b[0]-a[0],b[1]-a[1])
        assert (Q(delta)-1)%p == 0
        return ((Q(delta)-1)//p)%p
    carries = {e:carry(e) for e in all_edges}
    for cycle in C:
        assert sum(carries[e] for e in cycle)%p == (-2)%p
    for H in families:
        for v in V:
            ring=[add(v,d,p) for d in H]
            rim=[edge(ring[i],ring[(i+1)%6]) for i in range(6)]
            spokes=[edge(v,r) for r in ring]
            assert sum(e in selected for e in rim)==sum(e in selected for e in spokes)
            assert (sum(carries[e] for e in rim)-sum(carries[e] for e in spokes))%p==0
    # Independently integrate the actual full support on the two axis families.
    def rhs(a,d):
        e=edge(a,add(a,d,p))
        return (-2*(e in selected)-carries[e])%p
    A,B={},{}
    for y in range(p):
        A[0,y]=0
        for x in range(p-1):
            A[x+1,y]=(A[x,y]+rhs((x,y),(1,0)))%p
        assert (A[0,y]-A[p-1,y]-rhs((p-1,y),(1,0)))%p==0
    for x in range(p):
        B[x,0]=0
        for y in range(p-1):
            B[x,y+1]=(B[x,y]+rhs((x,y),(0,1)))%p
        assert (B[x,0]-B[x,p-1]-rhs((x,p-1),(0,1)))%p==0
    rows, augmented = [],[]
    for d in ds:
        for (x,y),row in zip(V,reduced_rows(p,[d])):
            w=add((x,y),d,p)
            residual=(rhs((x,y),d)-d[0]*(A[w]-A[x,y])-d[1]*(B[w]-B[x,y]))%p
            rows.append(row)
            augmented.append(row+[residual])
    rr, ar=rank(rows,p),rank(augmented,p)
    assert rr==2*p-3 and ar==rr+1
    periods=None
    if not full:
        T={}
        for x,y in V:
            w=add((x,y),(-1,1),p)
            residual=(rhs((x,y),(-1,1))+A[w]-A[x,y]-B[w]+B[x,y])%p
            assert residual==(-3*(x==0)+2*(x==1)+2*(y==0)-(y==p-1))%p
            T[x,y]=residual
        xp={sum(T[x,y] for x in range(p))%p for y in range(p)}
        yp={sum(T[x,y] for y in range(p))%p for x in range(p)}
        assert xp=={p-1} and yp=={1}
        periods=[p-1,1]
    return dict(p=p,full_graph=full,vertices=p*p,edges=len(all_edges),
                triangular_families=len(reps),selected_edges=len(selected),cycles=len(C),
                wheel_support_and_carry_checks=p*p*len(reps),
                reduced_rank=rr,augmented_rank=ar,periods=periods,
                representatives=reps)


def check_wheel_span(p):
    V=vertices(p)
    ds=[(1,0),(0,1),(-1,1)]
    H=[(1,0),(0,1),(-1,1),(-1,0),(0,-1),(1,-1)]
    E=sorted({edge(v,add(v,d,p)) for v in V for d in ds})
    ix={e:i for i,e in enumerate(E)}
    W=[]
    for v in V:
        ring=[add(v,d,p) for d in H]
        row=[0]*len(E)
        for i in range(6):
            row[ix[edge(ring[i],ring[(i+1)%6])]]+=1
            row[ix[edge(v,ring[i])]]-=1
        # Independently check that this is an actual rigidity stress.
        equilibrium={z:[0,0] for z in V}
        for i,a in enumerate(row):
            if not a:
                continue
            l,r=E[i]
            dx,dy=r[0]-l[0],r[1]-l[1]
            gradient=(2*dx+dy,dx+2*dy)
            for j in (0,1):
                equilibrium[l][j]-=a*gradient[j]
                equilibrium[r][j]+=a*gradient[j]
        assert all(a%p==0 for z in equilibrium.values() for a in z)
        W.append(row)
    for f in [lambda v:1,lambda v:v[0],lambda v:v[1]]:
        assert all(sum(f(v)*W[j][i] for j,v in enumerate(V))%p==0
                   for i in range(len(E)))
    C=[]
    for d in ds:
        for cycle in lines(p,d):
            row=[0]*len(E)
            for e in cycle:
                row[ix[e]]=1
            C.append(row)
    wr=rank(W,p)
    cr=rank(C,p)
    joined=rank(W+C,p)
    assert wr==p*p-3 and cr==3*p and joined==p*p+2
    return dict(p=p,scope='three-direction triangular graph only',
                edges=3*p*p,full_stress_dimension=p*p+3,
                wheel_rank=wr,cycle_rank=cr,cycle_and_wheel_rank=joined,
                missing_dimension=p*p+3-joined)


if __name__ == '__main__':
    out={'rigidity':[check_rigidity(p) for p in [3,5,7,11,13,17,19,23,29]],
         'positive_controls':[check_fixed_error(p) for p in [3,5]],
         'seams':[triangular(p) for p in [5,11,17,23,29]],
         'full_q29':triangular(29,full=True),
         'wheel_spans':[check_wheel_span(p) for p in [5,7,11,17,29]]}
    Path(__file__).with_name('independent-results.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))
    print('ALL INDEPENDENT EXACT CONTROLS PASSED')
