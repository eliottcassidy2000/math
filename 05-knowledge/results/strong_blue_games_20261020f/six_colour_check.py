"""Exact controls for pinned stress and minimum palette-blocker theorems.

No repository imports. Rational star; biquadratic Moser geometry and a
good-prime rank certificate; every pinned subset and proper four-colour
boundary assignment; all 32768 pair families and nine equality banks.
"""
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations, permutations, product
from math import ceil
from pathlib import Path
import json


def need(ok, why):
    if not ok:
        raise RuntimeError(why)


def rank_mod(rows, columns, p=17):
    a = [[x % p for x in row] for row in rows]
    r = 0
    for j in range(columns):
        pivot = next((i for i in range(r,len(a)) if a[i][j]), None)
        if pivot is None:
            continue
        a[r], a[pivot] = a[pivot], a[r]
        inv = pow(a[r][j],-1,p)
        a[r] = [(x*inv)%p for x in a[r]]
        for i in range(len(a)):
            if i != r and a[i][j]:
                b = a[i][j]
                a[i] = [(x-b*y)%p for x,y in zip(a[i],a[r])]
        r += 1
    return r


def pinned_rows(points, edges, ports):
    interior = sorted(set(range(len(points)))-set(ports))
    pos = {v:i for i,v in enumerate(interior)}
    rows = []
    for u,v in edges:
        if u in ports and v in ports:
            continue
        dx,dy = points[u][0]-points[v][0], points[u][1]-points[v][1]
        row = [0]*(2*len(interior))
        for w,sign in ((u,1),(v,-1)):
            if w in pos:
                row[2*pos[w]] = sign*dx
                row[2*pos[w]+1] = -3*sign*dy
        rows.append(row)
    return rows


def peel(edges, ports, interior):
    left = set(interior)
    order = []
    while left:
        vertex = next((v for v in sorted(left)
                       if sum(v in edge and bool(set(edge)-{v} <= left|set(ports))
                              for edge in edges) <= 3), None)
        need(vertex is not None, 'relative three-degeneracy failed')
        left.remove(vertex)
        order.append(vertex)
    return order


def field(a=0,b=0,c=0,d=0):
    return tuple(map(F,(a,b,c,d)))


def add(a,b):
    return tuple(x+y for x,y in zip(a,b))


def neg(a):
    return tuple(-x for x in a)


def mul(a,b):
    out = [F(0)]*4
    for i,x in enumerate(a):
        for j,y in enumerate(b):
            out[i^j] += x*y*(-3 if i&j&1 else 1)*(-11 if i&j&2 else 1)
    return tuple(out)


def conjugate(a):
    return tuple((-x if i.bit_count()%2 else x) for i,x in enumerate(a))


def modf(x):
    return x.numerator*pow(x.denominator,-1,17)%17


def moser_control():
    one,zero = field(1),field()
    alpha,beta = field(F(1,2),F(1,2)),field(F(5,6),0,F(1,6))
    z = [zero,one,alpha,add(one,alpha),beta,mul(alpha,beta),
         mul(beta,add(one,alpha))]
    edges = []
    for u,v in combinations(range(7),2):
        delta = add(z[u],neg(z[v]))
        if mul(delta,conjugate(delta)) == one:
            edges.append((u,v))
    need(len(edges)==11,'exact Euclidean spindle edge count')
    # At 17 take s=sqrt(-3)->4 theta, t=sqrt(-11)->6 theta,
    # theta^2=3. The norm is x^2-3y^2; all denominators are units.
    points = [(modf(a+4*d),modf(4*b+6*c)) for a,b,c,d in z]
    residue_edges = [(u,v) for u,v in combinations(range(7),2)
                     if ((points[u][0]-points[v][0])**2
                         -3*(points[u][1]-points[v][1])**2)%17==1]
    need(edges==residue_edges,'accidental residue edge')
    full = {q:[word for word in product(range(q),repeat=7)
               if all(word[u]!=word[v] for u,v in edges)] for q in (3,4)}
    need(not full[3] and bool(full[4]),'spindle chromatic controls')
    checked = Counter()
    for mask in range(128):
        ports = {u for u in range(7) if mask>>u&1}
        interior = set(range(7))-ports
        rows = pinned_rows(points,edges,ports)
        rank = rank_mod(rows,2*len(interior))
        if rank != len(rows):
            continue
        checked['independent_pin_sets'] += 1
        peel(edges,ports,interior)
        port_order = sorted(ports)
        projections = {tuple(word[u] for u in port_order) for word in full[4]}
        for word in product(range(4),repeat=len(ports)):
            colour = dict(zip(port_order,word))
            if all(colour[u]!=colour[v] for u,v in edges if u in ports and v in ports):
                need(word in projections,'four-colour boundary extension failed')
                checked['proper_boundary_assignments'] += 1
    one_pin_rows = pinned_rows(points,edges,{0})
    need(len(one_pin_rows)==rank_mod(one_pin_rows,12)==11,
         'one-pin sharpness rank')
    diamond_edges = [(u,v) for u,v in edges if u<4 and v<4]
    diamond_rows = pinned_rows(points[:4],diamond_edges,{0,3})
    diamond_rank = rank_mod(diamond_rows,4)
    need((len(diamond_rows),diamond_rank)==(5,4),'diamond pinned stress')
    need(not any(c[0]==0 and c[3]==1 and all(c[u]!=c[v] for u,v in diamond_edges)
                 for c in product(range(3),repeat=4)), 'diamond unequal tips extend')
    return dict(exact_edges=edges, residue_points=points,
                proper_four_colourings=len(full[4]), four_colour_checks=dict(checked),
                one_pin=dict(rows=11,rank=11,three_colourings=0),
                diamond=dict(rows=5,rank=4,stress=1,critical_lower_bound=1))


def rational_star():
    ports = [(F(1-t*t,1+t*t),F(2*t,1+t*t)) for t in range(1,7)]
    need(len(set(ports))==6,'port collision')
    need(all(x*x+y*y==1 for x,y in ports),'nonunit radial edge')
    need(all((x-u)**2+(y-v)**2!=1 for (x,y),(u,v) in combinations(ports,2)),
         'unit edge between ports')
    need(any(x*v-y*u for (x,y),(u,v) in combinations(ports,2)), 'rank below two')
    return dict(ports=[[str(x),str(y)] for x,y in ports],
                interior_vertices=1,port_edges=6,distinct_forbidden_colours=6,
                pinned_rank=2,pinned_stress=4,critical_lower_bound=4,
                whole_graph_chromatic_number=2)


def palette_control():
    pairs = list(combinations(range(6),2))
    options = sorted({tuple(tuple(sorted(p[2*i:2*i+2])) for i in range(3))
                      for p in permutations(range(6))})
    need(len(options)==90,'option count')
    blocks = {(i,p):{j for j,o in enumerate(options) if o[i]==p}
              for i in range(3) for p in pairs}
    need({len(v) for v in blocks.values()}=={6},'block sizes')
    full = (1<<15)-1
    neighbours = [sum(1<<j for j,q in enumerate(pairs) if set(p)&set(q))
                  for p in pairs]
    maxima = []
    large = []
    for a in range(1,1<<15):
        common = full
        for i in range(15):
            if a>>i&1:
                common &= neighbours[i]
        if common:
            total = a.bit_count()+common.bit_count()
            need(total<=10,'cross-intersection bound')
            if total==10:
                maxima.append(a)
        if a.bit_count()>=5 and common.bit_count()>=5:
            large.append(a)
    stars = [sum(1<<j for j,p in enumerate(pairs) if c in p) for c in range(6)]
    need(sorted(large)==sorted(stars),'five-by-five rigidity')
    banks = [[(i,p) for p in pairs] for i in range(3)]
    banks += [[(i,p) for i in range(3) for p in pairs if c in p] for c in range(6)]
    for bank in banks:
        covered = Counter(j for b in bank for j in blocks[b])
        need(len(bank)==15 and covered==Counter({j:1 for j in range(90)}),
             'minimum bank does not partition options')
    # A concentrated minimum bank does not force failure of actual colouring:
    # fifteen independent exterior vertices can each choose outside its pair.
    lift = [next(c for c in range(6) if c not in p) for p in pairs]
    need(all(c not in p for c,p in zip(lift,pairs)),'palette-failure hostile')
    # Common-star bank: complete tripartite exterior with five vertices per
    # part. Colour0 is forbidden everywhere and each part sees all exclusions1..5.
    # Each part needs at least two colours, while distinct parts share none.
    allowed = [set(c for c in range(1,6) if mask>>(c-1)&1)
               for mask in range(32) if mask.bit_count()>=2]
    need(not any(not(a&b or a&c or b&c) for a,b,c in product(allowed,repeat=3)),
         'tripartite failed-precolour control unexpectedly extends')
    exterior = [(i,j,6+5*i+j-1) for i in range(3) for j in range(1,6)]
    graph = {(min(u,v),max(u,v)) for i,j,u in exterior for k,l,v in exterior if i!=k}
    graph |= {(0,v) for i,j,v in exterior}
    graph |= {(j,v) for i,j,v in exterior}
    four = [3]*6+[i for i,j,v in exterior]
    need(all(four[u]!=four[v] for u,v in graph),'tripartite four-colouring')
    clique = (0,6,11,16)
    need(all((u,v) in graph for u,v in combinations(clique,2)), 'tripartite K4')
    return dict(options=90,blocker_types=45,options_per_blocker=6,
                pair_families=32768,cross_extremal_first_families=len(maxima),
                large_cross_families=len(large),minimum_cover_banks=9,
                concentrated_actual_lift=lift,
                star_graph=dict(exterior_vertices=15,classes=[5,5,5],
                                common_neighbours_per_port_pair=3,
                                six_precolour_extensions=0,ordinary_chromatic=4))


if __name__=='__main__':
    source=Path(__file__).resolve()
    result=dict(status='PASS',source_sha256=sha256(source.read_bytes().replace(b'\r\n',b'\n')).hexdigest(),
                rational_star=rational_star(),moser=moser_control(),palette=palette_control())
    payload=json.dumps(result,indent=2,sort_keys=True)+'\n'
    source.with_suffix('.json').write_text(payload,encoding='utf-8',newline='\n')
    print(payload,end='')
