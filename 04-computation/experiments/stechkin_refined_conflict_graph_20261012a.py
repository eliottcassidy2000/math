"""Exact region-conflict graphs of the fixed-scale Stechkin trihexagonal tiles."""
from dataclasses import dataclass
from fractions import Fraction as F
from itertools import product, combinations
from functools import lru_cache
import finite_algebraic_witness_20261007 as proof

CHECKS=0

def need(ok, why):
    if not ok: raise ValueError(why)

def check(ok, why):
    global CHECKS
    CHECKS+=1
    if not ok: raise RuntimeError(why)

def point(p):
    need(type(p) is tuple and len(p)==2 and all(type(v) in (int,F) for v in p),
         "two exact rational coordinates")
    return tuple(F(v) for v in p)

def add(a,b): return (a[0]+b[0],a[1]+b[1])
def sub(a,b): return (a[0]-b[0],a[1]-b[1])
def scale(c,a): return (c*a[0],c*a[1])
def dot(a,b): return a[0]*b[0]+a[1]*b[1]+F(1,2)*(a[0]*b[1]+a[1]*b[0])
def norm(a): return dot(a,a)
def cross(o,a,b):
    a,b=sub(a,o),sub(b,o)
    return a[0]*b[1]-a[1]*b[0]

def hull(points):
    points=sorted(set(point(p) for p in points))
    need(len(points)>=3,"nondegenerate polygon source")
    lo=[]; hi=[]
    for p in points:
        while len(lo)>1 and cross(lo[-2],lo[-1],p)<=0: lo.pop()
        lo.append(p)
    for p in reversed(points):
        while len(hi)>1 and cross(hi[-2],hi[-1],p)<=0: hi.pop()
        hi.append(p)
    answer=tuple(lo[:-1]+hi[:-1])
    need(len(answer)>=3,"two-dimensional polygon")
    return answer

BASE=(
    tuple((F(x,2),F(y,2)) for x,y in ((1,0),(1,-1),(0,-1),(-1,0),(-1,1),(0,1))),
    tuple((F(x,2),F(y,2)) for x,y in ((1,0),(1,1),(0,1))),
    tuple((F(x,2),F(y,2)) for x,y in ((1,1),(2,1),(1,2))),
)

@dataclass(frozen=True,order=True)
class Tile:
    kind:int
    i:int
    j:int

def audit_tile(t):
    need(type(t) is Tile and all(type(v) is int for v in (t.kind,t.i,t.j))
         and 0<=t.kind<3,"exact tile kind0/1/2 and lattice address")
    return t

def vertices(t):
    audit_tile(t)
    return tuple(add(p,(t.i,t.j)) for p in BASE[t.kind])

def distance_window(left,right):
    audit_tile(left);audit_tile(right)
    return _window(left.kind,right.kind,right.i-left.i,right.j-left.j)

@lru_cache(None)
def _window(a,b,i,j):
    return polygon_window(BASE[a],tuple(add(y,(i,j)) for y in BASE[b]))


def polygon_window(left,right):
    """Distance infimum/supremum for interiors of exact full-dimensional convex hulls."""
    need(type(left) is tuple and type(right) is tuple,"two exact polygon tuples")
    left,right=hull(left),hull(right)
    polygon=hull(sub(y,x) for x in left for y in right)
    upper=max(norm(p) for p in polygon)
    lower=F(0)
    if not all(cross(x,y,(0,0))>=0 for x,y in zip(polygon,polygon[1:]+polygon[:1])):
        candidates=[]
        for x,y in zip(polygon,polygon[1:]+polygon[:1]):
            d=sub(y,x)
            t=max(F(0),min(F(1),-dot(x,d)/norm(d)))
            candidates.append(norm(add(x,scale(t,d))))
        lower=min(candidates)
    return lower,upper

def conflict(left,right):
    lo,hi=distance_window(left,right)
    return lo<1<hi

def offsets(a,b):
    need(type(a) is int and type(b) is int and 0<=a<3 and 0<=b<3,"two tile kinds")
    return _offsets(a,b)

@lru_cache(None)
def _offsets(a,b):
    return frozenset((i,j) for i,j in product(range(-2,3),repeat=2)
                     if conflict(Tile(a,0,0),Tile(b,i,j)))

def region_graph(tiles):
    need(type(tiles) is tuple,"ordered source tiles")
    for tile in tiles: audit_tile(tile)
    need(len(set(tiles))==len(tiles),"distinct source tiles")
    return frozenset((i,j) for i,j in combinations(range(len(tiles)),2)
                     if conflict(tiles[i],tiles[j]))

def period_graph(period):
    need(type(period) is int and period>=1,"positive exact period")
    tiles=tuple(Tile(s,i,j) for s in range(3)
                for i,j in product(range(period),repeat=2))
    ids={t:i for i,t in enumerate(tiles)}
    edges=set(); loops=set()
    for v,t in enumerate(tiles):
        for s in range(3):
            for i,j in offsets(t.kind,s):
                u=ids[Tile(s,(t.i+i)%period,(t.j+j)%period)]
                if u==v:loops.add(v)
                else:edges.add(tuple(sorted((u,v))))
    return tiles,frozenset(edges),frozenset(loops)


EIGHT_TILES=tuple(Tile(*p) for p in (
    (0,0,0),(0,0,1),(0,1,0),(0,1,1),
    (1,0,0),(1,0,1),(1,1,0),
    (2,-1,0),(2,-1,1),(2,0,-1),(2,0,0),
    (2,0,1),(2,1,-1),(2,1,0),
))
TRIANGLE_FOUR=tuple(Tile(2,*p) for p in
    ((-1,-1),(-1,0),(-1,1),(0,-1),(0,0),(1,-1)))+(Tile(1,0,0),)
UNIVERSAL=(1,2,3)
RIM=(7,8,11,13,12,9)
HUB=10
SMALL=(4,5,6)
SPECIAL=(0,HUB)


def tile_colour(tile):
    audit_tile(tile)
    return (0 if tile.kind==0 else 4)+2*(tile.i%2)+(tile.j%2)


def triangle_colour(tile):
    audit_tile(tile)
    need(tile.kind in (1,2),"triangle tile")
    return tile_colour(tile)-4


HEX=hull(BASE[0])


def hex_contains(p,centre):
    p=point(p)
    need(type(centre) is tuple and len(centre)==2 and
         all(type(v) is int for v in centre),"integer centre")
    p=sub(p,centre)
    return all(cross(a,b,p)>=0 for a,b in zip(HEX,HEX[1:]+HEX[:1]))


def locate(p):
    """Exact rational point owner; triangles open, shared hex vertices lex-first."""
    p=point(p)
    i,j=p[0].numerator//p[0].denominator,p[1].numerator//p[1].denominator
    a,b=p[0]-i,p[1]-j
    if 0<a<F(1,2) and 0<b<F(1,2) and a+b>F(1,2):
        return Tile(1,i,j)
    if F(1,2)<a<1 and F(1,2)<b<1 and a+b<F(3,2):
        return Tile(2,i,j)
    centres=tuple((x,y) for x in range(i-1,i+2) for y in range(j-1,j+2)
                  if hex_contains(p,(x,y)))
    need(bool(centres),"partition coverage invariant")
    x,y=min(centres)
    return Tile(0,x,y)


def point_colour(p):
    return tile_colour(locate(p))


def clip_polygon(polygon,origin,direction,sign):
    """Closed rational half-plane clipping; interiors give the split-region shapes."""
    need(type(polygon) is tuple,"polygon tuple")
    polygon=hull(polygon);origin=point(origin);direction=point(direction)
    need(direction!=(0,0) and type(sign) is int and sign in (-1,1),"oriented nonzero split")
    out=[]
    for a,b in zip(polygon,polygon[1:]+polygon[:1]):
        fa=sign*cross(origin,add(origin,direction),a)
        fb=sign*cross(origin,add(origin,direction),b)
        if fa>=0:out.append(a)
        if (fa<0<fb) or (fb<0<fa):
            out.append(add(a,scale(fa/(fa-fb),sub(b,a))))
    return hull(out)


SPLIT_COLOURS=(0,5,6,4,5,0,3,2,1,3,3,2,1,1,0)


def split_packet():
    """Same 14 source regions, with only H(0,0) split along a+b=0."""
    original=tuple(hull(vertices(t)) for t in EIGHT_TILES)
    parts=tuple(clip_polygon(original[0],(0,0),(1,-1),s) for s in (1,-1))
    polygons=parts+original[1:]
    edges=frozenset((i,j) for i,j in combinations(range(15),2)
                    if (lambda w:w[0]<1<w[1])(polygon_window(polygons[i],polygons[j])))
    return polygons,edges,SPLIT_COLOURS


def periodic_split_packet():
    """Complete period-two quotient: split exactly even/even hexes along a+b=0."""
    labels=[];polygons=[]
    for kind in range(3):
        for i,j in product(range(2),repeat=2):
            tile=Tile(kind,i,j);poly=hull(vertices(tile))
            if tile==Tile(0,0,0):
                for sign in (1,-1):
                    labels.append((tile,sign))
                    polygons.append(clip_polygon(poly,(0,0),(1,-1),sign))
            else:
                labels.append((tile,0));polygons.append(poly)
    edges=set();loops=set();lifts={}
    for a,b in combinations(range(13),2):
        for i,j in product(range(-1,2),repeat=2):
            shifted=tuple(add(p,(2*i,2*j)) for p in polygons[b])
            lo,hi=polygon_window(polygons[a],shifted)
            if lo<1<hi:
                edges.add((a,b));lifts.setdefault((a,b),(i,j,lo,hi))
    for a in range(13):
        for i,j in product(range(-1,2),repeat=2):
            shifted=tuple(add(p,(2*i,2*j)) for p in polygons[a])
            lo,hi=polygon_window(polygons[a],shifted)
            if lo<1<hi:loops.add(a)
    return tuple(labels),tuple(polygons),frozenset(edges),frozenset(loops),lifts


def verify_eight_obstruction(tiles):
    """Authenticate a translated, marked finite capacity certificate against geometry."""
    need(type(tiles) is tuple and len(tiles)==14,"fourteen marked source regions")
    edges=region_graph(tiles)
    def edge(a,b): return tuple(sorted((a,b))) in edges
    need(all(edge(a,b) for a in UNIVERSAL for b in range(14) if a!=b),"universal K3")
    need(all(edge(a,b) for a,b in combinations(SMALL,2)),"small triangle clique")
    need(all(edge(HUB,a) for a in RIM),"rim hub")
    need(all(edge(a,b)==(((RIM.index(a)-RIM.index(b))%6) in (1,5))
             for a,b in combinations(RIM,2)),"exact six-cycle rim")
    alternating=(RIM[::2],RIM[1::2])
    need(all(any(edge(a,b) for b in triple) for a in SMALL for triple in alternating),
         "each small vertex meets both maximum rim-independent sets")
    remainder=tuple(i for i in range(14) if i not in UNIVERSAL)
    need(edge(*SPECIAL),"special classes must differ")
    for a in SPECIAL:
        nonneighbors=tuple(b for b in remainder if b!=a and not edge(a,b))
        need(all(edge(b,c) for b,c in combinations(nonneighbors,2)),
             "special independent classes have size at most two")
    # Remainder: 0, hub, three-clique, six-rim. An independent set avoiding
    # specials has at most one small vertex and at most three rim vertices;
    # the only rim triples are alternating and cannot accompany a small vertex.
    need(len(remainder)==11 and set(remainder)==set((0,HUB)+SMALL+RIM),"partition")
    need(2+2+3+3<11,"four-class capacity contradiction")
    return True


def main():
    # Full finite offset universe; outside it the coordinate lower bound is 27/16.
    for a,b in product(range(3),repeat=2):
        off=offsets(a,b)
        for i,j in product(range(-4,5),repeat=2):
            window=distance_window(Tile(a,0,0),Tile(b,i,j))
            check(all(type(v) is F for v in window),"exact distance window")
            check(((i,j) in off)==(window[0]<1<window[1]),"complete offset support")
            check(window==distance_window(Tile(b,i,j),Tile(a,0,0)),"distance symmetry")
            if max(abs(i),abs(j))>=3:
                check(window[0]>=F(27,16),"exterior offset bound")
        for i,j in off:
            check(tile_colour(Tile(a,0,0))!=tile_colour(Tile(b,i,j)),"global parity upper")
    counts=tuple(tuple(len(offsets(a,b)) for b in range(3)) for a in range(3))
    check(counts==((12,12,12),(12,6,3),(12,3,6)),"all directed local counts")
    check(distance_window(Tile(0,0,0),Tile(0,0,0))==(0,1),"open hex diameter endpoint")
    check(distance_window(Tile(0,0,0),Tile(0,2,0))==(1,9),"facing boundary endpoint")
    check(norm((1,1))==3 and distance_window(Tile(0,0,0),Tile(0,1,1))==(F(3,4),7),
          "region conflict is not centre unit adjacency")

    edges=region_graph(EIGHT_TILES)
    check(len(edges)==65,"finite obstruction incidence count")
    check(verify_eight_obstruction(EIGHT_TILES),"marked capacity certificate")
    check(proof.verify_colouring(14,edges,8,tuple(map(tile_colour,EIGHT_TILES))),"eight-colour witness")
    for i,j in product(range(-2,3),repeat=2):
        translated=tuple(Tile(t.kind,t.i+i,t.j+j) for t in EIGHT_TILES)
        check(region_graph(translated)==edges,"translation preserves every conflict")
        check(verify_eight_obstruction(translated),"source-bound translated obstruction")
    remainder=(0,)+tuple(range(4,14))
    independent={k:tuple(c for c in combinations(remainder,k)
                        if not any(tuple(sorted(e)) in edges for e in combinations(c,2)))
                 for k in range(1,5)}
    check(tuple(len(independent[k]) for k in range(1,5))==(11,26,14,0),"independent-set cross-check")
    check(all(not any(v in c for c in independent[3]) for v in SPECIAL),"special capacity cross-check")
    four_edges=region_graph(TRIANGLE_FOUR)
    check(len(four_edges)==12,"seven-triangle-region witness")
    colouring,refutation=proof.colour_or_refute(7,four_edges,3)
    check(colouring is None,"triangle graph non-three-colourability")
    refutation_nodes=proof.verify_refutation(7,four_edges,3,refutation)
    check(refutation_nodes>0,"complete source-bound three-colour refutation")
    check(proof.verify_colouring(7,four_edges,4,tuple(map(triangle_colour,TRIANGLE_FOUR))),"triangle four-colour upper")
    split_polygons,split_edges,split_colours=split_packet()
    check(len(split_edges)==70,"complete split-source incidence")
    check(proof.verify_colouring(15,split_edges,7,split_colours),"finite split seven-colour witness")
    for a,b in combinations(range(15),2):
        window=polygon_window(split_polygons[a],split_polygons[b])
        check(((a,b) in split_edges)==(window[0]<1<window[1]),"every split cross-context retained")
    check(split_colours[0]!=split_colours[1],"internal seam separates distinct colours")
    split_labels,split_period_polys,split_period_edges,split_loops,split_lifts=periodic_split_packet()
    check(len(split_labels)==13 and len(split_period_edges)==63 and not split_loops,"periodic split complete quotient")
    clique=(1,2,3,4,5,6,7,8)
    check(all(e in split_period_edges for e in combinations(clique,2)),"periodic split K8")
    for a,b in combinations(clique,2):
        i,j,lo,hi=split_lifts[a,b]
        shifted=tuple(add(p,(2*i,2*j)) for p in split_period_polys[b])
        check((lo,hi)==polygon_window(split_period_polys[a],shifted) and lo<1<hi,
              "source-complete periodic clique edge lift")
    check(proof.verify_colouring(13,split_period_edges,8,
                                tuple(tile_colour(t) for t,s in split_labels)),"periodic split eight upper")
    # Complete 2-period quotient provides a second, finite positive certificate.
    periodic,period_edges,loops=period_graph(2)
    check(not loops,"two-period quotient no self-conflicts")
    check(proof.verify_colouring(len(periodic),period_edges,8,tuple(map(tile_colour,periodic))),"period upper certificate")

    # Rational coverage and exact unit-pair controls including boundary ownership.
    for a,b in product(range(-12,25),repeat=2):
        p=(F(a,12),F(b,12)); owner=locate(p)
        check(0<=point_colour(p)<8,"total rational partition colour")
        if owner.kind==0:check(hex_contains(p,(owner.i,owner.j)),"owned hex boundary")
    directions=tuple((F(a,7),F(b,7)) for a,b in product(range(-9,10),repeat=2)
                     if a*a+a*b+b*b==49)
    check(len(directions)==18,"rational unit direction family")
    for a,b in product(range(-4,5),repeat=2):
        p=(F(a,6),F(b,6))
        for d in directions:
            check(norm(d)==1 and point_colour(p)!=point_colour(add(p,d)),"actual unit-pair colour")
    for v in BASE[0]:
        check((locate(v)==Tile(0,0,0))!=(locate(scale(-1,v))==Tile(0,0,0)),"one vertex per antipodal pair")
    p,q=(F(1,2),F(0)),(F(3,2),F(0))
    check(norm(sub(p,q))==1 and point_colour(p)!=point_colour(q),"facing same-parity hex ownership")
    check(locate(p)==Tile(0,0,0) and locate(q)==Tile(0,1,0),"lexicographic facing vertex labels")

    bad=[lambda: offsets(True,0),lambda: offsets(0,1.0),
         lambda: distance_window(Tile(False,0,0),Tile(0,0,0)),
         lambda: region_graph((Tile(0,0,0),Tile(0,0,0))),
         lambda: locate((0.5,0)),lambda: locate((0,0,0)),
         lambda: period_graph(True),lambda: period_graph(0),
         lambda: verify_eight_obstruction(EIGHT_TILES[:-1]),
         lambda: verify_eight_obstruction((Tile(0,30,30),)+EIGHT_TILES[1:])]
    for action in bad:
        try:action()
        except ValueError:check(True,"malformed or wrong-source packet rejected")
        else:check(False,"bad packet accepted")
    print("fixed-scale Stechkin region-conflict graph: exact arithmetic")
    print("directed offset counts:",counts)
    print("triangle-only: 7 regions, 12 edges, chi=4; refutation nodes:",refutation_nodes)
    print("whole partition: 14 regions, 65 edges, chi=8")
    print("K3 join remainder independent-set counts:",tuple(len(independent[k]) for k in range(1,5)))
    print("capacity certificate: 11 > 2+2+3+3; arbitrary colour period excluded")
    print("two-period upper quotient:",len(periodic),"regions",len(period_edges),"edges")
    print("one hexagon split: 15 pieces, 70 edges; verified seven-colouring:",split_colours)
    print("even/even periodic split: 13 regions, 63 edges, no loops; K8 ports:",clique)
    print("global 8-colour rule: separate hex/triangle parity palettes; lexicographic hex vertices")
    print("scope: interiors monochromatic; no eight-chromatic unit-point graph")
    print("checks:",CHECKS)


if __name__=="__main__":main()
