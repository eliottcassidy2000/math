"""Source-bound colour-difference tables, palette gauges and shade repairs.

Only the rational translate constructor claims geometric completeness itself.
An abstract layered edge bank must be authenticated by its geometric caller.
"""
from collections import Counter
from dataclasses import dataclass, replace
from fractions import Fraction
from itertools import combinations, permutations, product
import json

import snark_boundary_transfer_20261018a as native
import plane_six_reframes_parity_20261011b as parity
import exceptional_factor_support_20261018a as source8

Graph = native.Graph
require, nat = native.require, native.nat


@dataclass(frozen=True)
class Layers:
    base: Graph
    count: int
    whole: Graph


def audit_layers(p):
    require(type(p) is Layers, "Layers required")
    native.audit_graph(p.base); nat(p.base.n,1); nat(p.count,1)
    native.audit_graph(p.whole)
    n = p.base.n
    require(p.whole.n == n*p.count, "complete layer vertex bank")
    for i in range(p.count):
        actual = tuple((a%n,b%n) for a,b in p.whole.edges if a//n == b//n == i)
        require(actual == p.base.edges, "each source copy must retain its exact induced graph")
    return p


def layered(base, count, cross):
    native.audit_graph(base); nat(base.n,1); nat(count,1)
    require(type(cross) is tuple, "cross-edge tuple")
    for e in cross:
        require(type(e) is tuple and len(e)==2, "cross edge")
        nat(e[0]); nat(e[1])
        require(e[0]//base.n != e[1]//base.n, "cross edge is internal")
    edges = tuple((i*base.n+a,i*base.n+b) for i in range(count) for a,b in base.edges)+cross
    return audit_layers(Layers(base,count,native.graph(base.n*count,edges)))


@dataclass(frozen=True)
class Pattern:
    layers: Layers
    colours: tuple
    events: tuple
    matrices: tuple


def compile_pattern(p, colours):
    audit_layers(p); native.audit_colouring(p.base,6,colours)
    n = p.base.n
    events = tuple((a//n,b//n,colours[a%n],colours[b%n],a,b)
                   for a,b in p.whole.edges if a//n != b//n)
    rows = {}
    for i,j,a,b,_,_ in events:
        rows.setdefault((i,j),[0]*6)[a] |= 1<<b
    matrices = tuple((i,j,tuple(row)) for (i,j),row in sorted(rows.items()))
    return Pattern(p,colours,events,matrices)


def audit_pattern(p):
    require(type(p) is Pattern, "Pattern required")
    audit_layers(p.layers); native.audit_colouring(p.layers.base,6,p.colours)
    for name in ("events","matrices"):
        require(type(getattr(p,name)) is tuple, "tuple event bank")
    for row in p.events:
        require(type(row) is tuple and len(row)==6, "event row")
        for x in row: nat(x)
    for row in p.matrices:
        require(type(row) is tuple and len(row)==3, "matrix row")
        nat(row[0]); nat(row[1])
        require(type(row[2]) is tuple and len(row[2])==6, "six matrix masks")
        for x in row[2]: nat(x)
    require(p == compile_pattern(p.layers,p.colours), "pattern/source reconstruction mismatch")
    return p


def qpoint(x,y):
    require(type(x) in (int,Fraction) and type(y) in (int,Fraction), "exact rational point")
    return Fraction(x),Fraction(y)


def audit_points(points):
    require(type(points) is tuple and len(points)>0, "nonempty point tuple")
    for p in points:
        require(type(p) is tuple and len(p)==2 and all(type(c) is Fraction for c in p),
                "canonical Fraction point")
    require(len(set(points))==len(points), "coincident source points")


def sub(a,b): return a[0]-b[0],a[1]-b[1]
def add(a,b): return a[0]+b[0],a[1]+b[1]
def norm(a): return a[0]*a[0]+a[1]*a[1]


def translated(points, offsets):
    audit_points(points); audit_points(offsets)
    require(offsets[0] == qpoint(0,0), "source layer must remain literally fixed")
    allpoints = tuple(add(p,t) for t in offsets for p in points)
    require(len(set(allpoints)) == len(allpoints), "coincident translated vertices")
    base = native.graph(len(points),((a,b) for a,b in combinations(range(len(points)),2)
                                    if norm(sub(points[b],points[a]))==1))
    whole = native.graph(len(allpoints),((a,b) for a,b in combinations(range(len(allpoints)),2)
                                        if norm(sub(allpoints[b],allpoints[a]))==1))
    return audit_layers(Layers(base,len(offsets),whole))


def difference_bank(points, colours):
    """Each vector retains all ordered anchor pairs; not merely a scalar norm."""
    audit_points(points)
    require(type(colours) is tuple and len(colours)==len(points), "source colour tuple")
    for c in colours: nat(c); require(c<6,"six source colours")
    rows = {}
    for u,v in product(range(len(points)),repeat=2):
        key = colours[u],colours[v],sub(points[v],points[u])
        rows.setdefault(key,[]).append((u,v))
    return tuple((a,b,d,tuple(pairs)) for (a,b,d),pairs in sorted(rows.items()))


def geometric_pattern(points, offsets, colours):
    layers = translated(points,offsets)
    p = compile_pattern(layers,colours)
    bank = difference_bank(points,colours)
    events = []
    n = len(points)
    for i,j in combinations(range(len(offsets)),2):
        delta = sub(offsets[j],offsets[i])
        for a,b,d,pairs in bank:
            if norm(add(d,delta))==1:
                events.extend((i,j,a,b,i*n+u,j*n+v) for u,v in pairs)
    require(tuple(sorted(events,key=lambda e:(e[4],e[5]))) == p.events,
            "difference bank differs from complete geometric pair census")
    return p, bank


IDENTITY = tuple(range(6))
FULL = tuple(permutations(range(6)))
PAIRED = tuple(sorted(tuple(2*order[c//2]+((c%2)^bits[c//2]) for c in range(6))
                       for order in permutations(range(3)) for bits in product((0,1),repeat=3)))


def audit_gauges(p, gauges):
    audit_pattern(p)
    require(type(gauges) is tuple and len(gauges)==p.layers.count,"one gauge per layer")
    for row in gauges:
        require(type(row) is tuple and len(row)==6,"six permutation entries")
        for x in row: nat(x)
        require(tuple(sorted(row))==IDENTITY,"actual palette permutation")
    require(gauges[0]==IDENTITY,"source palette must remain fixed")
    require(all(gauges[i][a]!=gauges[j][b] for i,j,a,b,_,_ in p.events),"gauge edge conflict")
    result = tuple(gauge[c] for gauge in gauges for c in p.colours)
    native.audit_colouring(p.layers.whole,6,result)
    require(result[:len(p.colours)]==p.colours,"lost source colouring")
    return result


def candidates(p, mode):
    require(type(mode) is str and mode in ("full","paired"),"gauge mode")
    used = tuple(sorted(set(p.colours)))
    seen = {}
    for g in FULL if mode=="full" else PAIRED:
        seen.setdefault(tuple(g[c] for c in used),g)
    return tuple(seen.values())


def find_gauges(p, mode="full"):
    """Complete for fixed-layer-pattern gauges only, not all six-colourings."""
    audit_pattern(p)
    options = candidates(p,mode)
    chosen = [IDENTITY]
    def visit():
        if len(chosen)==p.layers.count:
            return tuple(chosen)
        j = len(chosen)
        relevant = tuple((i,a,b) for i,k,a,b,_,_ in p.events if k==j)
        for g in options:
            if all(chosen[i][a]!=g[b] for i,a,b in relevant):
                chosen.append(g)
                answer = visit()
                if answer is not None: return answer
                chosen.pop()
        return None
    answer = visit()
    if answer is not None: audit_gauges(p,answer)
    return answer


def count_two_gauges(p, mode="full"):
    audit_pattern(p); require(p.layers.count==2,"two layers required")
    return sum(all(p.colours[u%p.layers.base.n]!=g[b] for _,_,_,b,u,_ in p.events)
               for g in candidates(p,mode))


@dataclass(frozen=True)
class Repair:
    pattern: Pattern
    status: str
    colours: tuple
    count: int
    obstruction: tuple


def compile_repair(p):
    """Exact within the retained vertex-pair labels, allowing class splitting."""
    audit_pattern(p)
    n=p.layers.base.n; g=p.layers.whole
    pg=parity.Graph(tuple(range(g.n)),g.edges)
    labels=tuple(p.colours[v%n]//2 for v in range(g.n))
    boundary=tuple(range(n)); bits=tuple(c%2 for c in p.colours)
    packet=parity.compile_boundary(pg,labels,boundary)
    if packet.odd_cycle:
        return Repair(p,"odd_cycle",(),0,packet.odd_cycle)
    for equation in parity.constraints(packet):
        if bits[equation.start]^bits[equation.end] != equation.parity:
            return Repair(p,"source_path",(),0,equation.path)
    colours=parity.lift(packet,bits)
    free=set(packet.roots)-{packet.roots[v] for v in boundary}
    native.audit_colouring(g,6,colours)
    require(colours[:n]==p.colours,"source shade lost")
    return Repair(p,"lifted",colours,2**len(free),())


def audit_repair(r):
    require(type(r) is Repair,"Repair required")
    audit_pattern(r.pattern)
    require(type(r.status) is str and r.status in ("lifted","odd_cycle","source_path"),"repair status")
    nat(r.count)
    for row in (r.colours,r.obstruction):
        require(type(row) is tuple,"typed repair tuple")
        for x in row: nat(x)
    require(r==compile_repair(r.pattern),"repair differs from source reconstruction")
    return r


PAIRINGS=tuple(native.perfect_pairings(tuple(range(6))))


def audit_pairing(matching):
    require(type(matching) is tuple and len(matching)==3,"three marked colour pairs")
    for pair in matching:
        require(type(pair) is tuple and len(pair)==2,"colour pair")
        for c in pair: nat(c)
        require(pair[0]<pair[1],"ordered shade orientation")
    require(tuple(sorted(matching))==matching and tuple(sorted(c for pair in matching for c in pair))==IDENTITY,
            "canonical perfect matching of six colours")
    return matching


@dataclass(frozen=True)
class PairRepair:
    pattern: Pattern
    matching: tuple
    status: str
    colours: tuple
    count: int
    obstruction: tuple


def pairing_repair(p,matching):
    audit_pattern(p); audit_pairing(matching)
    decode=tuple(c for pair in matching for c in pair)
    encode={c:i for i,c in enumerate(decode)}
    encoded=compile_pattern(p.layers,tuple(encode[c] for c in p.colours))
    result=compile_repair(encoded)
    colours=tuple(decode[c] for c in result.colours)
    if colours:
        native.audit_colouring(p.layers.whole,6,colours)
        require(colours[:len(p.colours)]==p.colours,"pairing transport changed source")
    return PairRepair(p,matching,result.status,colours,result.count,result.obstruction)


def audit_pair_repair(r):
    require(type(r) is PairRepair,"PairRepair required")
    audit_pattern(r.pattern); audit_pairing(r.matching)
    require(type(r.status) is str,"status string")
    nat(r.count)
    for row in (r.colours,r.obstruction):
        require(type(row) is tuple,"receipt tuple")
        for x in row: nat(x)
    require(r==pairing_repair(r.pattern,r.matching),"pairing receipt source mismatch")
    return r


def pairing_sweep(p):
    """Complete fifteen-chart cover only; overlapping fibre counts are not summed."""
    audit_pattern(p)
    return tuple(pairing_repair(p,m) for m in PAIRINGS)


def copy_flip(p):
    """Signed XOR gauge: one global shade flip per copy, source flip zero."""
    audit_pattern(p)
    equations=tuple((i,j,1^(a%2)^(b%2),u,v) for i,j,a,b,u,v in p.events if a//2==b//2)
    adj=[[] for _ in range(p.layers.count)]
    for index,(i,j,s,_,_) in enumerate(equations):
        adj[i].append((j,s,index)); adj[j].append((i,s,index))
    bits=[-1]*len(adj); paths={}; components=0
    for seed in range(len(adj)):
        if bits[seed]>=0: continue
        components+=1; bits[seed]=0; paths[seed]=()
        queue=[seed]
        for i in queue:
            for j,s,index in adj[i]:
                if bits[j]<0:
                    bits[j]=bits[i]^s; paths[j]=paths[i]+(index,); queue.append(j)
                elif bits[i]^bits[j]!=s:
                    # A symmetric difference of rooted paths and this edge is
                    # an inconsistent equation set. It is not an actual unit cycle.
                    bad=tuple(sorted(set(paths[i])^set(paths[j])^{index}))
                    degrees=Counter(); charge=0
                    for k in bad:
                        a,b,t,_,_=equations[k]; degrees[a]+=1; degrees[b]+=1; charge^=t
                    require(charge==1 and all(v%2==0 for v in degrees.values()),"signed holonomy receipt")
                    return {"status":"gauge_blocked","equations":equations,"inconsistent_rows":bad}
    colours=tuple(c^bits[i] for i in range(p.layers.count) for c in p.colours)
    native.audit_colouring(p.layers.whole,6,colours)
    require(colours[:len(p.colours)]==p.colours,"source flip")
    return {"status":"lifted","bits":tuple(bits),"colours":colours,
            "count":2**(components-1),"equations":equations}


def canonical_patterns(g):
    native.audit_graph(g)
    edges=set(g.edges); out=[]
    def visit(row):
        if len(row)==g.n: out.append(tuple(row)); return
        v=len(row)
        for c in range(min(6,max(row,default=-1)+2)):
            if all(row[u]!=c for u in range(v) if (u,v) in edges): visit(row+[c])
    visit([])
    return tuple(out)


def complete_two_layer_count(p):
    """Full source-fixed six-colour count via every child partition + injection."""
    audit_pattern(p); require(p.layers.count==2,"two copies")
    n=p.layers.base.n; total=0; first=None
    cross=tuple((u%n,v%n) for u,v in p.layers.whole.edges if u//n != v//n)
    for pattern in canonical_patterns(p.layers.base):
        rank=1+max(pattern)
        for colours in permutations(range(6),rank):
            child=tuple(colours[c] for c in pattern)
            if all(p.colours[a]!=child[b] for a,b in cross):
                total+=1
                if first is None: first=p.colours+child
    if first is not None: native.audit_colouring(p.layers.whole,6,first)
    return total,first


CHECKS=0
def check(ok,message):
    global CHECKS
    CHECKS+=1
    if not ok: raise RuntimeError(message)


def rejects(fn):
    try: fn()
    except ValueError: check(True,"reject")
    else: check(False,"malformed input accepted")


def main():
    check(len(FULL)==720 and len(set(PAIRED))==48,"exact palette groups")
    # Six-cell sharp distinction between full permutations and paired gauges.
    base=native.graph(6,())
    cross=tuple((a,6) for a in range(1,6))+((1,7),)
    sharp=compile_pattern(layered(base,2,cross),IDENTITY)
    check(find_gauges(sharp,"paired") is None,"six constraints block paired subgroup")
    full=find_gauges(sharp)
    check(full is not None and count_two_gauges(sharp)==96,"full S6 survives")
    for i,a,b in ((1,a,b) for a in range(6) for b in range(6)):
        check(sum(g[a]==b for g in PAIRED)==8,"one constraint covers exactly eight gauges")
    # A one-colour child has six distinct images, not 720 or48 colourings.
    unused=compile_pattern(layered(native.graph(1,()),2,()),(0,))
    check(count_two_gauges(unused)==6 and count_two_gauges(unused,"paired")==6,"unused palette quotient")
    # Real rational translate: shade-fixed failure is not full-gauge failure.
    rational,bank=geometric_pattern((qpoint(0,0),qpoint(2,0)),(qpoint(0,0),qpoint(1,0)),(0,1))
    blocked=compile_repair(rational)
    check(blocked.status=="source_path" and len(blocked.obstruction)==3,"fixed source parity obstruction")
    check(find_gauges(rational,"paired") is not None,"changing pair label repairs actual geometry")
    rational_sweep=pairing_sweep(rational)
    check(Counter(r.status for r in rational_sweep)==Counter(lifted=12,source_path=3),
          "all fifteen matching charts on the actual rational path")
    for phi in canonical_patterns(rational.layers.base):
        p,b=geometric_pattern((qpoint(0,0),qpoint(2,0)),(qpoint(0,0),qpoint(1,0)),phi)
        check(len(b)==sum(1 for _ in difference_bank((qpoint(0,0),qpoint(2,0)),phi)),"difference provenance")
    square=(qpoint(0,0),qpoint(1,0),qpoint(1,1),qpoint(0,1))
    geometric_rows=[]
    for phi in canonical_patterns(native.graph(4,((0,1),(1,2),(2,3),(0,3)))):
        p,bank=geometric_pattern(square,(qpoint(0,0),qpoint(2,0)),phi)
        check(find_gauges(p) is not None,"rational square gauge witness")
        geometric_rows.append({"source":phi,"events":len(p.events),"difference_vectors":len(bank)})
    # Uniform copy flips can fail while independent vertex shades succeed.
    tiny=compile_pattern(layered(native.graph(2,()),3,((0,2),(3,4),(1,5))),(0,0))
    check(copy_flip(tiny)["status"]=="gauge_blocked","odd copy-quotient holonomy")
    check(compile_repair(tiny).status=="lifted","three disjoint edges still repair")
    # Even full S6 fixed-pattern gauges need not describe all colourings.
    n=6; cross=[]
    neighbours={i:tuple(j for j in range(7) if j!=i) for i in range(7)}
    for i,j in combinations(range(7),2):
        cross.append((6*i+neighbours[i].index(j),6*j+neighbours[j].index(i)))
    matching=compile_pattern(layered(native.graph(6,()),7,tuple(cross)),(0,)*6)
    check(find_gauges(matching) is None,"K7 obstruction applies only to monochromatic-copy gauges")
    repaired=compile_repair(matching)
    check(repaired.status=="lifted" and repaired.count==2**15,"42-vertex matching has full source-fixed shade lift")
    # Conversely all fifteen charts can fail while independent full gauges succeed.
    forced=compile_pattern(layered(native.graph(6,()),3,
                 tuple((c,6) for c in range(6) if c!=1)+
                 tuple((c,12) for c in range(6) if c!=2)),IDENTITY)
    failed_sweep=pairing_sweep(forced)
    check(all(r.status=="source_path" for r in failed_sweep),"all pair charts fail to send0 to both1and2")
    check(find_gauges(forced) is not None,"full S6 independently repairs the two copied layers")
    # An overlap control: fifteen size-two fibres cover only six child colours.
    chart_union=set()
    for m in PAIRINGS:
        chart_union.update(next(pair for pair in m if 0 in pair))
    check(sum(r.count for r in pairing_sweep(unused))==30 and len(chart_union)==6,
          "pairing fibre counts cannot be added as union mass")
    # Complete small universe: every source graph n<=3, every cross-edge subset.
    universes=0
    for n in range(1,4):
        possible=tuple(combinations(range(n),2))
        for mask in range(1<<len(possible)):
            base=native.graph(n,(e for i,e in enumerate(possible) if mask>>i&1))
            patterns=canonical_patterns(base)
            possible_cross=tuple((a,n+b) for a,b in product(range(n),repeat=2))
            for crossmask in range(1<<len(possible_cross)):
                layers=layered(base,2,tuple(e for i,e in enumerate(possible_cross) if crossmask>>i&1))
                for phi in patterns:
                    p=compile_pattern(layers,phi); universes+=1
                    count,witness=complete_two_layer_count(p)
                    direct=0; fixed_pairs=0
                    labels=tuple(c//2 for c in phi)
                    for child in product(range(6),repeat=n):
                        colours=phi+child
                        if all(colours[a]!=colours[b] for a,b in layers.whole.edges):
                            direct+=1
                            fixed_pairs+=tuple(c//2 for c in child)==labels
                    check(count==direct and count>0,"complete partition/injection union versus actual colourings")
                    repair=compile_repair(p)
                    check(repair.count==fixed_pairs,"exact source-fixed shade fibre")
                    if repair.status!="lifted":
                        audit_repair(repair)
    # Sharp projected13-pair control, abstract structure matching parent's geometric packet.
    base=native.graph(8,source8.BASE)
    clique=set(range(8))-{1}
    missing=(2,4)
    extra=tuple(e for e in combinations(sorted(clique),2) if e not in base.edges and e!=missing)
    check(len(extra)==13 and (2,6) in extra,"sharp projected bank")
    oriented=tuple((6,2) if e==(2,6) else e for e in extra)
    cross=tuple((a,8*(i+1)+b) for i,(a,b) in enumerate(oriented))
    star=layered(base,14,cross)
    rows=[]; sweep_status=Counter(); sweep_separating=Counter()
    for phi in source8.partitions():
        p=compile_pattern(star,phi)
        flip=copy_flip(p); fine=compile_repair(p)
        check(flip["status"]=="lifted" and fine.status=="lifted","all570 source patterns survive without retraction")
        same=sum(phi[a]//2==phi[b]//2 for a,b in oriented)
        check(flip["count"]==2**(13-same),"exact forest flip count")
        check(fine.count>=flip["count"],"component freedom contains copy flips")
        rows.append({"source":phi,"copy_flip_count":flip["count"],"shade_count":fine.count,
                     "copy_bits":flip["bits"]})
        for repair in pairing_sweep(p):
            check(repair.status=="lifted","all fifteen chart witnesses on the bridge forest")
            sweep_status[repair.status]+=1
            if phi[2]!=phi[4]: sweep_separating[repair.status]+=1
    check(len(rows)==570,"complete eight-port source support")
    # Hom(C5,B) has its vertex2 image in the old C5; none is adjacent to6.
    images=[]
    for word in product(range(8),repeat=5):
        if all(tuple(sorted((word[i],word[(i+1)%5]))) in base.edges for i in range(5)):
            images.append(word)
    check(len(images)==10 and all(set(w)==set(range(5)) for w in images),"complete C5 homomorphism barrier")
    check(not any(tuple(sorted((6,w[2]))) in base.edges for w in images),"no source-fixing retraction")
    # Empty source support / odd-cycle controls are kept separate.
    odd=compile_pattern(layered(native.graph(1,()),3,((0,1),(1,2),(0,2))),(0,))
    check(compile_repair(odd).status=="odd_cycle","true retained-pair odd-cycle rejection")
    check(all(r.status=="odd_cycle" for r in pairing_sweep(odd)),
          "all matching charts reject an actual equilateral triangle")
    check(find_gauges(odd) is not None,"three-colour triangle survives all-chart failure")
    rejects(lambda:audit_pattern(replace(sharp,colours=(False,1,2,3,4,5))))
    rejects(lambda:audit_pattern(replace(sharp,events=((False,1,0,1,0,7),)+sharp.events[1:])))
    rejects(lambda:audit_pattern(replace(sharp,matrices=())))
    rejects(lambda:audit_gauges(sharp,(IDENTITY,(False,1,2,3,4,5))))
    rejects(lambda:audit_gauges(sharp,((1,0,2,3,4,5),full[1])))
    rejects(lambda:audit_repair(replace(blocked,count=False)))
    rejects(lambda:audit_repair(replace(blocked,obstruction=(0,2,0))))
    rejects(lambda:qpoint(0.0,0))
    rejects(lambda:translated((qpoint(0,0),qpoint(1,0)),(qpoint(0,0),qpoint(1,0))))
    rejects(lambda:find_gauges(sharp,False))
    rejects(lambda:pairing_repair(sharp,((False,1),(2,3),(4,5))))
    rejects(lambda:audit_pair_repair(replace(rational_sweep[0],count=False)))
    rejects(lambda:audit_pair_repair(replace(rational_sweep[0],matching=PAIRINGS[-1])))
    return {"checks":CHECKS,"status":"PROVED typed finite interfaces; FINITE-EXACT controls; plane target OPEN",
            "small_complete_source_instances":universes,
            "paired_gauge_hostile":{"events":sharp.events,"full_fixed_pattern_count":96,
                                    "full_gauge":full,"minimum_forbidden_cells":6},
            "rational_path":{"edges":rational.layers.whole.edges,"source":rational.colours,
                              "shade_obstruction":blocked.obstruction,"paired_gauge":find_gauges(rational,"paired"),
                              "pairing_chart_status":dict(Counter(r.status for r in rational_sweep))},
            "rational_squares":geometric_rows,
            "full_gauge_false_negative":{"vertices":42,"edges":matching.layers.whole.edges,
                    "actual_shade_witness":repaired.colours,"source_fixed_shade_count":repaired.count},
            "thirteen_pair_star":{"vertices":star.whole.n,"edge_count":len(star.whole.edges),
                   "extra_projected_pairs":extra,"oriented_spokes":oriented,"rows":rows,
                   "all_pairing_status":dict(sweep_status),"separating_pairing_status":dict(sweep_separating)},
            "all_pairings_false_negative":{"edges":forced.layers.whole.edges,"full_gauge":find_gauges(forced),
                   "failed_gauges":[{"matching":r.matching,"status":r.status,"path":r.obstruction}
                                     for r in failed_sweep]},
            "scope":"Gauge failure is not full-colouring failure; difference existence matrices retain their actual edge provenance."}


if __name__=="__main__":
    print(json.dumps(main(),sort_keys=True,indent=2))
