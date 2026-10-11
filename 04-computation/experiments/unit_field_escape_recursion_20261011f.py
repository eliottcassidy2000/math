"""Exact geometric spindle closure, field escape, and retained cross edges.

No six-colour equality fragment or seven-chromatic plane graph is claimed.
"""
from dataclasses import dataclass, replace
from fractions import Fraction as F
from functools import lru_cache
from itertools import combinations, product
from math import isqrt, prod
import json
import finite_algebraic_witness_20261007 as proof

CHECKS = 0
ZERO = (F(0),)*4
ONE = (F(1),F(0),F(0),F(0))
A = (F(1,2),F(1,2),F(0),F(0))


def need(x, why):
    if not x:
        raise ValueError(why)


def check(x, why):
    global CHECKS
    CHECKS += 1
    if not x:
        raise RuntimeError(why)


def natural(x, minimum=0):
    need(type(x) is int and x >= minimum, "exact integer in stated domain")
    return x


def point(z):
    need(type(z) is tuple and len(z) == 4
         and all(type(c) in (int,F) for c in z), "four exact coefficients")
    return tuple(F(c) for c in z)


def add(x,y):
    return tuple(a+b for a,b in zip(x,y))


def scale(c,x):
    return tuple(c*a for a in x)


def sub(x,y):
    return add(x,scale(-1,y))


@dataclass(frozen=True)
class Field:
    D: int

    def audit(self):
        natural(self.D,1)
        need(isqrt(self.D)**2 != self.D and isqrt(3*self.D)**2 != 3*self.D,
             "independent quadratic generators")

    def mul(self,x,y):
        self.audit()
        x,y=point(x),point(y)
        z=[F(0)]*4
        for i,a in enumerate(x):
            for j,b in enumerate(y):
                z[i^j] += a*b*(-3 if i&j&1 else 1)*(-self.D if i&j&2 else 1)
        return tuple(z)

    def norm(self,x):
        x=point(x)
        return self.mul(x,(x[0],-x[1],-x[2],x[3]))


def coordinates(points):
    need(type(points) is tuple, "ordered exact point tuple")
    points=tuple(point(z) for z in points)
    need(len(set(points))==len(points), "distinct geometric vertices")
    return points


def unit_edges(field,points):
    need(type(field) is Field, "typed quadratic field")
    field.audit()
    points=coordinates(points)
    return _unit_edges(field,points)


@lru_cache(None)
def _unit_edges(field,points):
    return frozenset((i,j) for i,j in combinations(range(len(points)),2)
                     if field.norm(sub(points[i],points[j]))==ONE)


def edge_set(n,edges):
    need(type(edges) is frozenset, "canonical edge set")
    for e in edges:
        need(type(e) is tuple and len(e)==2 and all(type(i) is int for i in e)
             and 0<=e[0]<e[1]<n, "exact ordered edge")
    return edges


@dataclass(frozen=True)
class Seed:
    field: Field
    points: tuple
    tip: int


def audit_seed(seed):
    need(type(seed) is Seed, "typed geometric seed")
    pts=coordinates(seed.points)
    need(type(seed.field) is Field, "typed field")
    seed.field.audit()
    natural(seed.tip,1)
    need(seed.tip<len(pts) and pts[0]==ZERO, "retained root zero and distinct tip")
    norm=seed.field.norm(pts[seed.tip])
    need(norm[1:]==(0,0,0) and norm[0]>F(1,4), "rational squared tip distance >1/4")
    need(4*norm[0]-1==seed.field.D, "field matches supplied rotation radicand")
    return unit_edges(seed.field,pts),norm[0]


@dataclass(frozen=True)
class Spindle:
    seed: Seed
    points: tuple
    maps: tuple
    closing: tuple
    skeleton: frozenset
    full_edges: frozenset


def build_spindle(seed):
    edges,N=audit_seed(seed)
    rho=(F(2*N-1,2*N),F(0),F(1,2*N),F(0))
    rotated=tuple(seed.field.mul(rho,z) for z in seed.points)
    pts=tuple(dict.fromkeys(seed.points+rotated))
    index={z:i for i,z in enumerate(pts)}
    maps=(tuple(index[z] for z in seed.points),tuple(index[z] for z in rotated))
    closing=tuple(sorted((maps[0][seed.tip],maps[1][seed.tip])))
    need(closing[0]!=closing[1], "distinct closing tips")
    skeleton=frozenset(tuple(sorted((mp[u],mp[v]))) for mp in maps for u,v in edges)
    skeleton|={closing}
    result=Spindle(seed,pts,maps,closing,skeleton,unit_edges(seed.field,pts))
    audit_spindle(result)
    return result


def audit_spindle(packet):
    need(type(packet) is Spindle, "typed spindle packet")
    edges,N=audit_seed(packet.seed)
    pts=coordinates(packet.points)
    need(type(packet.maps) is tuple and len(packet.maps)==2, "two retained vertex maps")
    for mp in packet.maps:
        need(type(mp) is tuple and len(mp)==len(packet.seed.points)
             and all(type(i) is int and 0<=i<len(pts) for i in mp)
             and len(set(mp))==len(mp), "exact injective copy map")
    rho=(F(2*N-1,2*N),F(0),F(1,2*N),F(0))
    need(packet.seed.field.norm(rho)==ONE, "rotation has exact unit norm")
    expected=(packet.seed.points,tuple(packet.seed.field.mul(rho,z) for z in packet.seed.points))
    need(all(tuple(pts[i] for i in mp)==zs for mp,zs in zip(packet.maps,expected)),
         "same-source actual copy coordinates")
    need(set().union(*(set(mp) for mp in packet.maps))==set(range(len(pts))),
         "no unrelated or missing vertices")
    edge_set(len(pts),packet.skeleton);edge_set(len(pts),packet.full_edges)
    need(type(packet.closing) is tuple and len(packet.closing)==2
         and all(type(i) is int for i in packet.closing), "exact closing edge")
    close=tuple(sorted(mp[packet.seed.tip] for mp in packet.maps))
    need(packet.closing==close and close[0]<close[1], "retain exact closing ports")
    expected_edges=frozenset(tuple(sorted((mp[u],mp[v]))) for mp in packet.maps for u,v in edges)|{close}
    need(packet.skeleton==expected_edges, "complete inherited skeleton")
    need(packet.full_edges==unit_edges(packet.seed.field,pts), "all actual unit edges retained")
    need(packet.skeleton<=packet.full_edges, "skeleton is an actual unit subgraph")
    return packet.full_edges-packet.skeleton


def faithful(packet):
    extra=audit_spindle(packet)
    need(set(packet.maps[0])&set(packet.maps[1])=={packet.maps[0][0]},
         "count product requires otherwise disjoint copies")
    need(not extra, "induced count requires every cross edge")
    return True


def close_count(k,e,d):
    natural(k,2);natural(e);natural(d)
    return k*(k-1)*d*(2*e+(k-2)*d)


def conditional_counts(n,edges,k,ports):
    natural(k,2);edge_set(n,edges)
    need(type(ports) is tuple and len(ports)==2 and all(type(i) is int for i in ports)
         and 0<=ports[0]<n and 0<=ports[1]<n and ports[0]!=ports[1], "exact ports")
    free=tuple(i for i in range(n) if i not in ports)
    counts=[]
    for last in (0,1):
        count=0
        for tail in product(range(k),repeat=len(free)):
            colours=[0]*n;colours[ports[1]]=last
            for i,c in zip(free,tail):colours[i]=c
            count+=all(colours[u]!=colours[v] for u,v in edges)
        counts.append(count)
    return tuple(counts)


def path_seed(length):
    natural(length,2)
    N=length*length
    pts=tuple(scale(j,ONE) for j in range(length+1))
    return Seed(Field(4*N-1),pts,length)


def diamond_chain(m):
    natural(m,1)
    step=add(ONE,A)
    pts=tuple(dict.fromkeys(add(scale(j,step),z)
                           for j in range(m) for z in (ZERO,ONE,A,step)))
    return Seed(Field(12*m*m-1),pts,len(pts)-1)


def chain_counts(k,m):
    natural(k,2);natural(m,1)
    row=(k-1)*(k-2)**2
    transverse=2*(k-2)
    e=(row**m+(k-1)*transverse**m)//k
    d=(row**m-transverse**m)//k
    need(k*e==row**m+(k-1)*transverse**m and k*d==row**m-transverse**m,
         "integer matrix powers")
    return e,d


def chain_four_colouring(packet,m):
    expected=diamond_chain(m)
    need(packet.seed==expected, "same canonical chain source")
    faithful(packet)
    colours=[-1]*len(packet.points)
    for mp in packet.maps:
        for i in range(len(expected.points)):
            c=0 if i%3==0 else i%3
            need(colours[mp[i]] in (-1,c), "shared root colouring compatibility")
            colours[mp[i]]=c
    colours[packet.maps[1][expected.tip]]=3
    colours=tuple(colours)
    proof.verify_colouring(len(packet.points),packet.full_edges,4,colours)
    return colours


def chain_three_refutation(packet,m):
    """Checked structural proof: a diamond's opposite ports agree in 3 colours."""
    need(packet.seed==diamond_chain(m), "same diamond-chain seed")
    audit_spindle(packet)
    for mp in packet.maps:
        for j in range(m):
            u,a,b,v=(mp[3*j+i] for i in range(4))
            required={tuple(sorted(e)) for e in ((u,a),(u,b),(a,b),(a,v),(b,v))}
            need(required<=packet.full_edges, "actual diamond constraints")
    need(packet.closing==tuple(sorted(mp[3*m] for mp in packet.maps)), "tips close")
    # The local equality implications run from the same root to both tips.
    return True


def deleted_three_colour(packet,m,edge):
    """Authenticate an actual three-colouring after deleting any one unit edge."""
    chain_three_refutation(packet,m)
    faithful(packet)
    need(type(edge) is tuple and len(edge)==2 and all(type(i) is int for i in edge)
         and edge in packet.full_edges, "actual deleted unit edge")
    left,right=packet.maps
    blocks=[(left[3*j+3],left[3*j+1],left[3*j+2],left[3*j])
            for j in range(m-1,-1,-1)]
    blocks += [tuple(right[3*j+i] for i in range(4)) for j in range(m)]
    colours=[-1]*len(packet.points)
    state=0
    for u,a,b,v in blocks:
        block_edges={tuple(sorted(e)) for e in ((u,a),(u,b),(a,b),(a,v),(b,v))}
        next_state=1 if edge in block_edges else state
        answer=None
        for ca,cb in product(range(3),repeat=2):
            assignment={u:state,a:ca,b:cb,v:next_state}
            if all(assignment[x]!=assignment[y] for x,y in block_edges-{edge}):
                answer=assignment
                break
        need(answer is not None,"checked local deleted-diamond extension")
        for vertex,colour in answer.items():
            need(colours[vertex] in (-1,colour),"shared source vertices agree")
            colours[vertex]=colour
        state=next_state
    colours=tuple(colours)
    proof.verify_colouring(len(packet.points),packet.full_edges-{edge},3,colours)
    need(colours[edge[0]]==colours[edge[1]],"positive equality witness")
    return colours


def separated_omission(packet,colours,edge,k=6):
    audit_spindle(packet)
    natural(k,4)
    need(type(edge) is tuple and len(edge)==2 and all(type(i) is int for i in edge)
         and edge in packet.full_edges, "actual omitted unit edge")
    proof.verify_colouring(len(packet.points),packet.full_edges,k,colours)
    need(colours[edge[0]]!=colours[edge[1]], "same-source separating witness")
    return packet.full_edges-{edge},colours


def ring_seed():
    pts=tuple(sorted({add(scale(i,ONE),scale(j,A)) for i in range(-2,3)
                      for j in range(-2,3) if max(abs(i),abs(j),abs(i+j))<=2}))
    pts=(ZERO,)+tuple(z for z in pts if z!=ZERO)
    return Seed(Field(11),pts,pts.index(add(ONE,A)))


def rejects(call):
    try:call()
    except (ValueError,TypeError):return True
    return False


def main():
    global CHECKS
    CHECKS=0
    # All local diamond colourings independently authenticate the proof rule.
    edges=frozenset(((0,1),(0,2),(1,2),(1,3),(2,3)))
    local=0
    for c in product(range(3),repeat=4):
        if all(c[u]!=c[v] for u,v in edges):
            check(c[0]==c[3], "diamond opposite ports")
            local+=1
    check(local==6,"nonvacuous local rule")
    # Fixed-assignment algebra: independent literal closing-tip sum.
    for k in range(2,9):
        for e,d in product(range(5),repeat=2):
            literal=k*sum((e if x==0 else d)*(e if y==0 else d)
                          for x in range(k) for y in range(k) if x!=y)
            check(close_count(k,e,d)==literal,"closed two-copy colour sum")
            if k>=3 and e+d>0:
                check((literal==0)==(d==0),"exact k>=3 support criterion")
    paths=[]
    for ell in (2,3):
        packet=build_spindle(path_seed(ell));faithful(packet)
        check(len(packet.points)==2*ell+1 and len(packet.full_edges)==2*ell+1,"induced odd cycle")
        for k in range(2,7):
            e,d=conditional_counts(len(packet.seed.points),audit_seed(packet.seed)[0],k,(0,ell))
            actual=conditional_counts(len(packet.points),packet.full_edges,k,(0,packet.closing[0]))
            total=k*actual[0]+k*(k-1)*actual[1]
            check(total==close_count(k,e,d)==(k-1)**(2*ell+1)-(k-1),"odd-cycle exact counts")
        paths.append({"path_length":ell,"D":packet.seed.field.D,"vertices":len(packet.points),
                      "edges":len(packet.full_edges),"six_colourings":5**(2*ell+1)-5})
    odd_path=path_seed(3)
    check(conditional_counts(4,audit_seed(odd_path)[0],2,(0,3))==(0,1)
          and close_count(2,0,1)==0,"two-colour positive-D exception")
    table=[];omissions=0
    for m in range(1,13):
        packet=build_spindle(diamond_chain(m));faithful(packet)
        check(len(packet.points)==6*m+1 and len(packet.full_edges)==10*m+1,"complete all-pair incidence")
        check(chain_three_refutation(packet,m),"source-bound structural non-three proof")
        colours=chain_four_colouring(packet,m)
        check(max(colours)==3,"explicit four-colouring")
        for edge in packet.full_edges:
            # Incidence and colouring were checked once above; all omitted ports separate.
            check(colours[edge[0]]!=colours[edge[1]],"six-palette omitted edge separates")
            eq=deleted_three_colour(packet,m,edge)
            check(eq[edge[0]]==eq[edge[1]],"three-equality exported for every unit edge")
            omissions+=1
        for k in range(3,9):
            e,d=chain_counts(k,m)
            check((d==0)==(k==3),"no recursive support disappearance")
            check((close_count(k,e,d)==0)==(k==3),"exact chromatic transition")
            row=(k-1)*(k-2)**2;transverse=2*(k-2)
            check(close_count(k,e,d)==(k-1)*(row**(2*m)-transverse**(2*m)),
                  "closed chromatic polynomial")
            if m<=2 and k<=4:
                check((e,d)==conditional_counts(len(packet.seed.points),audit_seed(packet.seed)[0],
                                                k,(0,packet.seed.tip)),"literal seed counts")
        table.append({"m":m,"D":12*m*m-1,"vertices":len(packet.points),
                      "unit_edges":len(packet.full_edges),"chromatic_number":4,
                      "six_colourings":close_count(6,*chain_counts(6,m))})
    packet=build_spindle(diamond_chain(2))
    # Recheck the exact six-colour omitted-edge receipt on the first outside-field chain.
    c=chain_four_colouring(packet,2)
    for e in sorted(packet.full_edges):
        edges2,c2=separated_omission(packet,c,e)
        check(e not in edges2 and c2[e[0]]!=c2[e[1]],"same-source omission packet")
    # Additional unit edges invalidate an induced product-count claim.
    ring=build_spindle(ring_seed())
    extra=audit_spindle(ring)
    check(len(ring.points)==37 and len(ring.full_edges)==90 and len(extra)==5,
          "six geometrical closings, only one marked")
    check(rejects(lambda:faithful(ring)),"do not discard accidental cross edges")
    def seedcolour(z):
        # z=i+j*a, so j=2*s coefficient and i=constant-j/2.
        j=2*z[1];i=z[0]-j/2
        need(i.denominator==j.denominator==1,"lattice source")
        return int(i-j)%3
    badcolour=[-1]*len(ring.points)
    for mp in ring.maps:
        for i,z in enumerate(ring.seed.points):badcolour[mp[i]]=seedcolour(z)
    badcolour[ring.maps[1][ring.seed.tip]]=3
    proof.verify_colouring(len(ring.points),ring.skeleton,4,tuple(badcolour))
    check(any(badcolour[u]==badcolour[v] for u,v in extra),
          "proper skeleton witness fails full induced graph")
    # Joint field/geometry hostile: selected primes split in quadratic subfield.
    split_controls=[]
    for S in ((),(2,),(2,3),(2,3,5),(2,3,5,7),(2,3,5,7,11)):
        P=prod(p for p in S if p!=2);m=2*P;D=12*m*m-1
        check((-D)%8==1,"two splits")
        for p in S:
            if p!=2:check((-D)%p==1,"selected odd prime splits")
        check(isqrt(D)**2!=D and isqrt(3*D)**2!=3*D,"faithful biquadratic field")
        split_controls.append({"selected_primes":S,"m":m,"D":D,
                               "symbolic_vertices":6*m+1,"chromatic_number":4})
    D=3*1237*10193*11329*165713
    check(D==71013535489456671 and D%8==7,"external equation rotation field marker")
    check(pow(-3%11,5,11)==pow(-D%11,5,11)==10,"remaining prime-eleven residue screen")
    bads=(
        lambda:diamond_chain(True),lambda:path_seed(1),lambda:Field(True).audit(),
        lambda:Field(3).audit(),lambda:point((0.0,0,0,0)),
        lambda:close_count(6,True,0),lambda:chain_counts(6,0),
        lambda:audit_spindle(replace(packet,closing=(False,packet.closing[1]))),
        lambda:audit_spindle(replace(packet,full_edges=packet.full_edges-{min(packet.full_edges)})),
        lambda:audit_spindle(replace(packet,maps=((False,)+packet.maps[0][1:],packet.maps[1]))),
        lambda:separated_omission(packet,c,(False,1)),
        lambda:deleted_three_colour(packet,2,(False,1)),
    )
    for f in bads:check(rejects(f),"malformed/source-loss hostile")
    print(json.dumps({"scope":"exact field escape and geometric closure; no six-equality construction",
                      "checks":CHECKS,"path_controls":paths,"diamond_chain_bank":table,
                      "omitted_edges_separated":omissions,
                      "extra_cross_edge_hostile":{"vertices":37,"unit_edges":90,"unrecorded_by_skeleton":5},
                      "finite_prime_escape_hostiles":split_controls,
                      "plane_six_or_seven":"not decided"},indent=2))


if __name__=="__main__":
    main()
