"""Exact cyclotomic contacts and source-bound six-colour extensions.

No floating distances, SAT assumption, or unknown chromatic lower bound.
Run with Python -B [-O] -X utf8; stdout is deterministic JSON.
"""
from dataclasses import dataclass, replace
from fractions import Fraction as F
from functools import lru_cache
from itertools import combinations, product
import json
import sympy as s


def natural(n, lo=0):
    if type(n) is not int or n < lo:
        raise ValueError('exact integer required')
    return n


class Cyclotomic:
    def __init__(self, n):
        self.n = n
        x = s.Symbol('x')
        p = s.Poly(s.cyclotomic_poly(n, x), x)
        self.d = p.degree()
        self.mod = tuple(int(p.nth(i)) for i in range(self.d))
        self.zero = (F(0),)*self.d
        self.one = (F(1),)+(F(0),)*(self.d-1)
        self.roots = tuple(self.reduce([0]*j+[1]) for j in range(n))

    def reduce(self, a):
        a = [F(v) for v in a]
        a += [F(0)]*max(0, self.d-len(a))
        for k in range(len(a)-1, self.d-1, -1):
            c = a[k]
            if c:
                for j, v in enumerate(self.mod):
                    if v:
                        a[k-self.d+j] -= c*v
        return tuple(a[:self.d])

    def scalar(self, a):
        return (F(a),)+(F(0),)*(self.d-1)

    def audit(self, a):
        if type(a) is not tuple or len(a) != self.d:
            raise ValueError('field vector')
        if any(type(v) not in (int, F) for v in a):
            raise ValueError('exact rational field coefficients')
        return tuple(F(v) for v in a)

    def add(self, a, b):
        return tuple(x+y for x, y in zip(a,b))

    def neg(self, a):
        return tuple(-x for x in a)

    def sub(self, a, b):
        return self.add(a,self.neg(b))

    def scale(self, a, t):
        return tuple(t*x for x in a)

    def mul(self, a, b):
        out = [F(0)]*(2*self.d-1)
        aa = tuple((i,x) for i,x in enumerate(a) if x)
        bb = tuple((i,x) for i,x in enumerate(b) if x)
        for i,x in aa:
            for j,y in bb:
                out[i+j] += x*y
        return self.reduce(out)

    def auto(self, a, k):
        out = self.zero
        for i,v in enumerate(a):
            if v:
                out = self.add(out,self.scale(self.roots[i*k % self.n],v))
        return out

    def bar(self,a):
        return self.auto(a,self.n-1)

    def norm(self,a):
        return self.mul(a,self.bar(a))

    def power(self,a,n):
        natural(n)
        out = self.one
        while n:
            if n & 1:
                out = self.mul(out,a)
            a = self.mul(a,a)
            n //= 2
        return out


L = Cyclotomic(55)
C = Cyclotomic(11)
ZETA = L.roots[5]
J = L.zero
SQRT5 = L.zero
for i in range(1,11):
    J = L.add(J,L.scale(L.roots[5*i],1 if pow(i,5,11)==1 else -1))
for i in range(1,5):
    SQRT5 = L.add(SQRT5,L.scale(L.roots[11*i],1 if pow(i,2,5)==1 else -1))
K55 = L.mul(J,SQRT5)


def in_old_field(a):
    """Fixed field of zeta55 -> zeta55^26 is Q(zeta5,sqrt(-11))."""
    a = L.audit(a)
    return L.auto(a,26) == a


def points_tuple(points):
    if type(points) is not tuple or not 1 <= len(points) <= 128:
        raise ValueError('literal point tuple of size 1..128')
    p = tuple(L.audit(x) for x in points)
    if len(set(p)) != len(p):
        raise ValueError('coincident source points')
    return p


@lru_cache(None)
def _unit_edges(points):
    return frozenset((i,j) for i,j in combinations(range(len(points)),2)
                     if L.norm(L.sub(points[i],points[j])) == L.one)


def unit_edges(points):
    return _unit_edges(points_tuple(points))


def pair(e,n):
    if type(e) is not tuple or len(e)!=2:
        raise ValueError('marked edge pair')
    a,b=e; natural(a); natural(b)
    if not a < b < n:
        raise ValueError('ordered distinct source edge')
    return e


def proper(colours, edges, n, q=6):
    natural(q,1)
    if type(colours) is not tuple or len(colours)!=n:
        raise ValueError('complete labelled colouring')
    for c in colours:
        natural(c)
        if c>=q:
            raise ValueError('palette label')
    return all(colours[a]!=colours[b] for a,b in edges)


@dataclass(frozen=True)
class Exterior:
    points: tuple
    old_count: int
    edges: frozenset
    peel: tuple


def compile_exterior(points, old_count):
    p=points_tuple(points); natural(old_count,1)
    if not old_count<=len(p) or len(p)-old_count>17:
        raise ValueError('proved exterior universe: at most seventeen vertices')
    if not all(in_old_field(x) for x in p[:old_count]):
        raise ValueError('old source leaves K')
    if any(in_old_field(x) for x in p[old_count:]):
        raise ValueError('exterior must be outside K')
    es=_unit_edges(p)
    pending=set(range(old_count,len(p))); peel=[]
    while pending:
        possible=[v for v in sorted(pending)
                  if sum((min(v,w),max(v,w)) in es for w in pending if w!=v)<=4]
        if not possible:
            raise ValueError('missing four-degeneracy witness')
        v=possible[0]; pending.remove(v); peel.append(v)
    for v in range(old_count,len(p)):
        if sum((u,v) in es for u in range(old_count))>1:
            raise ValueError('two old contacts violate field contract')
    return Exterior(p,old_count,es,tuple(peel))


def audit_exterior(r):
    if type(r) is not Exterior or type(r.edges) is not frozenset or type(r.peel) is not tuple:
        raise ValueError('exterior receipt type')
    p=points_tuple(r.points)
    natural(r.old_count,1)
    for e in r.edges: pair(e,len(p))
    for i in r.peel: natural(i)
    if r!=compile_exterior(p,r.old_count):
        raise ValueError('source or peel changed')
    return r


def extend_boundary(r, colours):
    audit_exterior(r)
    old_edges=frozenset(e for e in r.edges if e[1]<r.old_count)
    if not proper(colours,old_edges,r.old_count):
        raise ValueError('old six-colouring invalid')
    result=list(colours)+( [None]*(len(r.points)-r.old_count))
    for v in reversed(r.peel):
        forbidden={result[w] for w in range(len(result)) if result[w] is not None
                   and (min(v,w),max(v,w)) in r.edges}
        result[v]=next(c for c in range(6) if c not in forbidden)
    result=tuple(result)
    if not proper(result,r.edges,len(r.points)):
        raise ValueError('extension failed')
    return result


@dataclass(frozen=True)
class Cap:
    source: tuple
    edge: tuple
    exterior: Exterior
    cycle: tuple


def compile_cap(points, edge):
    p=points_tuple(points); pair(edge,len(p))
    if not all(in_old_field(x) for x in p) or edge not in _unit_edges(p):
        raise ValueError('unit edge in old K source required')
    a,b=edge; u=L.sub(p[b],p[a]); partial=L.zero; polys=[]
    for j in range(11):
        polys.append(partial)
        partial=L.add(partial,L.roots[5*j % 55])
    new=tuple(L.add(p[a],L.mul(u,polys[j])) for j in range(2,11))
    ex=compile_exterior(p+new,len(p))
    cycle=(a,b)+tuple(range(len(p),len(p)+9))
    expected=_unit_edges(p)|frozenset(tuple(sorted((cycle[i],cycle[(i+1)%11]))) for i in range(11))
    if ex.edges!=expected:
        raise ValueError('extra or missing cap incidence')
    return Cap(p,edge,ex,cycle)


def audit_cap(r):
    if type(r) is not Cap or type(r.cycle) is not tuple:
        raise ValueError('cap receipt type')
    p=points_tuple(r.source); pair(r.edge,len(p))
    for v in r.cycle: natural(v)
    audit_exterior(r.exterior)
    if r!=compile_cap(p,r.edge):
        raise ValueError('source or cap changed')
    return r


def cap_lift_count(q):
    natural(q,2)
    return ((q-1)**10-1)//q


def contact_certificates():
    """Rational row certificates for the relative E=Q(sqrt(-11)) systems."""
    j=C.zero
    for i in range(1,11):
        j=C.add(j,C.scale(C.roots[i],1 if pow(i,5,11)==1 else -1))
    def vec(a): return s.Matrix([s.Rational(v.numerator,v.denominator) for v in a])
    rows=[]; p=C.zero
    for t in range(11):
        if t>=2:
            b=C.bar(p)
            cols=(C.neg(p),C.neg(C.mul(j,p)),C.neg(b),C.neg(C.mul(j,b)),C.one,j)
            M=s.Matrix.hstack(*(vec(x) for x in cols))
            rhs=vec(C.sub(C.one,C.mul(p,b)))
            jrhs=vec(C.mul(j,tuple(F(x) for x in rhs)))
            rank=M.rank(); augmented=M.row_join(rhs).row_join(jrhs).rank()
            data=dict(vertex=t,rank_over_E=rank//2,augmented_rank_over_E=augmented//2)
            if rank==augmented:
                solution=M.gauss_jordan_solve(rhs)[0]
                if M*solution!=rhs: raise ValueError('contact solution check')
                data['solution']=list(map(str,solution))
            else:
                witness=next(w for w in M.T.nullspace() if (w.T*rhs)[0]!=0)
                witness=witness/(witness.T*rhs)[0]
                if witness.T*M!=s.zeros(1,6) or (witness.T*rhs)[0]!=1:
                    raise ValueError('contact obstruction check')
                data['rational_annihilator']=list(map(str,witness))
            rows.append(data)
        p=C.add(p,C.roots[t])
    return rows


def elliptic_patch():
    alpha=L.scale(L.add(L.one,J),F(1,2)); points=[alpha]
    rows=[]
    for r,y in ((F(1,5),F(8,5)),(F(91,135),F(10552,3645))):
        if y*y!=-275*r**4+150*r*r-3: raise ValueError('quartic source')
        radius=(1+55*r*r)/4
        b=L.scale(L.mul(L.add(L.scalar(radius+2),L.scale(J,y/4)),alpha),F(1,6))
        delta=L.add(L.scalar((55*r*r-1)/(1+55*r*r)),L.scale(K55,2*r/(1+55*r*r)))
        if L.norm(delta)!=L.one or L.norm(b)!=L.scalar(radius) or L.norm(L.sub(alpha,b))!=L.one:
            raise ValueError('mixed contact')
        points.extend((L.mul(L.bar(delta),b),b,L.mul(delta,b)))
        rows.append(dict(r=str(r),y=str(y),radius2=str(radius)))
    return tuple(points),rows


def verify():
    checks=0
    def check(v):
        nonlocal checks
        checks+=1
        if not v: raise AssertionError(checks)
    def hostile(f):
        try: f()
        except ValueError: check(True)
        else: check(False)
    check(L.mul(J,J)==L.scalar(-11)); check(L.mul(SQRT5,SQRT5)==L.scalar(5))
    check(L.mul(K55,K55)==L.scalar(-55)); check(pow(26,5,55)==1 and pow(26,1,55)!=1)
    check(all(in_old_field(a) for a in (J,SQRT5,K55,L.roots[11])))
    check(not in_old_field(ZETA)); check(L.norm(ZETA)==L.one)
    # Formal sixteen-vertex equality core: its five first-layer unit vectors
    # would have Gram matrix of rank five, impossible in the plane.
    pairs=tuple(combinations(range(5),2))
    gs=s.symbols('g0:10'); G=s.eye(5)
    for (i,j),g in zip(pairs,gs): G[i,j]=G[j,i]=g
    eq=[]
    for a,b,c,d in combinations(range(5),4):
        for plus in ((a,b),(a,c),(a,d)):
            v=s.zeros(5,1)
            for i in (a,b,c,d): v[i]=1 if i in plus else -1
            eq.append((v.T*G*v)[0]-1)
    Am,bm=s.linear_eq_to_matrix(eq,gs)
    check(Am.rank()==10 and Am.row_join(bm).rank()==10)
    sol=Am.gauss_jordan_solve(bm)[0]
    check(tuple(sol)==(s.Rational(3,4),)*10)
    forced=G.subs(dict(zip(gs,sol)))
    check(forced.rank()==5 and forced.det()==s.Rational(1,64))
    check(1+6+(6*4+1)//2==19 and (17*5)%2==1)
    # A positive higher-dimensional control: the abstract equality core
    # really has an exact unit realization in R^5. No minimal-dimension
    # claim follows, since the planar second-intersection identity was used.
    units=tuple(tuple(F(3,10)+F(i==j,2) for j in range(5)) for i in range(5))
    r5=((F(0),)*5,)+units+tuple(tuple(units[i][k]+units[j][k] for k in range(5)) for i,j in pairs)
    check(s.Matrix([[sum(a*b for a,b in zip(u,v)) for v in units] for u in units])==forced)
    r5_edges=set()
    for i,j in combinations(range(16),2):
        distance=sum((a-b)**2 for a,b in zip(r5[i],r5[j]))
        expected=(i==0 and 1<=j<=5) or (1<=i<=5 and j>=6 and i-1 in pairs[j-6]) or (i>=6 and j>=6 and not set(pairs[i-6])&set(pairs[j-6]))
        check((distance==1)==expected)
        if distance==1: r5_edges.add((i,j))
    check(len(r5_edges)==40)
    check(all(sum(v in e for e in r5_edges)==5 for v in range(16)))
    check(not any(all(tuple(sorted(e)) in r5_edges for e in combinations(t,2)) for t in combinations(range(16),3)))
    certs=contact_certificates()
    check([r['vertex'] for r in certs if r['rank_over_E']==r['augmented_rank_over_E']]==[2,10])
    check(certs[0]['solution']==['1','0','1','0','1','0'])
    check(certs[-1]['solution']==['0']*6)
    p,quartic=elliptic_patch(); es=unit_edges(p)
    check(len(p)==7 and len(es)==6)
    cap=compile_cap(p,(0,2)); audit_cap(cap)
    check(len(cap.exterior.points)==16 and len(cap.exterior.edges)==16)
    check(cap_lift_count(6)==1627604); check(cap_lift_count(2)==0)
    # Independent labelled path DP, not the closed spectral count.
    for q in range(2,8):
        for start in range(q):
            a=[int(i==start) for i in range(q)]
            for _ in range(10): a=[sum(a)-v for v in a]
            for end in range(q):
                if start!=end: check(a[end]==cap_lift_count(q))
    # Every proper old 3-colouring extends at q=6, retaining the actual labels.
    old_three=0
    for colours in product(range(3),repeat=len(p)):
        if proper(colours,es,len(p),3):
            out=extend_boundary(cap.exterior,colours)
            check(out[:len(p)]==colours and proper(out,cap.exterior.edges,len(out)))
            old_three+=1
    check(old_three==192)
    for shift in range(6):
        six=tuple((c+shift)%6 for c in (0,1,2,3,4,5,0))
        out=extend_boundary(cap.exterior,six)
        check(out[:len(p)]==six and proper(out,cap.exterior.edges,len(out)))
    # Actual 3-colour witness for the odd-cycle attachment.
    old_colour=(0,0,1,0,0,1,0)
    check(proper(old_colour,es,len(p),3))
    colour=list(old_colour)+[None]*9
    cycle=cap.cycle
    for idx in range(2,11):
        v=cycle[idx]; forbidden={colour[cycle[idx-1]]}
        if idx==10: forbidden.add(colour[cycle[0]])
        colour[v]=next(c for c in range(3) if c not in forbidden)
    check(proper(tuple(colour),cap.exterior.edges,len(colour),3))
    # Eleven same-centre rotations of a nonzero marked K segment: disjoint edges.
    a=L.scalar(2); b=L.scalar(3)
    fan=tuple(L.mul(L.roots[5*i%55],x) for i in range(11) for x in (a,b))
    fan_edges=unit_edges(fan)
    check(fan_edges==frozenset((2*i,2*i+1) for i in range(11)))
    # The exceptional central point is retained: rotations of {0,1} give a star.
    star=(L.zero,)+tuple(L.roots[5*i] for i in range(11))
    check(unit_edges(star)==frozenset((0,i) for i in range(1,12)))
    # Outside-K is essential: one new OLD-field centre can block a fixed
    # six-colour boundary, although its entire star is two-colourable.
    gamma=L.scale(L.add(L.scalar(5),J),F(1,6))
    leaves=tuple(L.power(gamma,i) for i in range(6))
    check(unit_edges(leaves)==frozenset())
    old_star=leaves+(L.zero,)
    check(unit_edges(old_star)==frozenset((i,6) for i in range(6)))
    check(proper(tuple(range(6)),frozenset(),6))
    check(all(not proper(tuple(range(6))+(c,),unit_edges(old_star),7) for c in range(6)))
    hostile(lambda:compile_exterior(old_star,6))
    hostile(lambda:compile_exterior((L.zero,)+tuple(L.roots[i] for i in range(18)),1))
    hostile(lambda:compile_cap(p,(False,2)))
    hostile(lambda:compile_cap(p,(0,1)))
    hostile(lambda:compile_cap((p[0],p[0]),(0,1)))
    hostile(lambda:compile_exterior(cap.exterior.points,True))
    hostile(lambda:compile_exterior(cap.exterior.points,6))
    hostile(lambda:extend_boundary(cap.exterior,(False,)+old_colour[1:]))
    hostile(lambda:extend_boundary(cap.exterior,(0,)*7))
    hostile(lambda:audit_exterior(replace(cap.exterior,peel=cap.exterior.peel[::-1])))
    hostile(lambda:audit_cap(replace(cap,cycle=(False,)+cap.cycle[1:])))
    hostile(lambda:cap_lift_count(6.0))
    hostile(lambda:points_tuple((tuple([0.0]+[0]*39),)))
    return dict(status='PROVED scoped interfaces; FINITE-EXACT controls; six-colour forcing OPEN',
                checks=checks,fields=dict(old='Q(zeta5,sqrt(-11))',new='Q(zeta55)',relative_degree=5),
                contact_certificates=certs,quartic_points=quartic,
                source_edges=sorted(es),cap_edges=sorted(cap.exterior.edges),cap_cycle=cap.cycle,
                cap_coordinates=[[str(x) for x in a] for a in cap.exterior.points],
                cap_three_colouring=colour,old_three_colourings=old_three,
                six_lifts_per_old_colouring=cap_lift_count(6),
                old_six_count=6*5**6,cap_six_count=5**6*(5**10-1),
                same_centre=dict(segment_copies=11,vertices=len(fan),edges=len(fan_edges)),
                sixteen_core=dict(disjoint_pair_equations=15,forced_gram_rank=5,forced_gram_determinant='1/64',
                                  positive_realization_dimension=5,positive_coordinates=[[str(x) for x in p] for p in r5],
                                  positive_edges=sorted(r5_edges),minimal_realization_dimension='not claimed'),
                exterior_bound=17,next_obligation='At least18 coupled exterior vertices; genuine support loss, not only new contacts or multiplicities.')


if __name__=='__main__':
    print(json.dumps(verify(),indent=2,sort_keys=True))
