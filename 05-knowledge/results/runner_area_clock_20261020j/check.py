"""Exact polygon clocks and source-bound piecewise-linear loneliness controls."""
from dataclasses import dataclass
from fractions import Fraction as F
from itertools import combinations, product
from math import isqrt, gcd
from pathlib import Path
import hashlib
import json
import sys

HERE = Path(__file__).resolve().parent
CHECKS = 0

def need(ok, msg):
    global CHECKS
    CHECKS += 1
    if not ok: raise ValueError(msg)

def rational(x):
    if type(x) not in (int, F): raise ValueError('exact rational required')
    return F(x)

def dist(x, y=0):
    q = (rational(x)-rational(y)) % 1
    return min(q, 1-q)

@dataclass(frozen=True)
class Lift:
    knots: tuple

    def __post_init__(self):
        if type(self.knots) is not tuple or len(self.knots)<2:
            raise ValueError('immutable clock knots required')
        for p in self.knots:
            if type(p) is not tuple or len(p)!=2: raise ValueError('clock pair')
            for x in p: rational(x)
        if self.knots[0]!=(0,0) or self.knots[-1][0]!=1:
            raise ValueError('marked period endpoints')
        if rational(self.knots[-1][1]).denominator!=1:
            raise ValueError('integral winding')
        if self.knots[-1][1]<0: raise ValueError('nonnegative winding')
        for (x,y),(u,v) in zip(self.knots,self.knots[1:]):
            if not x<u or not y<=v: raise ValueError('monotone clock')
        if self.winding and any(y==v for (_,y),(_,v) in zip(self.knots,self.knots[1:])):
            raise ValueError('positive winding requires strict clock')

    @property
    def winding(self): return int(self.knots[-1][1])

    def __call__(self, t):
        t=rational(t); lap=t//1; u=t-lap
        for (x,y),(a,b) in zip(self.knots,self.knots[1:]):
            if x<=u<=a:
                return lap*self.winding+rational(y)+(b-y)*(u-x)/(a-x)
        raise RuntimeError('uncovered clock')

    def inverse(self, y):
        if self.winding!=1: raise ValueError('inverse requires degree one')
        return Lift(tuple((rational(b),rational(a)) for a,b in self.knots))(y)

def linear(v):
    if type(v) is not int or v<0: raise ValueError('exact nonnegative speed')
    return Lift(((F(0),F(0)),(F(1),F(v))))

def compose(f, clock):
    if type(f) is not Lift or type(clock) is not Lift or clock.winding!=1:
        raise ValueError('one shared degree-one clock required')
    xs={F(x) for x,_ in clock.knots}
    xs.update(clock.inverse(x) for x,_ in f.knots)
    return Lift(tuple((x,f(clock(x))) for x in sorted(xs)))

def wrapped_double(g):
    if type(g) is not Lift or g.winding!=1: raise ValueError('marked spatial circle lift')
    xs={F(x+j,2) for j in (0,1) for x,_ in g.knots}
    return Lift(tuple((x,g(2*x)) for x in sorted(xs)))

def maximum_clearance(phases):
    """Exact lower-envelope maximum for the supplied literal PL histories."""
    if type(phases) is not tuple or not phases or any(type(f) is not Lift for f in phases):
        raise ValueError('literal finite phase tuple')
    breaks={F(0),F(1)}
    for f in phases:
        for (x,y),(u,v) in zip(f.knots,f.knots[1:]):
            breaks.update((F(x),F(u)))
            if y==v: continue
            for k in range(int(2*y)-1,int(2*v)+2):
                h=F(k,2)
                if y<h<v: breaks.add(F(x)+(h-y)*(u-x)/(v-y))
    points=set(breaks)
    for x,y in zip(sorted(breaks),sorted(breaks)[1:]):
        lines=[]
        for f in phases:
            a,b=dist(f(x)),dist(f(y)); slope=(b-a)/(y-x)
            lines.append((slope,a-slope*x))
        for (a,b),(c,d) in combinations(lines,2):
            if a!=c and x<(z:=(d-b)/(a-c))<y: points.add(z)
    values={x:min(dist(f(x)) for f in phases) for x in points}
    best=max(values.values())
    return {'maximum':best,'witnesses':tuple(x for x in sorted(points) if values[x]==best),
            'breaks':tuple(sorted(breaks)),'candidates':len(points)}

def sqrt_rat(x):
    x=rational(x); a=isqrt(x.numerator); b=isqrt(x.denominator)
    if a*a!=x.numerator or b*b!=x.denominator: raise ValueError('nonrational length in this exact reader')
    return F(a,b)

def sub(a,b): return tuple(x-y for x,y in zip(a,b))
def cross(a,b): return a[0]*b[1]-a[1]*b[0]
def dot(a,b,D): return a[0]*b[0]+D*a[1]*b[1]

@dataclass(frozen=True)
class Polygon:
    vertices: tuple
    D: int=1

    def __post_init__(self):
        if type(self.D) is not int or self.D<1: raise ValueError('metric field parameter')
        if type(self.vertices) is not tuple or len(self.vertices)<3: raise ValueError('literal polygon')
        for p in self.vertices:
            if type(p) is not tuple or len(p)!=2: raise ValueError('point dimension')
            for x in p: rational(x)
        if len(set(self.vertices))!=len(self.vertices): raise ValueError('distinct vertices')
        for i,a in enumerate(self.vertices):
            b=self.vertices[(i+1)%len(self.vertices)]; e=sub(b,a)
            if cross(a,b)<=0: raise ValueError('origin must be strictly interior')
            for j,c in enumerate(self.vertices):
                if j not in (i,(i+1)%len(self.vertices)) and cross(e,sub(c,a))<=0:
                    raise ValueError('strict convex orientation')
            sqrt_rat(dot(e,e,self.D))

    def data(self):
        lengths=[]; areas=[]
        for i,a in enumerate(self.vertices):
            b=self.vertices[(i+1)%len(self.vertices)]; e=sub(b,a)
            lengths.append(sqrt_rat(dot(e,e,self.D))); areas.append(cross(a,b))
        return tuple(lengths),tuple(areas)

    def area_clock(self):
        lengths,areas=self.data(); P=sum(lengths); A=sum(areas)
        x=y=F(0); knots=[(x,y)]
        for l,a in zip(lengths,areas):
            x+=l/P;y+=a/A;knots.append((x,y))
        return Lift(tuple(knots))

    def point(self,s):
        s=rational(s)%1; lengths,_=self.data(); x=F(0);P=sum(lengths)
        for i,l in enumerate(lengths):
            y=x+l/P
            if x<=s<=y:
                f=(s-x)/(y-x);a=self.vertices[i];b=self.vertices[(i+1)%len(self.vertices)]
                return tuple(a[j]+f*(b[j]-a[j]) for j in (0,1))
            x=y
        raise RuntimeError('polygon interval')

def angular_safe(a,b,D,n):
    """Exact angle threshold for n=2,3,4,6, no floating angles."""
    if type(n) is not int or n not in (2,3,4,6): raise ValueError('declared exact cosine alphabet')
    if type(D) is not int or D<1: raise ValueError('positive exact metric parameter')
    for p in (a,b):
        if type(p) is not tuple or len(p)!=2: raise ValueError('exact planar point')
        for x in p: rational(x)
        if p==(0,0): raise ValueError('angular origin has no direction')
    z=dot(a,b,D); p=dot(a,a,D)*dot(b,b,D)
    if n==2: return z<0 and z*z==p
    if n==3: return z<0 and 4*z*z>=p
    if n==4: return z<=0
    return z<=0 or 4*z*z<=p

def encode(x):
    if type(x) is F:return [x.numerator,x.denominator]
    if type(x) is dict:return {str(k):encode(v) for k,v in x.items()}
    if type(x) in (tuple,list):return [encode(v) for v in x]
    return x

def verify():
    global CHECKS
    CHECKS=0
    square=Polygon(((F(1),F(1)),(F(-1),F(1)),(F(-1),F(-1)),(F(1),F(-1))))
    triangle=Polygon(((F(-1),F(-1,3)),(F(1),F(-1,3)),(F(0),F(2,3))),3)
    hexagon=Polygon(((F(1),F(0)),(F(1,2),F(1,2)),(F(-1,2),F(1,2)),
                     (F(-1),F(0)),(F(-1,2),F(-1,2)),(F(1,2),F(-1,2))),3)
    rectangle=Polygon(((F(2),F(1)),(F(-2),F(1)),(F(-2),F(-1)),(F(2),F(-1))))
    grid=tuple(F(j,60) for j in range(60)); records=[]
    for m,p in ((3,triangle),(4,square),(6,hexagon)):
        area=p.area_clock()
        need(all(x==y for x,y in area.knots),'regular polygon equal-area arclength')
        lengths,areas=p.data()
        need(len(set(p.D*a*a/(l*l) for a,l in zip(areas,lengths)))==1,'common squared side distance')
        for n in (2,3,4,6):
            misses=[]
            for s,t in product(grid,repeat=2):
                same=(dist(s,t)>=F(1,n))==angular_safe(p.point(s),p.point(t),p.D,n)
                if m%n==0:need(same,'matched-count critical predicate')
                elif not same and not misses:misses.append((s,t))
            if m%n:need(bool(misses),'mismatched count has exact hostile')
            records.append({'sides':m,'runners':n,'compatible':m%n==0,'first_hostile':misses})
    a=rectangle.area_clock();slopes=tuple((v-y)/(u-x) for (x,y),(u,v) in zip(a.knots,a.knots[1:]))
    need(slopes==(F(3,4),F(3,2),F(3,4),F(3,2)),'non-tangential rectangle density')
    for s,t in product(grid,repeat=2):
        need(min(slopes)*dist(s,t)<=dist(a(s),a(t))<=max(slopes)*dist(s,t),'two-arc distortion bounds')
    need(a(F(1,3))==F(1,4),'area and perimeter differ')
    # The same target angle need not mean the same physical chord.
    corner=dot(sub(square.point(0),square.point(F(1,4))),sub(square.point(0),square.point(F(1,4))),1)
    middle=dot(sub(square.point(F(1,8)),square.point(F(3,8))),sub(square.point(F(1,8)),square.point(F(3,8))),1)
    need((corner,middle)==(4,2),'radial chord versus physical chord')
    g=Lift(((F(0),F(0)),(F(1,3),F(1,6)),(F(2,3),F(5,6)),(F(1),F(1))))
    altered=wrapped_double(g)
    native=maximum_clearance((linear(1),linear(2)))
    hostile=maximum_clearance((linear(1),altered))
    need(native['maximum']==F(1,3) and hostile['maximum']==F(3,10),'individual-clock global hostile')
    need(hostile['witnesses']==(F(3,10),F(7,10)),'complete hostile maximizers')
    for t in grid:
        need(altered(t+1)==altered(t)+2,'same winding and average speed')
    common=[]
    for speeds in ((1,2),(1,3),(1,2,3),(1,2,3,4),(2,5,7),(1,3,8,11)):
        phases=tuple(linear(v) for v in speeds)
        before=maximum_clearance(phases);after=maximum_clearance(tuple(compose(f,g) for f in phases))
        need(before['maximum']==after['maximum'],'shared clock preserves exact maximum')
        need(all(min(dist(f(g(t))) for f in phases)==after['maximum'] for t in after['witnesses']), 'literal shared-clock witnesses')
        common.append({'speeds':speeds,'maximum':before['maximum'],'new_witnesses':after['witnesses']})
    ap=[]
    for n in range(2,17):
        result=maximum_clearance(tuple(linear(v) for v in range(1,n)))
        expected=tuple(F(j,n) for j in range(1,n) if gcd(j,n)==1)
        need(result['maximum']==F(1,n) and result['witnesses']==expected,'AP marked primitive clocks')
        ap.append({'runners':n,'maximum':result['maximum'],'witnesses':result['witnesses']})
    # Same honest half-turn relation, different loneliness. No tie is invented.
    configs=((F(0),F(1,10),F(1,5)),(F(0),F(7,20),F(2,5)))
    def directions(xs):
        return tuple(0 if (xs[j]-xs[i])%1 in (0,F(1,2)) else
                     (1 if (xs[j]-xs[i])%1<F(1,2) else -1) for i,j in combinations(range(len(xs)),2))
    need(directions(configs[0])==directions(configs[1])==(1,1,1),'orientation loses target clearance')
    need(min(dist(x) for x in configs[0][1:])<F(1,3)<=min(dist(x) for x in configs[1][1:]),'same orientation different loneliness')
    need(directions((F(0),F(1,2)))==(0,),'antipodal tie retained')
    cycles=[]; unseen=set(range(1,7))
    while unseen:
        x=min(unseen); cycle=[]
        while x not in cycle:
            cycle.append(x); unseen.remove(x); x=2*x%7
        cycles.append(tuple(cycle))
    need(cycles==[(1,2,4),(3,6,5)],'unit permutation is not necessarily one cycle')
    def rejects(fn):
        try:fn()
        except (ValueError,TypeError):need(True,'rejected malformed input');return
        need(False,'type or source hostile accepted')
    rejects(lambda:linear(True));rejects(lambda:dist(.5));rejects(lambda:angular_safe((1,0),(0,1),1,True))
    rejects(lambda:angular_safe((1.,0),(0,1),1,4))
    rejects(lambda:angular_safe((1,0),(0,0),1,4))
    rejects(lambda:angular_safe((1,0),(0,1),True,4))
    rejects(lambda:Lift(((0,0),(F(1,2),F(1,2)),(F(1,2),F(3,4)),(1,1))))
    rejects(lambda:Polygon(((1,1),(-1,1),(-1,-1),(1,-1)),True))
    rejects(lambda:Polygon(((1,1,0),(-1,1),(-1,-1),(1,-1))))
    return {'status':'PROVED clock and threshold reductions; FINITE-EXACT controls; general LRC OPEN',
            'checks':CHECKS,'grid_denominator':60,'regular_polygon_controls':records,
            'rectangle':{'area_knots':a.knots,'density_bounds':(min(slopes),max(slopes))},
            'chord_squared':{'corner':corner,'side_midpoints':middle},
            'independent_clock':{'G':g.knots,'second_phase':altered.knots,'native':native,'changed':hostile},
            'common_clock_controls':common,'AP_controls':ap,
            'orientation_hostiles':configs,'doubling_mod7_cycles':cycles}

if __name__=='__main__':
    data=json.dumps(encode(verify()),sort_keys=True,indent=2)+'\n';p=HERE/'receipt.json'
    if '--write' in sys.argv:p.write_text(data,encoding='utf-8',newline='\n')
    else:need(p.read_text(encoding='utf-8')==data,'saved receipt mismatch')
    print(json.dumps({'checks':json.loads(data)['checks'],'LF_SHA256':hashlib.sha256(data.encode()).hexdigest()}))
