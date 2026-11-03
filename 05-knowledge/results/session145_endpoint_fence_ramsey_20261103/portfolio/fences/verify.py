#!/usr/bin/env python3
"""Exact geometry and rational-interval controls; no optimization-removable asserts."""
from dataclasses import dataclass
from fractions import Fraction as F
from functools import cmp_to_key
from itertools import combinations
from pathlib import Path
import hashlib, json, sys

CHECKS = 0
def check(ok, label):
    global CHECKS
    CHECKS += 1
    if not ok:
        raise RuntimeError(label)

@dataclass(frozen=True)
class Q:
    a: F = F(0)
    b: F = F(0)
    def __post_init__(self):
        object.__setattr__(self, 'a', F(self.a))
        object.__setattr__(self, 'b', F(self.b))
    def __add__(self, other):
        z = coerce(other)
        return Q(self.a+z.a, self.b+z.b)
    __radd__ = __add__
    def __neg__(self): return Q(-self.a,-self.b)
    def __sub__(self,other): return self+-coerce(other)
    def __rsub__(self,other): return coerce(other)+-self
    def __mul__(self,other):
        z=coerce(other)
        return Q(self.a*z.a+13*self.b*z.b,self.a*z.b+self.b*z.a)
    __rmul__=__mul__
    def __truediv__(self,other):
        z=coerce(other); den=z.a*z.a-13*z.b*z.b
        if not den: raise ZeroDivisionError
        return self*Q(z.a/den,-z.b/den)
    def __pow__(self,n):
        if n<0: return (Q(1)/self)**(-n)
        out=Q(1)
        while n:
            if n&1: out=out*self
            self=self*self; n//=2
        return out
    def sign(self):
        sa=(self.a>0)-(self.a<0); sb=(self.b>0)-(self.b<0)
        if not sa: return sb
        if not sb or sa==sb: return sa
        c=self.a*self.a-13*self.b*self.b
        return sa*((c>0)-(c<0))
    def __lt__(self,other): return (self-other).sign()<0
    def text(self): return {'rational':str(self.a),'sqrt13':str(self.b)}

def coerce(x): return x if isinstance(x,Q) else Q(x)
def point(x,y): return (coerce(x),coerce(y))
def sub(p,q): return (p[0]-q[0],p[1]-q[1])
def cross(v,w): return v[0]*w[1]-v[1]*w[0]
def orient(a,b,c): return cross(sub(b,a),sub(c,a)).sign()
def between(p,a,b):
    if orient(a,b,p): return False
    return ((p[0]-a[0])*(p[0]-b[0])+(p[1]-a[1])*(p[1]-b[1])).sign()<=0

def intersection(a,b,c,d):
    o1,o2,o3,o4=orient(a,b,c),orient(a,b,d),orient(c,d,a),orient(c,d,b)
    if o1*o2<0 and o3*o4<0: return 'cross'
    contact=set()
    for p in (a,b):
        if between(p,c,d): contact.add(p)
    for p in (c,d):
        if between(p,a,b): contact.add(p)
    if len(contact)>1: return 'overlap'
    if not contact: return 'none'
    p=next(iter(contact))
    return 'endpoint' if p in (a,b) and p in (c,d) else 'T'

def graph_info(n,edges):
    adj=[set() for _ in range(n)]
    for a,b in edges: adj[a].add(b); adj[b].add(a)
    seen=set(); comp=[]
    for v in range(n):
        if v in seen: continue
        stack=[v]; seen.add(v); vertices=[]
        while stack:
            w=stack.pop(); vertices.append(w)
            for z in adj[w]:
                if z not in seen: seen.add(z); stack.append(z)
        e=sum(len(adj[w]) for w in vertices)//2
        comp.append((len(vertices),e))
    # Independent cactus criterion by enumerating all simple cycles once.
    cycles=set()
    for s in range(n):
        def walk(path):
            v=path[-1]
            for w in adj[v]:
                if w==s and len(path)>=3:
                    es=frozenset(tuple(sorted((path[i],path[(i+1)%len(path)]))) for i in range(len(path)))
                    cycles.add(es)
                elif w>s and w not in path: walk(path+[w])
        walk([s])
    counts={e:0 for e in edges}
    for cycle in cycles:
        for e in cycle: counts[e]+=1
    return {'vertices':n,'edges':len(edges),'components':len(comp),
            'cycle_rank':len(edges)-n+len(comp),
            'pseudoforest':all(e<=v for v,e in comp),
            'cactus':all(c<=1 for c in counts.values()),
            'minimum_degree':min((len(a) for a in adj),default=0),
            'cycles':sorted(map(len,cycles))},adj

def geometric_fixture(segments, metric_y2, mode):
    types=[]
    for i,(a,b) in enumerate(segments):
        dx,dy=sub(b,a)
        check(dx*dx+metric_y2*dy*dy==Q(1),'original segment unit length')
        for c,d in segments[:i]:
            kind=intersection(a,b,c,d); types.append(kind)
            check(kind not in ('cross','overlap'),'segments cannot cross or overlap')
            if mode=='product': check(kind!='T','product forbids T junctions')
    if mode=='fence':
        for i,(a,b) in enumerate(segments):
            for p in (a,b):
                check(any(j!=i and between(p,c,d) for j,(c,d) in enumerate(segments)), 'fence endpoint support')
    points=set(p for e in segments for p in e); vertex=list(sorted(points))
    ids={p:i for i,p in enumerate(vertex)}; edges=set()
    for a,b in segments:
        on=[p for p in vertex if between(p,a,b)]
        # x-order suffices unless vertical; either orientation gives the same undirected edges.
        on.sort()
        for c,d in zip(on,on[1:]): edges.add(tuple(sorted((ids[c],ids[d]))))
    info,adj=graph_info(len(vertex),edges)
    def polar(a,b):
        ha=0 if a[1].sign()>0 or (a[1].sign()==0 and a[0].sign()>0) else 1
        hb=0 if b[1].sign()>0 or (b[1].sign()==0 and b[0].sign()>0) else 1
        if ha!=hb: return ha-hb
        return -cross(a,b).sign()
    rotations=[]
    for v,ns in enumerate(adj):
        rotations.append(sorted(ns,key=cmp_to_key(lambda a,b: polar(sub(vertex[a],vertex[v]),sub(vertex[b],vertex[v])))))
    visited=set(); areas=[]
    for v,w in [(a,b) for a,b in edges]+[(b,a) for a,b in edges]:
        if (v,w) in visited: continue
        start=(v,w); cur=start; area2=Q(0)
        while cur not in visited:
            visited.add(cur); a,b=cur; area2+=cross(vertex[a],vertex[b])
            around=rotations[b]; c=around[(around.index(a)-1)%len(around)]; cur=(b,c)
        check(cur==start,'half-edge walk closes at its own origin')
        if area2.sign()>0: areas.append(area2/2)
    check(len(areas)==info['cycle_rank'],'connected/non-nested fixture face count')
    info['original_fences']=len(segments)
    info['T_contacts']=types.count('T')
    info['squared_face_areas']=[(metric_y2*a*a).text() for a in sorted(areas)]
    return info,areas

def boundary(points): return list(zip(points,points[1:]+points[:1]))

def geometry_tests():
    r=Q(0,1); dy2=(r-1)/6; a=(r+1)/4; b=Q(F(3,2)); h=(r+7)/4
    top=[point(-1,0),point(-a,b),point(F(-1,2),h),point(F(1,2),h),point(a,b),point(1,0)]
    bottom=[point(x,-y) for x,y in top]
    seg=boundary(top)[:-1]+boundary(bottom)[:-1]+[(point(-1,0),point(0,0)),(point(0,0),point(1,0))]
    info,areas=geometric_fixture(seg,dy2,'product')
    check(info['vertices']==11 and info['edges']==12 and info['cycle_rank']==2,'theta topology')
    check(not info['cactus'] and not info['pseudoforest'],'theta not cactus')
    check(info['cycles']==[7,7,10],'theta cycle lengths 5+2,5+2,5+5')
    check(len(areas)==2 and areas[0]==areas[1],'equal product fields')
    product=dy2*areas[0]*areas[1]
    check(product==(89+26*r)/16,'Friedman exact product identity')
    check((product-Q(1)).sign()>0,'product example does not obey fence cap')
    check(F(18,5)**2<13 and F(7,4)**2>3,'comparison radical bounds')
    check(26*F(18,5)-7-48*F(7,4)>0,'strict cactus gap')
    # Strict convexity of the physical six-corner polygon.
    check(all(orient(top[i-1],top[i],top[(i+1)%6])==-1 for i in range(6)), 'upper hexagon convex')
    # Fence counterexample to pseudoforest extension: adjacent unit squares.
    sq=[(point(x,y),point(x+1,y)) for y in (0,1) for x in (0,1)]
    sq += [(point(x,0),point(x,1)) for x in (0,1,2)]
    sqinfo,sqareas=geometric_fixture(sq,Q(1),'fence')
    check(sum(sqareas,Q())==Q(2) and len(sq)==7,'seven-fence area two')
    check((sum(sqareas,Q())-Q(F(7,4))).sign()>0 and not sqinfo['pseudoforest'],'fence obstruction needs graph scope')
    # Exact full-fence cactus with T junctions; y=coefficient*sqrt(35).
    left=[point(0,0),point(F(1,6),F(1,6)),point(F(1,3),0)]
    right=[point(F(2,3),0),point(F(5,6),F(1,6)),point(1,0)]
    cactus=boundary(left)[:-1]+boundary(right)[:-1]+[(point(0,0),point(1,0))]
    cinfo,careas=geometric_fixture(cactus,Q(35),'fence')
    check(cinfo['cactus'] and not cinfo['pseudoforest'],'T-bridge cactus topology')
    check(cinfo['original_fences']==5 and cinfo['edges']==7 and cinfo['T_contacts']==2,'split edge count differs')
    check(cinfo['cycles']==[3,3],'two triangular cycle blocks')
    check(all(35*z*z==Q(F(35,1296)) for z in careas),'T-bridge fields area sqrt35/36')
    # Deliberate invalid product transfer must be rejected.
    caught=False
    try: geometric_fixture(cactus,Q(35),'product')
    except RuntimeError as e: caught=str(e)=='product forbids T junctions'
    check(caught,'T-fence must fail endpoint-only product validation')
    # Positive no-through pseudoforest control.
    single=boundary([point(0,0),point(1,0),point(1,1),point(0,1)])
    sinfo,sareas=geometric_fixture(single,Q(1),'fence')
    check(sinfo['pseudoforest'] and sareas==[Q(1)],'square equality control')
    check(intersection(point(0,0),point(1,1),point(0,1),point(1,0))=='cross','proper crossing hostile')
    check(intersection(point(0,0),point(2,0),point(1,0),point(3,0))=='overlap','overlap hostile')
    return {'product_theta12':info,'adjacent_squares':sqinfo,'T_bridge_cactus':cinfo,'square':sinfo}

# Exact rational interval evaluation; Machin's identity and alternating series.
def atan_bounds(q,terms):
    total=F(0)
    for j in range(terms): total+=(-1)**j*q**(2*j+1)/(2*j+1)
    next_term=(-1)**terms*q**(2*terms+1)/(2*terms+1)
    return min(total,total+next_term),max(total,total+next_term)

def trig_bounds(x,cosine=False,terms=24):
    term=F(1) if cosine else x; total=term
    for j in range(1,terms):
        a=2*j-1 if cosine else 2*j
        term=-term*x*x/(a*(a+1)); total+=term
    a=2*terms-1 if cosine else 2*terms
    nxt=-term*x*x/(a*(a+1))
    return min(total,total+nxt),max(total,total+nxt)

def round_interval(lo,hi,digits=48):
    scale=10**digits
    l=lo*scale; h=hi*scale
    return F(l.numerator//l.denominator,scale),F(-((-h.numerator)//h.denominator),scale)

def polygon_intervals(n):
    a,b=atan_bounds(F(1,5),40); c,d=atan_bounds(F(1,239),12)
    pl,pu=round_interval(16*a-4*d,16*b-4*c)
    check(F(314159,100000)<pl<pu<F(314160,100000),'Machin pi interval sanity')
    result={}
    for k in range(3,n+1):
        xl,xu=pl/k,pu/k
        sl=trig_bounds(xl)[0]; su=trig_bounds(xu)[1]
        cl=trig_bounds(xu,True)[0]; cu=trig_bounds(xl,True)[1]
        result[k]=round_interval(k*cl/(4*su),k*cu/(4*sl))
        check(0<result[k][0]<result[k][1],'positive polygon interval')
        check(result[k][1]-result[k][0]<F(1,10**40),'narrow polygon interval')
    return result

def product_interval(part,a):
    lo=hi=F(1)
    for k in part: lo*=a[k][0]; hi*=a[k][1]
    return lo,hi

def parts(n,minimum=3):
    if n==0:
        yield ()
    for k in range(minimum,n+1):
        for tail in parts(n-k,k): yield (k,)+tail

def balanced(n,r):
    q,s=divmod(n,r)
    return (q,)*(r-s)+(q+1,)*s

def interval_decimal(x,places=12):
    # Human-readable outward decimal interval without binary floating point.
    scale=10**places; lo=x[0]*scale; hi=x[1]*scale
    flo=lo.numerator//lo.denominator; cei=-((-hi.numerator)//hi.denominator)
    def fmt(i): return str(i//scale)+'.'+str(i%scale).zfill(places)
    return [fmt(flo),fmt(cei)]

def partition_tests():
    a=polygon_intervals(50); table=[]; count=0
    for n in range(3,51):
        candidates=[balanced(n,r) for r in range(1,n//3+1)]
        best=max(candidates,key=lambda p:product_interval(p,a)[0])
        bi=product_interval(best,a)
        for p in candidates:
            if p!=best: check(product_interval(p,a)[1]<bi[0],'separated balanced maximum')
        if n<=30:
            for p in parts(n):
                count+=1
                if p!=best: check(product_interval(p,a)[1]<bi[0],'all-partition maximum agrees')
        table.append({'n':n,'maximizing_cycle_lengths':best,'maximum_interval':interval_decimal(bi)})
    check(table[9]['maximizing_cycle_lengths']==(12,), 'twelve-edge cactus maximum is dodecagon')
    check(table[13]['maximizing_cycle_lengths']==(8,8), 'sixteen-edge balancing transition')
    return {'universe':'all integer partitions with parts >=3 for 3<=n<=30; all balanced candidates through n=50',
            'partitions_tested':count,'table':table}

def graph_reduction_tests():
    total=hits=0
    for n in range(1,6):
        possible=list(combinations(range(n),2))
        for mask in range(1<<len(possible)):
            edges={e for i,e in enumerate(possible) if mask>>i&1}
            info,adj=graph_info(n,edges); total+=1
            if info['pseudoforest'] and info['minimum_degree']>=2:
                hits+=1
                check(all(len(a)==2 for a in adj),'minimum-degree-two pseudoforest is cycle union')
                check(info['edges']==n and info['cycle_rank']==info['components'],'cycle union Euler count')
            if info['cactus']:
                check(sum(info['cycles'])<=len(edges),'cactus cycle edge budget')
                check(len(info['cycles'])==info['cycle_rank'],'cactus cycle count')
    return {'universe':'all labelled simple graphs of orders 1 through 5','graphs':total,'minimum_degree_two_pseudoforests':hits}

def main():
    result={'status':'PASS','geometry':geometry_tests(),'partitions':partition_tests(),'graph_reduction':graph_reduction_tests()}
    result['checked_conditions']=CHECKS
    result['source_sha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    out=Path(sys.argv[1] if len(sys.argv)>1 else 'verification.json')
    out.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    print(json.dumps({'status':'PASS','checked_conditions':CHECKS,'receipt':str(out),'graphs':result['graph_reduction']['graphs'],'partitions':result['partitions']['partitions_tested']}))
if __name__=='__main__': main()
