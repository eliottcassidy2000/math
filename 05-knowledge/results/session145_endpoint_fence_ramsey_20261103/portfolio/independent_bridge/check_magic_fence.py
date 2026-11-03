#!/usr/bin/env python3
"""Independent audit, importing no Session145 discovery modules.
Exact SymPy linear algebra/polynomials/geometry and integer enumeration.
"""
from itertools import permutations, product, combinations
from fractions import Fraction
from pathlib import Path
import json, sys, hashlib, subprocess
import sympy as s

CHECKS=0
def must(ok,label):
    global CHECKS
    CHECKS+=1
    if not bool(ok): raise RuntimeError(label)

def magic_tests():
    lines=((0,1,2),(3,4,5),(6,7,8),(0,3,6),(1,4,7),(2,5,8),(0,4,8),(2,4,6))
    raw=[174,285,396,417,528,639,741,852,963]
    target=(108,148,122,140,126,112,130,104,144)
    ranks=(2,9,4,7,5,3,6,1,8)
    coords=[((a-1)%3-1,(a-1)//3-1) for a in ranks]
    basis=s.Matrix([[1,u,v] for u,v in coords])
    H=s.Matrix([[int(i in ln)-int(i in lines[0]) for i in range(9)] for ln in lines[1:]])
    must(H.rank()==6 and basis.rank()==3,'exact rank/completeness')
    must(H*basis==s.zeros(7,3),'independent basis contained in magic kernel')
    must(len(H.nullspace())==3,'kernel dimension')
    g=s.Matrix([raw[a-1] for a in ranks]); one=s.ones(9,1)
    u,v=basis[:,1],basis[:,2]; vv=s.Matrix([x*x for x in v])
    must(g==528*one+111*u+s.Rational(567,2)*v+s.Rational(81,2)*vv,'raw carry decomposition')
    must(H*vv==s.Matrix([0,0,0,0,0,0,-2]),'only anti-diagonal defect')
    p=basis*(basis.T*basis).inv()*basis.T*g
    must(p==555*one+111*u+s.Rational(567,2)*v,'projection by matrix normal equations')
    must((g-p).dot(g-p)==s.Rational(6561,2),'projection squared error')
    must(sum(raw)==4995 and 555 not in raw,'all raw arrangements excluded')
    # Full 9! census, independent of any center theorem/discovery routine.
    grids=[]; visited=0
    for a in permutations(target):
        visited+=1
        if any(sum(a[i] for i in ln)!=378 for ln in lines): continue
        grids.append(a)
    must(visited==362880 and len(grids)==8,'full nine factorial census')
    must(target in grids,'requested arrangement present')
    # The eight solutions form the dihedral orbit, not merely a common count.
    def rotate(a): return tuple(a[(2-j)*3+i] for i in range(3) for j in range(3))
    def reflect(a): return tuple(a[3*i+2-j] for i in range(3) for j in range(3))
    orbit=set(); x=target
    for _ in range(4): orbit.update((x,reflect(x))); x=rotate(x)
    must(set(grids)==orbit,'all magic target arrangements form D4 orbit')
    t=s.symbols('t'); U=t*(t*t-9)*(3*t**4-57*t*t+124)/560
    V=-t*(t*t-16)*(t**4-12*t*t+67)/840
    nodes=(-3,2,4,1,0,-1,-4,-2,3)
    must(s.Matrix([U.subs(t,x) for x in nodes])==u,'U trace interpolant')
    must(s.Matrix([V.subs(t,x) for x in nodes])==v,'V trace interpolant')
    must(s.expand(4*U+18*V+t**3-15*t)==0,'identity as exact polynomials, not point samples')
    Z=s.prod(t-j for j in range(-4,5))
    must(s.expand(Z-t*(t*t-1)*(t*t-4)*(t*t-9)*(t*t-16))==0,'node annihilator')
    dimensions=[]
    for degree in range(9):
        evaluation=s.Matrix([[x**j for j in range(degree+1)] for x in nodes])
        dimensions.append(degree+1-(H*evaluation).rank())
    must(dimensions==[1,1,1,2,2,2,2,3,3],'polynomial magic degree filtration')
    must(s.expand(U).coeff(t,7)==s.Rational(3,560),'leading U coefficient')
    must(s.expand(V).coeff(t,7)==s.Rational(-1,840),'leading V coefficient')
    values=sorted({Fraction(p,q) for p in range(-8,9) for q in range(1,5)})
    tested=0
    forbidden={Fraction(-2),Fraction(-1),Fraction(-1,2),Fraction(1,2),Fraction(1),Fraction(2)}
    for p,q in product(values,repeat=2):
        out=[p*((a-1)%3-1)+q*((a-1)//3-1) for a in range(1,10)]
        must((len(set(out))==9)==(p!=0 and q!=0 and p/q not in forbidden),'rational injectivity condition')
        must(all(x<y for x,y in zip(out,out[1:]))==(p>0 and q>2*p),'rational monotonicity condition')
        must(min(out)==-abs(p)-abs(q),'positivity extreme')
        tested+=1
    # Mathematical hostiles and boundaries, not altered receipt comparison alone.
    must(H*s.Matrix(nodes)!=s.zeros(7,1),'trace grid itself not magic')
    must(H*s.Matrix([x*x for x in target])!=s.zeros(7,1),'squaring does not transport magic')
    must(H*s.Matrix([x+1 for x in target])==s.zeros(7,1),'center is freely shiftable')
    rotated=tuple(126+4*((a-4)%9%3-1)+18*(((a-4)%9)//3-1) for a in ranks)
    must([sum(rotated[i] for i in ln) for ln in lines]==[378]*7+[324],'digit rotation breaks anti-diagonal')
    return {'full_permutations':visited,'target_magic_arrangements':len(grids),'nullity':3,
            'polynomial_magic_dimensions_degree_0_to_8':dimensions,'rational_parameter_pairs':tested}

R=s.sqrt(13)
def simp(x): return s.expand(x)
def eq(x,y=0): return s.simplify(x-y)==0
def sign(x):
    x=s.simplify(x); z=s.sign(x)
    must(z in (-1,0,1),'exact algebraic sign determined')
    return int(z)
def vec(a,b): return (b[0]-a[0],b[1]-a[1])
def det(a,b): return simp(a[0]*b[1]-a[1]*b[0])
def orient(a,b,c): return sign(det(vec(a,b),vec(a,c)))
def on(p,a,b):
    if orient(a,b,p): return False
    return sign((p[0]-a[0])*(p[0]-b[0])+(p[1]-a[1])*(p[1]-b[1]))<=0

def geometry_tests():
    # Independent field arithmetic: SymPy exact radical expressions. y is scaled by D>0.
    D2=(R-1)/6; A=(R+1)/4; B=s.Rational(3,2); H=(R+7)/4
    top=[(-1,0),(-A,B),(s.Rational(-1,2),H),(s.Rational(1,2),H),(A,B),(1,0)]
    bot=[(x,-y) for x,y in top]
    segments=list(zip(top,top[1:]))+list(zip(bot,bot[1:]))+[((-1,0),(0,0)),((0,0),(1,0))]
    must(sign(D2)>0,'positive coordinate scale')
    for a,b in segments:
        dx,dy=vec(a,b); must(eq(dx*dx+D2*dy*dy,1),'unit original segment')
    contacts=0
    for (a,b),(c,d) in combinations(segments,2):
        oa,ob,oc,od=orient(a,b,c),orient(a,b,d),orient(c,d,a),orient(c,d,b)
        must(not(oa*ob<0 and oc*od<0),'no proper crossing')
        contact={p for p in (a,b,c,d) if on(p,a,b) and on(p,c,d)}
        must(len(contact)<=1,'no collinear overlap')
        if contact:
            contacts+=1; p=next(iter(contact)); must(p in (a,b) and p in (c,d),'endpoint-only contacts')
    for polygon in (top,bot):
        turns=[orient(polygon[i-1],polygon[i],polygon[(i+1)%6]) for i in range(6)]
        must(len(set(turns))==1 and turns[0]!=0,'strict six-corner convexity')
    # Area as triangulation from (-1,0), separate from half-edge face enumeration.
    twice_area=sum(det(vec(top[0],top[i]),vec(top[0],top[i+1])) for i in range(1,5))
    area=-twice_area/2
    must(eq(area,(23+5*R)/8),'scaled triangle-sum area')
    product_area=s.simplify(D2*area**2)
    must(eq(product_area,(89+26*R)/16),'n12 exact product')
    must(sign(product_area-1)>0,'each equal field area exceeds one')
    vertices=set(p for e in segments for p in e)
    must(len(vertices)==11 and len(segments)==12,'theta E V')
    must(12-11+1==2,'theta cycle rank two')
    # Standalone radical comparison and exact n12 balanced competitors.
    must(s.Rational(18,5)**2<13 and s.Rational(7,4)**2>3,'strict radical bounds')
    must(26*s.Rational(18,5)-7-48*s.Rational(7,4)==s.Rational(13,5),'strict gap certificate')
    best=6+3*s.sqrt(3)
    for alternative in (s.Rational(27,4),s.Integer(1),s.Rational(9,256)):
        must(sign(best-alternative)>0,'n12 balanced competitors lower')
    # Verify the stated derivative identities symbolically.
    x=s.symbols('x',positive=True); t=s.pi/x; f=s.log(x*s.cot(t)/4)
    first=(1+2*t/s.sin(2*t))/x
    second=-(1+4*t/s.sin(2*t)-4*t*t*s.cos(2*t)/s.sin(2*t)**2)/x**2
    must(s.trigsimp(s.diff(f,x)-first)==0,'strictly increasing derivative formula')
    must(s.trigsimp(s.diff(f,x,2)-second)==0,'strict concavity derivative formula')
    # Actual T-junction cactus: two length-one sides plus one-third base per triangle.
    T=s.sqrt(35)/6
    triangles=[((0,0),(s.Rational(1,6),T),(s.Rational(1,3),0)),
               ((s.Rational(2,3),0),(s.Rational(5,6),T),(1,0))]
    for a,b,c in triangles:
        for p,q in ((a,b),(b,c)):
            dx,dy=vec(p,q); must(eq(dx*dx+dy*dy,1),'T fixture unit sloped fences')
        dx,dy=vec(a,c); must(eq(dx*dx+dy*dy,s.Rational(1,9)),'T cycle base not unit')
        twice=det(vec(a,b),vec(a,c)); must(eq(twice**2/4,s.Rational(35,1296)),'T fixture field area square')
    must(on((s.Rational(1,3),0),(0,0),(1,0)),'first interior base junction')
    must(on((s.Rational(2,3),0),(0,0),(1,0)),'second interior base junction')
    return {'unit_segments':12,'distinct_vertices':11,'endpoint_contacts':contacts,
            'cycle_rank':2,'product':str(product_area),'field_area_squared':str(product_area),
            'T_fixture_cycle_edge_lengths':['1','1','1/3']}

if __name__=='__main__':
    result={'status':'PASS','magic':magic_tests(),'fence':geometry_tests(),'checked_conditions':CHECKS,
            'method':'No discovery code imported; exact SymPy algebra and independent integer enumeration.'}
    result['checker_sha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    out=Path(sys.argv[1]); out.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    print(json.dumps(result,sort_keys=True))
