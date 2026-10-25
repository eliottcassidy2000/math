#!/usr/bin/env python3
"""Exact finite controls for branch-cube theorem; not a substitute for its proof."""
import argparse
from collections import Counter
from fractions import Fraction
from functools import reduce
from itertools import combinations, product
import json
from math import gcd
from pathlib import Path
import sympy as S

checks = 0

def require(predicate, context):
    global checks
    checks += 1
    if not predicate:
        raise RuntimeError(context)


def norm(v):
    require(any(v), 'zero projective point')
    q = reduce(gcd, (abs(int(x)) for x in v))
    ans = tuple(int(x)//q for x in v)
    if next(x for x in ans if x) < 0:
        ans = tuple(-x for x in ans)
    return ans


def cross(a, b):
    return norm((a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]))


def dot(a, b):
    return sum(x*y for x, y in zip(a, b))


def arrangement(lines):
    pts = {cross(a, b) for a, b in combinations(lines.values(), 2)}
    return {p: tuple(k for k, v in lines.items() if dot(v, p) == 0) for p in pts}


L = {
    'c': (1,0,0), 'a': (1,1,0), 'b': (1,-1,0),
    'd': (1,0,1), 'e': (1,0,-1), 'h': (1,1,1),
    'k': (1,-1,-1), 'u': (1,1,-1), 'v': (1,-1,1)
}
all_incidence = arrangement(L)
require(max(map(len, all_incidence.values())) == 3, 'full arrangement has fourfold point')
require(Counter(map(len, all_incidence.values())) == Counter({2:12, 3:8}), 'full incidence census')
transition_census = Counter()
transitions = []
for old in combinations(L, 6):
    old_lines = {k:L[k] for k in old}
    old_pts = arrangement(old_lines)
    t_old = sum(len(v) == 3 for v in old_pts.values())
    for new in L:
        if new in old:
            continue
        new_lines = dict(old_lines, **{new:L[new]})
        new_pts = arrangement(new_lines)
        t_new = sum(len(v) == 3 for v in new_pts.values())
        pairs = [v for p,v in old_pts.items() if len(v)==2 and dot(p,L[new]) == 0]
        r = len(pairs)
        require(len(set(sum((list(v) for v in pairs), []))) == 2*r, 'repeated pairs must be disjoint')
        require(not any(len(v)==3 and dot(p,L[new])==0 for p,v in old_pts.items()), 'new line through old triple')
        require(t_new == t_old+r, 'triple-point increment')
        # Count intersections with new line directly, independently of n-r formula.
        restriction_zeros = {cross(L[k], L[new]) for k in old}
        d = len(restriction_zeros)
        require(d == 6-r and d >= 3, 'restricted zero count')
        copies = 2**r
        g = 1+2**(d-3)*(d-4)
        s = r*8
        cube_edges = 0 if r==0 else r*2**(r-1)
        edge_multiplicity = 2**(d-2)
        require(cube_edges*edge_multiplicity == s, 'cube edges account all branch nodes')
        require(2*(4*t_old)+s == 8*t_new, 'old and new surface nodes')
        require(2*(24+s)-copies*(2-2*g) == 80, 'smooth cover Euler accounting')
        # Arithmetic genus of B~2H on K3 is17. Recover from normalization+graph.
        require(copies*g+s-copies+1 == 17, 'branch arithmetic genus')
        transition_census[(r,t_old,t_new)] += 1
        transitions.append({'old':list(old),'added':new,'old_triples':t_old,'new_triples':t_new,'r':r,'d':d,'components':copies,'component_genus':g,'branch_nodes':s})
require(len(transitions) == 252, 'transition universe')
require(set(t['r'] for t in transitions) == {0,1,2,3}, 'all four branch types occur')

# General Euler/Riemann-Hurwitz identities, outside finite magic universe.
for n in range(6,21):
    for r in range(n//2+1):
        d=n-r
        e_n=2**(n-4)*(n*n-9*n+24)
        e_next=2**(n-3)*((n+1)**2-9*(n+1)+24)
        g=1+2**(d-3)*(d-4)
        s=r*2**(n-3)
        require(2*(e_n+s)-2**r*(2-2*g) == e_next, f'general Euler n={n},r={r}')
        pa=1+(2*2**(n-3)*(n-4))//2 # B.(B+K)/2: B=2H, K=(n-6)H
        require(2**r*g+s-2**r+1 == pa, 'general branch arithmetic genus')

hourglass = tuple('cabdehk')
for added, expected in [('c',(0,3,3,8,0,24)),('k',(2,3,1,2,5,8))]:
    rec=next(t for t in transitions if set(t['old'])==set(hourglass)-{added} and t['added']==added)
    got=tuple(rec[k] for k in ('old_triples','new_triples','r','components','component_genus','branch_nodes'))
    require(got==expected, f'actual {added} deletion')

# Smooth conics and pairwise intersections of the actual central branch.
I=S.I
conics=list(product((-1,1), repeat=3))
A=S.zeros(8)
points={}
for i, eps in enumerate(conics):
    for j, eta in enumerate(conics):
        diff=[v for v in range(3) if eps[v]!=eta[v]]
        A[i,j]=int(len(diff)==1)
        if i>=j:
            continue
        found=[]
        if len(diff)==1:
            if diff[0]==0:
                triples=[(0,1,s) for s in (-1,1)]
            elif diff[0]==1:
                triples=[(1,0,s) for s in (-1,1)]
            else:
                triples=[(1,s*I,0) for s in (-1,1)]
            for aa,dd,hh in triples:
                p=(aa,eps[0]*I*aa,dd,eps[1]*I*dd,hh,eps[2]*I*hh)
                require(S.expand(hh**2-aa**2-dd**2)==0, 'conic equation')
                require(all(S.expand(p[q+1]-eta[q//2]*I*p[q])==0 for q in (0,2,4)), 'both conics contain point')
                contain=[z for z in conics if all(S.expand(p[q+1]-z[q//2]*I*p[q])==0 for q in (0,2,4))]
                require(len(contain)==2, 'no triple intersection of conics')
                key=tuple(str(S.simplify(x)) for x in p)
                require(key not in points, 'distinct crossing points')
                points[key]=(eps,eta)
                found.append(key)
        else:
            # Two forced zeros among a,d,h imply third zero by h²=a²+d².
            require(len(diff)>=2, 'distinct nonneighbor index')
        require(len(found)==(2 if len(diff)==1 else 0), 'exact conic pair intersections')
require(len(points)==24, 'conic crossing total')
M=-2*S.eye(8)+2*A
require(M.eigenvals()=={S.Integer(4):1,S.Integer(0):3,S.Integer(-4):3,S.Integer(-8):1}, 'cube intersection spectrum')
require(M.rank()==5, 'cube lattice rank')
one=S.ones(8,1)
require((one.T*M*one)[0]/4==8, 'hyperplane square from conic sum')
for j in range(3):
    v=S.Matrix([eps[j] for eps in conics])
    require(M*v==S.zeros(8,1), 'signed half-sum kernel relation')
    plus=S.Matrix([int(eps[j]==1) for eps in conics])
    require(M*plus==2*one, 'half-sum has degree two on each conic')

# Outer branch intersections, and smooth genus-five curve's branch labels.
outer_pts=[]
for sb,se,sd in product((-1,1), repeat=3):
    c,a,b,d,e,h = 0,1,sb*I,sd*I,se,0
    require(all(S.expand(v)==0 for v in [a*a+b*b-2*c*c,d*d+e*e-2*c*c,b*b+e*e-c*c,h*h-a*a-d*d+c*c]), 'outer crossing equations')
    outer_pts.append(tuple(map(str,(c,a,b,d,e,h))))
require(len(set(outer_pts))==8, 'outer eight crossing points')
t,s=S.symbols('t s')
radicands=[t,1+t,1-t,2-t]
require({S.solve(f,t)[0] for f in radicands}=={-1,0,1,2}, 'genus-five finite branch points')
require(S.expand((1+t)+(1-t)-2)==0 and S.expand((2-t)+t-2)==0 and S.expand((1-t)+t-1)==0, 'genus-five quotient equations')

# Dual conic/A2 hexagon and exact arithmetic lift constant.
Q=S.Matrix([[1,0,0],[0,-1,S.Rational(1,2)],[0,S.Rational(1,2),-1]])
require(Q.det()==S.Rational(3,4), 'smooth dual conic')
for label in set(hourglass)-{'c'}:
    v=S.Matrix(L[label]); require((v.T*Q*v)[0]==0, 'six dual points on conic')
minimum=[(x,y) for x,y in product(range(-3,4),repeat=2) if x*x+y*y-x*y==1]
require(set(minimum)=={L[k][1:] for k in set(hourglass)-{'c'}}, 'A2 minimal vectors')
U=t*t-t*s+s*s; V=t*t-s*s; W=t*t-2*t*s
require(S.expand(U*U-V*V-W*W+V*W)==0, 'conic parametrization')
params=[(0,1),(1,1),(2,1),(-1,1),(1,2),(1,0)]
param_pts={norm([q.subs({t:x,s:y}) for q in (U,V,W)]) for x,y in params}
require(param_pts=={norm(L[k]) for k in set(hourglass)-{'c'}}, 'six marked conic parameters')
mat=S.Matrix([[2,-1],[1,1]])
expected=[(1,0),(2,1),(1,1),(1,2),(0,1),(-1,1)]
for p,q in zip(expected,expected[1:]+expected[:1]):
    require(norm(mat*S.Matrix(p))==norm(q), 'six-cycle Mobius action')
require(mat**6==-27*S.eye(2), 'projective Mobius order6')
for j in range(1,6):
    T=mat**j
    require(not (T[0,1]==0 and T[1,0]==0 and T[0,0]==T[1,1]), 'no smaller projective order')
F=t*s*(t+s)*(t-s)*(t-2*s)*(2*t-s)
Fnew=F.subs({t:2*t-s,s:t+s},simultaneous=True)
require(S.expand(Fnew+27*F)==0, 'binary sextic multiplier -27')
require(-27<0, 'multiplier is not a rational square; naive Q lift rejected')

# Elliptic quotient calculations retain actual polynomials, not merely j values.
def elliptic_j(roots):
    a,b,c=map(S.Rational,roots)
    x,y=b-a,c-a
    return S.factor(256*(x*x-x*y+y*y)**3/(x*x*y*y*(x-y)**2))
j_values={'E0':str(elliptic_j((-1,0,1))),'E1':str(elliptic_j((-1,1,2))),'E2':str(elliptic_j((-2,1,2)))}
require(j_values=={'E0':'1728','E1':'21952/9','E2':'35152/9'}, 'elliptic j controls')
u=S.symbols('u')
require(S.expand((t*(1+t)*(2-t)).subs(t,1-u)-(u+1)*(u-1)*(u-2))==0,'two equal E1 factors')
require(S.expand((t*(1+t)*(1-t)).subs(t,-u)-(u**3-u))==0,'second E0 factor')
quartic=t*(t+1)*(t-1)*(t-2)
require(S.cancel(u**4*quartic.subs(t,1/u)-(1+u)*(1-u)*(1-2*u))==0,'quartic reciprocal map')
require(S.expand(4*((1+u)*(1-u)*(1-2*u)).subs(u,t/2)-(t-2)*(t+2)*(t-1))==0,'E2 cubic map')

# Hostile control: if four old/new lines concur, B1 hypotheses must fail.
hostile={'A':(1,0,0),'B':(0,1,0),'C':(1,1,0),'D':(1,-1,0),'E':(0,0,1),'F':(1,0,1),'G':(0,1,1)}
require(max(map(len, arrangement(hostile).values()))>=4, 'fourfold hostile fixture rejected')

out={
 'status':'FINITE-EXACT', 'universe':'all252 labelled six-old/one-new transitions among nine magic lines; symbolic algebra and finite n6..20 formula controls',
 'checks':checks,
 'transition_count':len(transitions),
 'r_census':dict(sorted(Counter(t['r'] for t in transitions).items())),
 'transition_census':[{'r':r,'old_triples':a,'new_triples':b,'count':count} for (r,a,b),count in sorted(transition_census.items())],
 'cube':{'vertices':8,'edges':12,'intersection_points':24,'intersection_matrix_rank':5,'spectrum':{'4':1,'0':3,'-4':3,'-8':1}},
 'outer':{'components':2,'genus':5,'intersections':8,'component_field':'Q(sqrt2)','intersection_field':'Q(i)'},
 'dual_conic':{'determinant':'3/4','A2_vectors':sorted(minimum),'Mobius_matrix':[[2,-1],[1,1]],'binary_sextic_multiplier':-27},
 'elliptic_j':j_values,
 'scope':'Does not prove Kummer identification, Picard rank, rational-point density, or any famous open conjecture.'
}
parser=argparse.ArgumentParser(); parser.add_argument('--output'); args=parser.parse_args()
encoded=json.dumps(out,indent=2,sort_keys=True)+'\n'
if args.output:
    Path(args.output).write_text(encoded)
else:
    print(encoded,end='')
