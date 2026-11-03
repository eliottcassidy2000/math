#!/usr/bin/env python3
"""Exact PG(2,3) puzzle-to-C5 tournament witnesses; Python stdlib only.
All checks remain active under -O. --mutate options deliberately MUST fail.
"""
import argparse, itertools, json, math
parser=argparse.ArgumentParser()
parser.add_argument('--mutate',choices=['bad-move','bad-conjugator','bad-adjacency'])
args=parser.parse_args()
checks=0

def require(value,message):
    global checks
    checks+=1
    if not value:raise SystemExit('FAIL: '+message)

points=tuple([(1,a,b) for a in range(3) for b in range(3)]+[(0,1,a) for a in range(3)]+[(0,0,1)])
N=len(points);identity=tuple(range(N))

def norm(v):
    for a in v:
        if a%3:
            inverse=pow(a%3,-1,3)
            return tuple(inverse*x%3 for x in v)
    raise ValueError('zero vector has no projective normalization')

def line_cross(a,b):
    x,y=points[a],points[b]
    n=((x[1]*y[2]-x[2]*y[1])%3,(x[2]*y[0]-x[0]*y[2])%3,(x[0]*y[1]-x[1]*y[0])%3)
    return tuple(i for i,p in enumerate(points) if sum(u*v for u,v in zip(n,p))%3==0)

def line_span(a,b):
    return tuple(sorted({points.index(norm(tuple(s*x+t*y for x,y in zip(points[a],points[b])))) for s,t in itertools.product(range(3),repeat=2) if (s,t)!=(0,0)}))

require(N==13 and len(set(points))==13,'13 distinct projective points')
normalized={norm(p) for p in itertools.product(range(3),repeat=3) if p!=(0,0,0)}
require(normalized==set(points),'normalization covers all nonzero F3 triples')
lines=set()
for a,b in itertools.combinations(range(N),2):
    L=line_cross(a,b)
    require(L==line_span(a,b),f'independent line constructions disagree at {a},{b}')
    require(len(L)==4 and a in L and b in L,f'line {a},{b} must have four points')
    lines.add(L)
require(len(lines)==13,'13 projective lines')
for i in range(N):require(sum(i in L for L in lines)==4,f'point {i} lies on four lines')
for a,b in itertools.combinations(range(N),2):
    require(sum(a in L and b in L for L in lines)==1,f'unique line through {a},{b}')

moves={}
for a,b in itertools.permutations(range(N),2):
    L=line_cross(a,b);c,d=[i for i in L if i not in (a,b)]
    p=list(identity);p[a],p[b]=b,a;p[c],p[d]=d,c
    if args.mutate=='bad-move' and (a,b)==(0,1):p[c],p[d]=c,d
    moves[a,b]=tuple(p)
    expected={a:b,b:a,c:d,d:c}
    require(all(p[i]==expected.get(i,i) for i in range(N)),f'move {a}->{b} violates double-transposition rule')
    require(sorted(p)==list(identity) and all(p[p[i]]==i for i in range(N)),f'move {a}->{b} is an involutory permutation')

def compose(p,q):return tuple(p[x] for x in q)

def power(p,k):
    result=identity
    for _ in range(k):result=compose(p,result)
    return result

def cycles(p):
    unseen=set(range(N));result=[]
    while unseen:
        i=min(unseen);c=[]
        while i in unseen:c.append(i);unseen.remove(i);i=p[i]
        result.append(tuple(c))
    return tuple(result)

def order(p):return math.lcm(*(len(c) for c in cycles(p)))

def word_permutation(path):
    require(path[0]==path[-1]==0,'word must be closed at hole 0')
    result=identity;hole=0;steps=[]
    for a,b in zip(path,path[1:]):
        require(hole==a and a!=b,'successive legal nontrivial hole moves')
        result=compose(moves[a,b],result);hole=b
        require(result[0]==hole,'permutation transports the hole along its declared path')
        steps.append({'from':a,'to':b,'line':line_cross(a,b),'other_swap':[i for i in line_cross(a,b) if i not in (a,b)]})
    require(result[0]==0,'closed word fixes original hole')
    return result,steps

gword=(0,1,3,2,4,0)
hword=(0,1,7,12,9,4,0)
g,gsteps=word_permutation(gword)
h,hsteps=word_permutation(hword)
require(order(g)==5,'g has order exactly 5')
require(order(h)==4,'h has order exactly 4')
if args.mutate=='bad-conjugator':h=power(h,3)
require(compose(h,g)==compose(power(g,2),h),'h*g*h^-1 must equal g^2')
require(cycles(g)==((0,),(1,4,2,12,3),(5,),(6,11,10,7,8),(9,)),'pinned g cycle decomposition')
require(cycles(h)==((0,),(1,4,12,2),(3,),(5,),(6,10,11,8),(7,),(9,)),'pinned h cycle decomposition')
V=(3,1,4,2,12)
require(all(g[V[i]]==V[(i+1)%5] for i in range(5)),'g is translation on the chosen five-point orbit')
require(all(h[V[i]]==V[(2*i)%5] for i in range(5)),'h is multiplication by 2 on the chosen orbit')
subgroup={compose(power(g,a),power(h,b)) for a in range(5) for b in range(4)}
require(len(subgroup)==20,'normal forms give 20 distinct puzzle permutations')
require(all(compose(x,y) in subgroup for x in subgroup for y in subgroup),'20-element set is closed')
require(len({tuple(p[v] for v in V) for p in subgroup})==20,'action on chosen orbit is faithful')
# Row-permutation convention: P_ij=1 exactly when g(V_i)=V_j.
P=[[int(g[V[i]]==V[j]) for j in range(5)] for i in range(5)]
P2=[[sum(P[i][k]*P[k][j] for k in range(5)) for j in range(5)] for i in range(5)]
A=[[P[i][j]+P2[i][j] for j in range(5)] for i in range(5)]
if args.mutate=='bad-adjacency':A=[[int((j-i)%5 in (1,4)) for j in range(5)] for i in range(5)]
require(all(A[i][i]==0 for i in range(5)),'loopless adjacency')
require(all(A[i][j]+A[j][i]==1 for i in range(5) for j in range(i)),'one directed edge per unordered pair')
require(all(sum(row)==2 for row in A),'regular outdegree two')
require(A==[[int((j-i)%5 in (1,2)) for j in range(5)] for i in range(5)],'A is precisely the C5 regular tournament')
S={1,2};family=[sorted({pow(2,k,5)*s%5 for s in S}) for k in range(4)]
require(len({tuple(t) for t in family})==4,'h cyclically permutes four tournament orientations')
require(set(family[2])=={-s%5 for s in S},'h^2 gives converse orientation')
require(set(family[1])!=S,'h is not a tournament automorphism')
# Independent characteristic polynomial, using the 120-term determinant expansion.
def mul(a,b):
    c=[0]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b):c[i+j]+=x*y
    return c
charpoly=[0]*6
for pi in itertools.permutations(range(5)):
    term=[(-1)**sum(pi[i]>pi[j] for i in range(5) for j in range(i+1,5))]
    for i,j in enumerate(pi):term=mul(term,[-A[i][j],1] if i==j else [-A[i][j]])
    for i,x in enumerate(term):charpoly[i]+=x
require(charpoly==[-2,-5,-5,0,0,1],'characteristic polynomial x^5-5x^2-5x-2')
require(charpoly==mul([-2,1],[1,3,4,2,1]),'cyclotomic quartic factorization')
# Concrete incidence obstruction, independent of PGL's abstract order.
L=(0,1,2,12);image=tuple(sorted(g[i] for i in L))
require(L in lines and image not in lines,'g is not a projective collineation')
pgl_order=((27-1)*(27-3)*(27-9))//(3-1)
require(pgl_order==5616 and pgl_order%5!=0,'PGL3(3) order also excludes order-five collineations')
print(json.dumps({'status':'PASS','checks_executed':checks,'points':points,'g_word':gword,'g_cycles':cycles(g),'g_steps':gsteps,'h_word':hword,'h_cycles':cycles(h),'h_steps':hsteps,'orbit_labels':V,'orbit_points':[points[i] for i in V],'orbit_actions':{'g':'i -> i+1 mod5','h':'i -> 2i mod5'},'subgroup_order':len(subgroup),'adjacency':A,'connection_set_family':family,'characteristic_polynomial_ascending':charpoly,'noncollineation_witness':{'line':L,'g_image':image},'boundary':'Puzzle action and C5 tournament construction only; no Baby Monster intertwiner or canonical correspondence.'},indent=2))
