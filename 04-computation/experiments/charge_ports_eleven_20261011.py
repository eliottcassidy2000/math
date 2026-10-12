"""Exact checks for charge receipts and preserved-port eleventh-root contexts.
Run with python3 (and independently python3 -O) from the repository root.
All declarations and checks survive optimization. PySAT and SymPy required.
"""
from fractions import Fraction as Q
from itertools import product,combinations
from math import lcm
from collections import defaultdict
from functools import lru_cache
import hashlib
import sympy as s,json
from pysat.solvers import Solver
checks=0
def check(ok,label):
 global checks
 checks+=1
 if not ok:raise RuntimeError(label)
z=s.symbols('z'); poly=s.Poly(s.cyclotomic_poly(55,z),z)
rel=[(i,int(poly.nth(i))) for i in range(40) if poly.nth(i)]
def red(a):
 a=list(a)+[0]*max(0,40-len(a))
 for d in range(len(a)-1,39,-1):
  if a[d]:
   for i,v in rel:a[d-40+i]-=a[d]*v
 return tuple(a[:40])
O=(0,)*40;I=(1,)+(0,)*39
powers=[red([0]*k+[1]) for k in range(55)]
def add(a,b):return tuple(x+y for x,y in zip(a,b))
def sc(a,n):return tuple(n*x for x in a)
def mul(a,b):
 c=[0]*79
 for i,x in enumerate(a):
  if x:
   for j,y in enumerate(b):
    if y:c[i+j]+=x*y
 return red(c)
def conj(a):
 v=O
 for i,c in enumerate(a):
  if c:v=add(v,sc(powers[-i%55],c))
 return v
a=powers[11];g=powers[5]
r5=add(I,sc(add(a,conj(a)),2));b=O
for j in range(1,11):b=add(b,sc(powers[5*j],int(s.legendre_symbol(j,11))))
w=mul(r5,b)
eta=sc(add(sc(I,-3),w),Q(1,8));rho=sc(add(r5,b),Q(1,4));v=sc(add(sc(I,-7),sc(b,4)),Q(1,15))
P=[O]
for j in range(4):P.append(add(P[-1],powers[11*j]))
R=[O,sc(eta,-1),sc(I,Q(3,4)),sc(conj(eta),-1)]
ports=list(dict.fromkeys(P+R))
def unit(a,j=1):
 # inverse of a-g^-j, via finite geometric sum, includes all 11 roots.
 inv=O
 for k in range(11):inv=add(inv,sc(powers[(-5*j*k)%55],Q(a**(10-k),a**11-1)))
 return mul(add(sc(I,a),sc(powers[5*j%55],-1)),inv)
def graph(points):
 den=lcm(*(x.denominator for p in points for x in p if isinstance(x,Q)))
 pts=[tuple(int(x*den) for x in p) for p in points]
 filters=[]
 for prime in [331,661]:
  root=pow(int(s.primitive_root(prime)),(prime-1)//55,prime)
  check(pow(root,55,prime)==1 and pow(root,5,prime)!=1 and pow(root,11,prime)!=1,
        'filter root has exact order 55')
  check(int(poly.eval(root))%prime==0,'filter is a power-basis ring homomorphism')
  rp=[pow(root,k,prime) for k in range(40)];rc=[pow(root,-k,prime) for k in range(40)]
  filters.append((prime,[(sum(x*y for x,y in zip(p,rp))%prime,sum(x*y for x,y in zip(p,rc))%prime) for p in pts]))
 edges=[];surv=0
 for i,j in combinations(range(len(pts)),2):
  if all((vv[i][0]-vv[j][0])*(vv[i][1]-vv[j][1])%q==den*den%q for q,vv in filters):
   surv+=1;d=add(pts[i],sc(pts[j],-1))
   if mul(d,conj(d))==sc(I,den*den):edges.append((i,j))
 return edges,surv,den

def v2(n):return (n&-n).bit_length()-1

@lru_cache(None)
def orbit(n):
 states=[n];word=[]
 while n!=1:
  check(len(word)<1000,'declared orbit control completed; timeout is not rejection')
  a=v2(3*n+1);word.append(a);n=(3*n+1)>>a;states.append(n)
 return tuple(states),tuple(word)

def odd_charge(word):return sum(a%2 for a in word)

@lru_cache(None)
def carry_floor(blocks):
 out=Q(1)
 for j in range(1,blocks+1):
  r=(j+1)//2;n=6*r-1 if j%2 else 6*r+1
  out*=Q(3*(n-1),3*n+1)
 return out

def deadline(n,B):
 if n==1:return 0
 t=0
 while (n-1)*2**B*3**t>=2*4**t:t+=1
 return B+t

def floor_log_ratio(a,b):
 k=a.bit_length()-b.bit_length()
 if k>=0:
  if a<b*2**k:k-=1
 else:
  if a*2**(-k)<b:k-=1
 check(k>=0 and b*2**k<=a<b*2**(k+1),'exact rational logarithm')
 return k

def precision(n,B):
 T=deadline(n,B)
 if not T:return 1
 return 1+floor_log_ratio((3*n+1)*16**(T-1),5**(T-1))

def charge_controls():
 rows=[]
 for n in range(1,10002,2):
  states,word=orbit(n);o=odd_charge(word)
  check(len(set(states))==len(states),'strict ROOT states do not repeat')
  carry=Q(1)
  for x in states[1:-1]:carry*=Q(3*x,3*x+1)
  check(carry>carry_floor(o+1),'distinct even-block endpoint floor')
  check(carry_floor(o+1)>Q(2,o+3),'simple rational floor')
  if n in [1,3,7,27,223,4591]:
   B=o;T=deadline(n,B)
   old=1+(3*(1+2**B*(n-1))+1).bit_length()*T
   rows.append({'n':n,'charge':o,'deadline':T,'old_precision':old,
                'new_precision':precision(n,B),
                'carry_product':str(carry),'distinct_floor':str(carry_floor(o+1))})
 for n in range(1,1002,2):
  states,word=orbit(n)
  for B in range(13):
   T=deadline(n,B);D=precision(n,B);q=0;A=0;stop=0
   if n==1:continue
   for t,a in enumerate(word,1):
    q+=a%2;A+=a
    if states[t]==1 or q>B:stop=t;break
   check(0<stop<=T and A<D,'all stop events lie inside the new precision')
 # A cheap hostile to extrapolating the declared bank's log budget.
 failures=[]
 for n in range(1,328,2):
  if odd_charge(orbit(n)[1])>((n+1).bit_length()-1)+24:failures.append(n)
 check(failures==[327],'least hostile to this particular all-source log budget')
 return {'universe':'odd n <= 10001; event budgets 0..12 at odd n <= 1001',
         'examples':rows,'least_log_budget_failure':327,'failure_charge':35}

def words(L,A,prefix=()):
 if not L:
  yield prefix;return
 for a in range(1,A-L+2):yield from words(L-1,A-a,prefix+(a,))

def carrier(word):
 A=0;P=1;B=0
 for a in word:B=3*B+2**A;P*=3;A+=a
 return P,2**A,B

def replay(n,word):
 for a in word:
  check(n!=1,'no ROOT padding in a controller')
  check(v2(3*n+1)==a,'exact valuation guard')
  n=(3*n+1)>>a
 return n

def controller_controls():
 bank={n:odd_charge(orbit(n)[1]) for n in range(1,30,2)}
 tables=defaultdict(list);flat=[];levels={}
 for L in range(1,13):
  P=3**L;Amax=(P-1).bit_length()-1
  for word in words(L,Amax):
   p,Qv,B=carrier(word);check(p==P and Qv<P,'complete contracting word universe')
   k=0
   while Qv*4**k<P:
    C=Qv*4**k;D=Qv*(4**k-1)//3-B
    residue=(-D*pow(C,-1,P))%P
    row=(C,D,odd_charge(word),word,k)
    tables[L,residue].append(row);flat.append((L,residue,row))
    # The least two positive odd sources with a strictly smaller positive child.
    n=residue if residue%2 else residue+P
    while n<=1 or C*n+D<=0 or (P-C)*n<=D:n+=2*P
    for x in [n,n+2*P]:
     h=(C*x+D)//P
     check(0<h<x and h%2==1,'guarded smaller source')
     endpoint=replay(h,word)
     check(endpoint==4**k*x+(4**k-1)//3,'inverse word reaches the exact sibling')
     if k:
      a=v2(3*x+1)
      check(replay(x,(a,))==replay(endpoint,(a+2*k,)),'source retained through sibling join')
    k+=1
  levels[L]=len(flat)
 check([levels[L] for L in [2,4,6,8,10,12]]==[4,24,137,1045,6577,77441],
       'frozen complete inverse-template counts')

 def options(n,depth,extra):
  out=[]
  if n%4==1:
   a=v2(3*n+1);out.append(((3*n+1)>>a,a%2,'forward'))
  if n%8==5:out.append(((n-1)//4,0,'sibling'))
  K=v2(n+1)
  if K>=2:
   t=(n+1)>>K;r=v2(3**K*t-1)
   if r>=2:out.append(((n-1)//2,1,'reset'))
   else:
    y=(3**K*t-1)//2;a=v2(3*y+1);y=(3*y+1)>>a
    if a==6:out.append(((n-1)//2,-1,'sporadic26'))
    if a==2 and v2(3*y+1)==10:out.append(((n-1)//2,-1,'sporadic2210'))
  for L in range(1,depth+1):
   P=3**L
   for C,D,o,word,k in tables.get((L,n%P),[]):
    h=(C*n+D)//P
    if 0<h<n:out.append((h,-o,'inverse'))
  if extra>=1 and n%1024==43:out.append(((81*n+101)//512,1,'43-cylinder'))
  if extra>=2 and n%1024==39:out.append(((243*n+251)//512,3,'39-cylinder'))
  return out

 results=[];accepted_sets=[]
 for depth,extra,rate in [(2,0,1),(4,0,1),(6,0,1),(8,0,1),(10,0,1),(12,0,1),
                          (12,1,1),(12,2,3)]:
  constant=max(o-rate*((n+1).bit_length()-1) for n,o in bank.items())
  ok=dict(bank);uses=defaultdict(int);pending=[];witness=[];receipts=[]
  for n in range(31,100002,2):
   candidates=options(n,depth,extra)
   for h,d,name in candidates:
    check(0<h<n and h%2==1,'every dispatcher arrow decreases the actual source')
    check(d<=rate*(((n+1).bit_length()-1)-((h+1).bit_length()-1)),
          'charge is paid by the declared source rank')
   hits=[r for r in candidates if r[0] in ok]
   if hits:
    h,d,name=min(hits,key=lambda row:(ok[row[0]]+row[1],row[0],row[1],row[2]))
    ok[n]=ok[h]+d;uses[name]+=1
    receipts.append([n,h,d,name,ok[n]])
    # Validation uses a separate actual orbit; recognition above never queries it.
    check(ok[n]==odd_charge(orbit(n)[1]),'independent replay of the compiled charge')
    check(0<=ok[n]<=rate*((n+1).bit_length()-1)+constant,'source-only budget')
    if n in [31,39,43,57,59,327,53803]:witness.append([n,h,d,name,ok[n]])
   elif len(pending)<12:pending.append(n)
  results.append({'depth':depth,'extra_cylinders':extra,'rate':rate,'constant':constant,
                  'accepted':len(ok),'universe':50001,'uses':dict(sorted(uses.items())),
                  'first_pending':pending,'examples':witness,
                  'receipt_hash':hashlib.sha256(json.dumps(receipts).encode()).hexdigest()})
  if extra==2:final_receipts=receipts
  accepted_sets.append(set(ok))
 for x,y in zip(accepted_sets,accepted_sets[1:]):check(x<=y,'declared policy extensions are monotone')
 check([r['accepted'] for r in results[:6]]==[3648,3652,3923,4074,4074,4074],
       'frozen finite coverage controls')
 # All inverse depths and all sibling exponents at one fixed source.
 pool=sorted({x for h in range(1,39,2) for x in orbit(h)[0]})
 siblings=[];x=39
 while x<=max(pool):siblings.append(x);x=4*x+1
 check(len(pool)==55 and max(pool)==3077 and siblings==[39,157,629,2517],
       'complete smaller-source orbit union at 39')
 check(not set(pool)&set(siblings) and x==10069,'all-depth inverse-sibling shield at 39')
 check(next(j for j,y in enumerate(orbit(39)[0]) if y in pool)==5,
       'first source39 common future with any smaller source is at step5')
 for t in range(201):
  n=39+1024*t;h=19+486*t
  check(replay(n,(1,1,2,1,4))==h and 2*(h+1)<=n+1,'paid all-height 39 controller')
  n=43+1024*t;h=7+162*t
  check(replay(n,(1,2,2,4))==h and 4*(h+1)<=n+1,'paid all-height 43 controller')
 return {'bank_charges':bank,'bank_words':{n:orbit(n)[1] for n in bank},
         'final_receipts':final_receipts,'inverse_template_count':len(flat),
         'template_hash':hashlib.sha256(json.dumps(flat).encode()).hexdigest(),
         'finite_policies':results,'shield39':{'lower_orbit_union':pool,'siblings':siblings,
                                            'first_sibling_above_union':x},
         'universal_entry':'OPEN'}

def proper_colour(count,edges,k):
 clauses=[[k*i+c+1 for c in range(k)] for i in range(count)]
 clauses += [[-(k*i+c+1),-(k*j+c+1)] for i,j in edges for c in range(k)]
 with Solver(name='g3',bootstrap_with=clauses) as solver:
  check(solver.solve(),'declared colouring exists')
  model=set(solver.get_model())
  colours=[next(c for c in range(k) if k*i+c+1 in model) for i in range(count)]
 for i,j in edges:check(colours[i]!=colours[j],'every complete geometric edge is coloured')
 return colours

def port_controls():
 # All labelled six-colour boundary assignments on the retained pentagon.
 lifted=0
 for colouring in product(range(6),repeat=5):
  if any(colouring[i]==colouring[(i+1)%5] for i in range(5)):continue
  # Greedy distinct representatives: before the last choice at most four
  # previously used colours plus this port's colour are forbidden.
  satellite=[]
  for c in colouring:
   satellite.append(next(j for j in range(6) if j!=c and j not in satellite))
  check(len(set(satellite))==5 and all(a!=b for a,b in zip(satellite,colouring)),
        'universal one-layer extension colours')
  lifted+=1
 check(lifted==3120,'all labelled proper C5 six-colourings')
 # Exact degree-five barrier and all unit norms, independent of finite graphs.
 check(s.totient(55)==40 and 40==5*8,'relative extension has odd prime degree')
 check(mul(w,w)==sc(I,-55),'marked rhombus carrier')
 for vv in [eta,rho,v]+[unit(2,j) for j in range(1,11)]:
  check(mul(vv,conj(vv))==I,'generator is a genuine unit direction')
 # The mixed generator is not a product chosen merely by field name.
 num=add(I,sc(add(a,g),-1))
 expr=sum(c*z**j for j,c in enumerate(conj(num)) if c)
 ip=s.Poly(s.invert(expr,poly.as_expr(),z),z)
 u=mul(num,tuple(Q(ip.nth(i)) for i in range(40)))
 check(mul(u,conj(u))==I,'mixed eleventh/fifth root generator norm')
 up=sum(s.Rational(c)*z**j for j,c in enumerate(u) if c)
 ub=sum(s.Rational(c)*z**j for j,c in enumerate(conj(u)) if c)
 check(s.rem(up*ub,poly.as_expr(),z)==1,'independent symbolic norm calculation')
 den=lcm(*(Q(x).denominator for x in u))
 check(den==60002251,'mixed generator exact denominator')
 check(add(add(I,mul(u,conj(a))),mul(u,conj(g)))==add(add(u,a),g),
       'six-unit additive relation')
 # sigma fixes alpha and beta, and sends gamma to gamma^3.
 e=next(e for e in range(55) if e%5==1 and e%11==3)
 sigma=O
 for j,c in enumerate(u):
  if c:sigma=add(sigma,sc(powers[e*j%55],c))
 check(sigma!=u,'mixed generator genuinely escapes the degree-eight field')
 ds=[powers[11*j] for j in range(5)]+[powers[5*j] for j in range(11)]
 ds += [u,mul(u,conj(a)),mul(u,conj(g))]
 ds=list(dict.fromkeys(ds+[sc(vv,-1) for vv in ds]))
 D=lcm(*(Q(x).denominator for p in ds+ports for x in p))
 dd=[tuple(int(x*D) for x in p) for p in ds]
 ball={O};shell={O}
 for _ in range(3):
  shell={add(x,y) for x,y in product(shell,dd)}-ball;ball|=shell
 points=list(dict.fromkeys(ports+[tuple(Q(x,D) for x in p) for p in sorted(ball)]))
 edges,survivors,common=graph(points)
 colours=proper_colour(len(points),edges,3)
 # A second, independently checked colour lift applies to every word length
 # for these directions, not just this finite ball. It is not a colouring
 # of every unit direction in the full field.
 linear=[1,2,1,0,2,1,1,1,0,0,1,1,0,0,0,1,0,0,0,0,
         1,0,1,0,0,1,0,0,0,0,1,0,0,1,0,1,0,0,0,0]
 def residue(point):
  return sum(Q(c).numerator*pow(Q(c).denominator,-1,3)*v
             for c,v in zip(point,linear))%3
 enlarged=ds+[eta,conj(eta),sc(eta,-1),sc(conj(eta),-1)]
 for d in enlarged:check(residue(d)!=0,'all-length generator colouring certificate')
 check(residue(powers[3])==0 and mul(powers[3],conj(powers[3]))==I,
       'one root direction defeats this functional without forcing more colours')
 linear_colours=[residue(p) for p in points]
 direction_set={tuple(int(c*common) for c in d) for d in enlarged}
 for i,j in edges:
  check(linear_colours[i]!=linear_colours[j],'independent additive lift on the full graph')
  difference=tuple(int((x-y)*common) for x,y in zip(points[i],points[j]))
  check(difference in direction_set,'every observed unit direction is in the screened alphabet')
 check((len(points),len(edges))==(8436,24238),'complete three-layer graph census')
 check(all(p in points for p in P+R),'original pentagon and rhombus ports retained')
 small=points[:8]
 small_edges,_,_=graph(small)
 for i,j in combinations(range(len(small)),2):
  diff=add(small[i],sc(small[j],-1))
  check(((i,j) in small_edges)==(mul(diff,conj(diff))==I),
        'independent small complete distance check without modular rejection')
 coordinates=[tuple(int(x*common) for x in p) for p in points]
 return {'one_layer_pentagon_boundary_assignments':lifted,
         'mixed_generator_denominator':den,'mixed_generator':[str(x) for x in u],
         'all_length_linear_colour_vector':linear,
         'graph':{'vertices':len(points),'edges':len(edges),'filter_survivors':survivors,
                  'coordinate_denominator':common,
                  'coordinate_hash':hashlib.sha256(json.dumps(coordinates).encode()).hexdigest(),
                  'edge_hash':hashlib.sha256(json.dumps(edges).encode()).hexdigest(),
                  'colouring':colours,'ports':list(range(len(ports))),'chi':3},
         'six_colour_forcing':'OPEN; one-layer quintic satellites cannot force it'}

if __name__=='__main__':
 result={'charge':charge_controls(),'controllers':controller_controls(),'ports':port_controls()}
 result['checks']=checks
 print(json.dumps(result,indent=2,sort_keys=True))
