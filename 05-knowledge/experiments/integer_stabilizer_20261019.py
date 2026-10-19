#!/usr/bin/env python3
"""Focused exact controls for moving-stabilizer integer coalescence."""
from fractions import Fraction as Q
from math import gcd
from itertools import product
from collections import Counter
import json
CHECKS=0
def check(test,label):
 global CHECKS
 CHECKS+=1
 if not test:raise RuntimeError(label)
def T(n):return (3*n+1)/2 if n.numerator%2 else n/2

def step(k,e,b):
 if e.numerator%2==0:return k, e/2 if b==0 else (3*e+1-Q(3)**k)/2
 if b==0:return k+1,(3*e+1)/2
 return k-1,(e-Q(3)**(k-1))/2

def traj(u,v,k,e,steps):
 rows=[]
 for t in range(steps+1):
  check(u==Q(3)**k*v+e,'literal affine representation')
  rows.append({'time':t,'u':str(u),'v':str(v),'k':k,'e':str(e),
               'literal_merge':u==v,'identity_state':k==0 and e==0})
  if t<steps:
   k,e=step(k,e,v.numerator%2);u,v=T(u),T(v)
 return rows
rows=traj(Q(45),Q(15),0,Q(30),16)
check(rows[6]['u']==rows[6]['v']=='20','time-six actual merge')
check(rows[6]['k']==-1 and rows[6]['e']=='40/3','nonidentity stabilizer target')
check(all(row['literal_merge'] and not row['identity_state'] for row in rows[6:]),'stable literal but nonidentity tail')
reparam=traj(Q(45),Q(15),1,Q(0),6)
check(reparam[6]['identity_state'],'different initial presentation reaches identity')

# Fixed-point fibre formula modulo an odd 11-power.
fibres=[]
for A in (1,2,3):
 M=11**A
 for k in (-55,-5,-1,0,1,4,5,10,55):
  a=pow(3,k,M)
  d=gcd(a-1,M)
  for b in (0,1,11,121,1331):
   fixed=[v for v in range(M) if ((a-1)*v+b)%M==0]
   expected=d if b%d==0 else 0
   check(len(fixed)==expected,'exact modular stabilizer fibre')
   for v in fixed:check(((a-1)*v+b)%M==0,'replayed fixed point')
   fibres.append({'A':A,'k':k,'b':b,'fixed_points':len(fixed)})

# Rational stabilizer classification for reduced denominators.
rational_cases=0
for denominator in range(1,51):
 if denominator%2==0:continue
 for numerator in (-7,-1,0,1,7):
  v=Q(numerator,denominator)
  d=v.denominator
  while d%3==0:d//=3
  for k in range(-12,13):
   e=(1-Q(3)**k)*v
   den=e.denominator
   while den%3==0:den//=3
   check((den==1)==(pow(3,k,d)==1%d),'rational stabilizer membership')
   rational_cases+=1

# Moving stabilizer differs from coefficient identity in every tested quotient.
k=-1;r=1;v=20;a=40
D=(3**(r+k)-3**r)*v+a
check(D==0,'cleared exact displacement')
for A in (1,2,3,4):
 M=11**A
 check(((pow(3,k,M)-1)*v+a*pow(3,-r,M))%M==0,'modular point stabilizer passes')
 check(pow(3,k,M)!=1,'modular group identity rejects')

# Height condition is necessary: nonzero M survives congruence.
for M in (11,121,1331):
 check(M%M==0 and not abs(M)<M,'strict height boundary hostile')
 for D in range(-M+1,M):
  if D%M==0:check(D==0,'height-aware exact reconstruction')

def affine(x,N):
 """Actual integer branch itinerary and its independently retained coefficients."""
 a=C=0;bits=[];start=x
 for j in range(N):
  bit=x%2;bits.append(bit)
  if bit:x=(3*x+1)//2;a+=1;C=3*C+2**j
  else:x//=2
 check(2**N*x==3**a*start+C,'actual itinerary affine coefficient')
 return a,C,x,bits

def direct(x,N):
 """Fraction path; numerator parity is valid only with odd denominator."""
 check(x.denominator%2==1,'actual 2-adic rational source')
 odds=0
 for _ in range(N):
  bit=x.numerator%2;odds+=bit;x=T(x)
 return x,odds

# Sharp maximum normalized branch offset, exhaustively checked independently
# over all bit words at depths 1..9, not only positive integer paths.
offset_controls=[]
for N in range(1,10):
 Hmax=Q(0)
 for bits in product((0,1),repeat=N):
  a=C=0
  for j,bit in enumerate(bits):
   if bit:a+=1;C=3*C+2**j
  check(C<=2**(N-a)*(3**a-2**a),'fixed odd-count maximal affine offset')
  H=Q(C,3**a);Hmax=max(Hmax,H)
  check(H<=Q(2**N,6),'sharp normalized affine offset bound')
 check(Hmax==Q(2**N,6),'normalized bound attained by last odd bit')
 offset_controls.append({'N':N,'maximum_normalized_offset':str(Hmax)})

# Complete isolated-positive-integer atlas at each depth, with all actual
# branch cylinders retained. A separate bounded-head scan verifies completeness.
atlas=[]
for N in range(1,13):
 M=2**N
 for d in (1,2,3,15,30,61):
  counts=Counter();integer_roots=[];rational_denominators=Counter()
  H=Q(d,2)+Q(M,12)
  for residue in range(M):
   a,Cu,u,_=affine(residue+d,N);b,Cv,v,_=affine(residue,N)
   A=3**a-3**b;B=3**a*d+Cu-Cv
   if A==0:
    kind='identity_cylinder' if B==0 else 'parallel_no_collision'
    counts[kind]+=1
    # Positive lifts cannot confuse an isolated coincidence with identity.
    for lift in (1,2):
     x=residue+lift*M
     U,_=direct(Q(x+d),N);V,_=direct(Q(x),N)
     check((U==V)==(B==0),'equal-slope positive-lift control')
    continue
   counts['unequal_slope_cylinder']+=1
   x=Q(-B,A)
   check(abs(A)%x.denominator==0,'isolated denominator divides slope difference')
   if x>0:
    check(x<H,'isolated positive source height bound')
    if a>b:check(x<Q(M,12),'upper odd-excess sharper bound')
   lawful=x.denominator%2==1 and (x.numerator-residue*x.denominator)%M==0
   if not lawful:continue
   counts['lawful_rational_root']+=1
   U,aa=direct(x+d,N);V,bb=direct(x,N)
   check(U==V and (aa,bb)==(a,b),'branch-valid rational root replay')
   if x<=0:continue
   counts['lawful_positive_rational_root']+=1
   rational_denominators[x.denominator]+=1
   if x.denominator!=1:continue
   integer_roots.append(int(x))
   U,_=direct(x+d+M,N);V,_=direct(x+M,N)
   check(U!=V,'one whole-cylinder step destroys isolated equality')
  independently_found=[]
  for x in range(1,(H.numerator+H.denominator-1)//H.denominator):
   U,a=direct(Q(x+d),N);V,b=direct(Q(x),N)
   if U==V and a!=b:independently_found.append(x)
  check(sorted(integer_roots)==independently_found,'complete bounded-head atlas control')
  check(sum(counts[key] for key in ('identity_cylinder','parallel_no_collision','unequal_slope_cylinder'))==M,
        'all branch cylinders classified exactly once')
  atlas.append({'N':N,'translation':d,'height_upper_bound':str(H),
                'cylinder_counts':dict(counts),'isolated_positive_integer_roots':sorted(integer_roots),
                'lawful_positive_rational_denominators':dict(rational_denominators)})

# The actual 223/233 clock repair is an isolated unequal-clock root.
A=15;B=7;d=10;source=223
a,Cu,u,bu=affine(source+d,A);b,Cv,v,bv=affine(source,B)
L=2**B*3**a-2**A*3**b
R=2**A*Cv-2**B*(3**a*d+Cu)
check(u==v==425,'actual unequal-clock 223/233 join')
check(L!=0 and Q(R,L)==source,'unique branch-valid unequal-clock root')
check(L//2**B==-3**5*(2**8-3**5),'unequal-clock resonance determinant')
clock_example={'u':233,'v':223,'u_clock':A,'v_clock':B,'u_odd_count':a,'v_odd_count':b,
               'u_affine_constant':Cu,'v_affine_constant':Cv,'common_endpoint':u,
               'root_equation_coefficient':L//2**B,'root_equation_rhs':R//2**B,
               'u_parity_word':bu,'v_parity_word':bv}

# Actual integer sharpness, not merely formal real branch words.  N=6m gives
# d<M, a last-odd upper word, two first-odd lower bits, and first merge at N.
sharpness=[]
for N in range(6,61,6):
 M=2**N;x=5*(M-1)//9;d=(17*M+10)//18
 check(9*x==5*(M-1) and 18*d==17*M+10,'sharpness source integrality')
 check(0<d<M,'sharpness within one source period')
 a,Cu,u,bu=affine(x+d,N);b,Cv,v,bv=affine(x,N)
 check(bu==[0]*(N-1)+[1] and bv==[1,1]+[0]*(N-2),'lawful extremal itinerary pair')
 check(u==v==5 and (a,b)==(1,2),'sharpness integer merge')
 check(Q(x)-Q(d,2)==Q(M,12)-Q(5,6),'sharp asymptotic coefficient one twelfth')
 U=Q(x+d);V=Q(x)
 for _ in range(N):
  check(U!=V,'sharpness is a first collision');U=T(U);V=T(V)
 sharpness.append({'N':N,'lower_source':x,'translation':d,'upper_source':x+d,
                   'first_merge_time':N,'common_endpoint':u})
print(json.dumps({'status':'FINITE-EXACT controls, not integer convergence',
                  'checks':CHECKS,'literal_merge_not_identity':rows,
                  'reparametrized_identity':reparam,'rational_membership_cases':rational_cases,
                  'modular_fixed_point_fibres':fibres,'normalized_offset_controls':offset_controls,
                  'branch_valid_isolated_collision_atlas':atlas,
                  'unequal_clock_example':clock_example,'integer_height_sharpness_family':sharpness},indent=2,sort_keys=True))
