"""Independent native pair-gcd entry audit; no producer imports."""
from math import gcd,prod,comb
from itertools import combinations
from fractions import Fraction as F
import sys
sys.stdout.reconfigure(newline="\n")
Q=91**6;CHECKS=0

def need(x,msg):
 global CHECKS
 CHECKS+=1
 if not x:raise ArithmeticError(msg)

def primes(n):return [p for p in range(2,n+1) if all(p%d for d in range(2,int(p**.5)+1))]
SUMS={1}
for p in primes(356):
 if p%3==2:SUMS={s*p**e for s in SUMS for e in range(3) if s*p**e<=356}

def components(row):
 groups=[{i} for i in range(len(row))]
 for i,j in combinations(range(len(row)),2):
  if (row[i]+row[j])//gcd(row[i],row[j]) in SUMS:
   A=next(z for z in groups if i in z);B=next(z for z in groups if j in z)
   if A is not B:A.update(B);groups.remove(B)
 return sorted(tuple(sorted(z)) for z in groups)

def ceil(a,b):return -((-a)//b)
def mixed_coeff(A,B,Y):
 # Use the coefficient of the larger pair entry as the modular variable.
 d=gcd(A,B);p,q=sorted((A//d,B//d));delta=gcd(d,Y)
 c=d//delta;x=Y//delta
 if c>Q:return None
 s=0 if p==1 else (x*pow(q,-1,p))%p
 r=(x-q*s)//p
 low=max(ceil(-Q-s,p),ceil(r-Q,q));high=min((Q-s)//p,(r+Q)//q)
 if low>high:return None
 return c,r-high*q,s+high*p

def norm(x):return min(x%1,(-x)%1)

def main():
 need(sum(gcd(a,s-a)==1 for s in SUMS for a in range(1,(s+1)//2))==5855,'complete atlas census')
 need(7 not in SUMS and 8 not in SUMS and 356 in SUMS,'split prime and exponent-three hostiles')
 rows=[]
 for a in range(2,7):
  b=13-a;G={11:2,10:4,9:9,8:30,7:90}[b];den=7*(b+1)
  cutoff=a*Q//(den*355**(a-2));Bmin=ceil(G*a,den)
  rows.append((a,b,G,Bmin,cutoff))
 need([r[-1] for r in rows]==[13520696477,62323312,257485,1007,3],'automatic no-unit max cutoffs')
 # A signed box must retain both the reduced pair and coefficient height.
 # Check its full central interval and its first absent positive integer.
 boxes=0
 for B in range(2,19):
  for q in range(2,B+1):
   for p in range(1,q):
    if gcd(p,q)!=1:continue
    C={p*r+q*s for r in range(-B,B+1) for s in range(-B,B+1)}
    R=B*(p+q)-(p-1)*(q-1)
    need(R>B*q,'strict central radius beyond Qq')
    need(all(x in C for x in range(-R,R+1)) and R+1 not in C,'complete signed-box interval')
    boxes+=1
 qs=(179,181,183,185,187);P=prod(qs)
 V=tuple(sorted([P]+[(356-q)*(P//q) for q in qs]))
 U=(2,3,4,6,18,486,9234);t=1768827685;g=1
 D,A,B=min((gcd(x,y),x,y) for x,y in combinations(V,2))
 physical=tuple(t*x for x in V)+U
 need(gcd(*V)==gcd(*U)==gcd(t,g)==1 and 1 not in V and 1 not in U,'primitive nonunit components')
 need(D==5929017 and min(V)==185370716505,'independent six-star and selected pair')
 need(min(V)>3*Q//28 and D*max(U)<=3*Q//28,'past former minimum cutoff, inside new pair gate')
 need((3*Q//28)//D==10261,'exact structured maxU cutoff')
 need(len(set(physical))==13 and sum(physical)<=Q**2,'native physical box')
 need(components(physical)==[tuple(range(6)),tuple(range(6,13))],'full correct decoder atlas gives6+7')
 for X in (V,U):
  for x,y in combinations(X,2):need(max(x,y)//gcd(x,y)<=Q,'all internal pair heights')
 mixed=0
 for X,Y in ((physical[:6],physical[6:]),(physical[6:],physical[:6])):
  for x,y in combinations(X,2):
   for z in Y:
    need(mixed_coeff(x,y,z) is None,'all bounded mixed supports absent')
    mixed+=1
 need(mixed==231,'both crossing orientations complete')
 maxima={}
 for size,cap in ((7,90),(8,30),(9,9),(10,4),(11,2),(12,1)):
  m=max(gcd(*x) for x in combinations(physical,size));maxima[size]=m
  need(m<=cap,'every actual subset scalar gcd cap')
 phase=F(1,5)+F(1,2*t)
 clearance=min(norm(x*phase) for x in physical)
 need(clearance==F(70752184,353765537)>F(1,14),'literal strict full-row phase')
 need(t*D>Q,'missing distinguished coefficient forces physical scale')
 need(F(6*t,56*max(U))>1,'larger-component safe arc exceeds grid mesh')
 # A whole six-component may not use a grow-to-seven gcd bound.
 need(t>31950 and components(physical)[0]==tuple(range(6)),'whole-component cap hostile')
 # Native normalization including shared factors of t and the nonunit label.
 scaled=0
 for tt in range(1,14):
  for gg in range(1,14):
   if gcd(tt,gg)!=1:continue
   for DD in range(1,11):
    for u in (2,3,5,7):
     delta=gcd(tt*DD,gg*u);c=tt*DD//delta;x=gg*u//delta
     need(c*gg*u==tt*DD*x,'minimal physical coefficient identity')
     need(c<=tt*DD and x<=gg*u,'normalization inequalities without unit assumption')
     scaled+=1
 print('INDEPENDENT_AUDIT: native no-unit pair-gcd criterion and concrete balanced entry PASS')
 print('ATLAS: generated from inert-prime products with exponents0,1,2;5855 ratios')
 for a,b,G,Bmin,L in rows:print('SPLIT',a,b,'failure_scale_cap',G,'pair_upper_min',Bmin,'automatic_maxU',L)
 print('SIGNED_BOX_CONTROLS:',boxes,'PHYSICAL_NORMALIZATIONS:',scaled)
 print('V:',V,'U:',U,'SCALES:',t,g,'PAIR_GCD:',D,'MAXU_CUTOFF:',10261)
 print('PHYSICAL_SUM:',sum(physical),'COMPONENTS:6+7 MIXED_SUPPORTS:',mixed,'SUBSET_GCD_MAXIMA:',maxima)
 print('PHASE:',phase,'CLEARANCE:',clearance)
 print('HOSTILES: whole six-component scale exceeds31950; prime7 and inert cube8 excluded from actual atlas')
 print('ACTIVE_CHECKS:',CHECKS)

if __name__=='__main__':main()
