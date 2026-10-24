import json,math
from pathlib import Path
N=120

def require(v,m):
 if not v:raise RuntimeError(m)
def add(*vs):return [sum(x) for x in zip(*vs)]
def scale(v,k):return [k*x for x in v]
def mul(a,b):
 r=[0]*(N+1)
 for i,x in enumerate(a):
  for j,y in enumerate(b[:N+1-i]):r[i+j]+=x*y
 return r
def power(a,k):
 r=[1]+[0]*N
 for _ in range(k):r=mul(r,a)
 return r
def inv(a):
 r=[1]+[0]*N
 require(a[0]==1,'unit')
 for i in range(1,N+1):r[i]=-sum(a[j]*r[i-j] for j in range(1,i+1))
 return r
def shift(a,k):return [0]*k+a[:N+1-k]
one=[1]+[0]*N
# Direct lattice, rather than divisor sums or the producer's Newton traces.
raw=[0]*(3*N+3)
B=math.isqrt(6*N)+3
for x in range(-B,B+1):
 for y in range(-B,B+1):
  z=x*x+x*y+2*y*y
  if z<len(raw):raw[z]+=1
th=raw[:N+1];th3=[th[j//3] if j%3==0 else 0 for j in range(N+1)]
A=raw[0::3][:N+1];Bv=raw[1::3][:N+1];C=raw[2::3][:N+1]
# Cubic norm identity in q^(1/3), not Newton's identities.
e1=scale(A,3);e2=scale(add(mul(A,A),scale(shift(mul(Bv,C),1),-1)),3)
e3=add(power(A,3),shift(power(Bv,3),1),shift(power(C,3),2),scale(shift(mul(mul(A,Bv),C),1),-3))
ith=inv(th)
actual=[scale(one,27),mul(add(scale(e1,9),scale(th3,-27)),ith),mul(add(scale(e2,3),scale(mul(th3,e1),-9)),power(ith,2)),mul(add(e3,scale(mul(th3,e2),-3)),power(ith,3)),scale(mul(mul(th3,e3),power(ith,4)),-1)]
# eta product t=E/theta^3; derive s from logarithmic differential.
E=one[:]
for n in range(1,N+1):
 for step in [n,7*n]:
  for _ in range(3):
   for j in range(N,step-1,-1):E[j]-=E[j-step]
E=shift(E,1);t=mul(E,power(ith,3));tunit=t[1:]+[0]
s=mul(mul([j*t[j+1] + t[j+1] if j<N else 0 for j in range(N+1)],inv(tunit)),power(ith,2))
# Last coefficient requires t at N+1; omit it for s comparisons.
expected=[scale(one,27),scale(one,0),scale(one,-18),scale(s,-8),add(scale(t,24),scale(one,-1))]
require(all(a[:N]==b[:N] for a,b in zip(actual,expected)),'T3 norm identity')
R=mul(th,inv(th3));R2=mul(R,R)
require(all(x%27==0 for x in add(R2,scale(add(one,scale(t,4),scale(power(t,2),21)),-1))),'mod27 lift')
# Match independently generated recurrence with theta(t(q)).
c=[1,2]
for n in range(1,40000):
 z,r=divmod((26*n*n+13*n+2)*c[-1]+3*(3*n-1)*(3*n-2)*c[-2],(n+1)**2)
 require(not r,'nonintegral recurrence');c.append(z)
composed=[0]*(N+1);tp=one[:]
for k in range(N+1):
 composed=add(composed,scale(tp,c[k]));tp=mul(tp,t)
require(composed==th,'theta vs recurrence')
rows=[];sums=[0]*4
for p in range(3,140,2):
 if any(p%d==0 for d in range(2,math.isqrt(p)+1)):continue
 z,r=divmod(c[p*p]-2,p*p);require(not r,'calibration');a=z*pow(2,-1,p)%p
 counts=[0]*4;first=None
 for n in range(2,40000//(p*p)+1):
  require((c[p*p*n]-c[n]-p*p*a*n*n*c[n])%(p*p*p)==0,'quadratic holdout')
  for j,v in enumerate([True,n%p==0,c[n]%p==0,n%p==0 and c[n]%p==0]):counts[j]+=int(v)
  linear=(c[p*p*n]-c[n]-p*p*a*n*c[n])%p**3
  if linear and first is None:first=[n,linear]
 require(first is not None,'linear negativecontrol');rows.append([p,a,counts,first]);sums=[x+y for x,y in zip(sums,counts)]
fixture=json.loads((Path(__file__).resolve().parent/'sources/next-digit-fixture.json').read_text())
require(len(rows)==33 and sums==[7995,1983,7162,1864],'counts')
for row,fr in zip(rows,fixture['rows']):
 require(row[0]==fr['p'] and row[1]==fr['a_p_calibrated_only_at_n1'],'calibrations')
 require(row[2]==list(fr['holdout_counts'].values()),'strata')
 require(row[3]==list(fr['first_linear_model_failure'].values()),'negativecontrol')
f3=sum(math.comb(3,k)**2*math.comb(3+k,k)*math.comb(2*k,3) for k in range(2,4))
require(f3==sum(c[k]*c[3-k] for k in range(4))==760,'binomial sharpness')
require((f3-4)%81==27,'sharpness mod81')
f5=sum(math.comb(5,k)**2*math.comb(5+k,k)*math.comb(2*k,5) for k in range(3,6))
require(f5==sum(c[k]*c[5-k] for k in range(6))==273504,'binomial naive-lift control')
require((f5-48*4)%27==18,'naive-lift rejection')
print(json.dumps({'status':'PASS','q_bound':N-1,'recurrence_bound':40000,'prime_count':len(rows),'holdout_counts':sums,'mod81_t3_residual':756%81,'rows':rows},indent=2))
