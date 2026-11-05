#!/usr/bin/env python3
import math,json,random
from fractions import Fraction

def need(b):
 if not b:raise RuntimeError('independent algebra check failed')
def basis(k,odd,x):return x**(1 if odd else 2)*math.prod(x*x-j*j for j in range(1,k+1))
def gamma(k,odd):return math.factorial(2*k+1)*(1 if odd else k+1)
def response(k,odd,J):
 values={x:basis(k,odd,x)for x in range(-J-2,J+4)};out=[]
 for j in range(J+1):
  out.append(values[1]);values={x:x*x*values[x]+x*(values[x+1]-values[x-1])for x in range(min(values)+1,max(values))}
 return out

def newton_values(f,n):
 values=[f(x)for x in range(n+1)];out=[]
 while values:
  out.append(values[0]);values=[b-a for a,b in zip(values,values[1:])]
 return out

def run():
 rng=random.Random(1971522);columns=compilers=formal=0
 for odd in (0,1):
  for d in range(1,11):
   p=[math.factorial(k)*gamma(k,odd)for k in range(d)];raw=[response(k,odd,2*d+3)for k in range(d)]
   L=[]
   for j in range(d):
    row=[]
    for k in range(d):
     need(raw[k][j]%p[k]==0 and(j>=k or raw[k][j]==0)and(j!=k or raw[k][j]==p[k]));row.append(raw[k][j]//p[k])
    L.append(row)
   for k in range(d):need(all(v%p[k]==0 for v in raw[k]));columns+=1
   n=2*d-(1 if odd else 0)
   for trial in range(35):
    u=[rng.randrange(-25,26)for k in range(d)];u[0]=2*rng.randrange(-20,21)+1
    a=[p[k]*u[k]for k in range(d)];h=[sum(L[j][k]*a[k]for k in range(d))for j in range(d)]
    need(math.gcd(*h)==math.gcd(*a));v1=u[0];v2=sum(u[k]*basis(k,odd,2)for k in range(d));hW=v1+v2
    strong=(h[1:]if d>1 else[hW])+[3*v2-v1]
    need(math.gcd(*strong)==math.gcd(*a));formal+=1
    ar=[0]+[math.factorial(r-1)*rng.randrange(-5,6)for r in range(1,n+1)]
    c=[]
    for k in range(d):
     f=newton_values(lambda x:basis(k,odd,x),n);z=sum(ar[r]*f[r]for r in range(1,n+1));need(z%p[k]==0);c.append(z//p[k])
    coeff=[0]*d
    for k in range(d-1,-1,-1):coeff[k]=c[k]-sum(coeff[j]*L[j][k]for j in range(k+1,d))
    need(sum(coeff[j]*h[j]for j in range(d))==sum(c[k]*a[k]for k in range(d)));compilers+=1
 P=lambda x:x**10+64*x**8+1000*x**6+4020*x**4+2772*x*x
 def normalized(j,k):
  choose=lambda a,b:math.comb(a,b)if 0<=b<=a else 0
  return choose(j+k+1,2*k+2)+choose(j+k,2*k+2)
 def contents(evaluate,d):
  A=[]
  for k in range(d):
   a=evaluate(k+1)-sum(normalized(k+1,j)*A[j]for j in range(k));need(a%gamma(k,0)==0);A.append(a)
  return math.gcd(*A),math.gcd(*(math.factorial(k)*a for k,a in enumerate(A))),math.gcd(*(math.factorial(k)**2*a for k,a in enumerate(A)))
 def v3(n):
  need(n>0);v=0
  while n%3==0:n//=3;v+=1
  return v
 powers=[]
 for m in range(1,25):
  Q=lambda x:x*x*P(x)**m+x*(P(x+1)**m-P(x-1)**m)
  raw=contents(lambda x:P(x)**m,5*m);wrapped=contents(Q,5*m+1)
  need(raw==(3**(2*m),3**(2*m+1),3**(2*m+2)))
  need(tuple(map(v3,wrapped))==(2*m,2*m+1,2*m+1))
  z=4*Q(3)+5*Q(2)-11*Q(1)
  need(z==12*P(4)**m+46*P(3)**m-3*P(2)**m-21*P(1)**m and v3(z)==2*m+1)
  powers.append([m,*map(v3,raw),*map(v3,wrapped),v3(z)])
 # The strong bank does not give the same integral response lattice.
 den=9*33-15*23;alpha=Fraction(33-3*23,den);beta=Fraction(9*3-15,den)
 need((alpha,beta)==(Fraction(3,4),Fraction(-1,4)))
 print(json.dumps({'status':'PASS','evaluation_response_columns':columns,'formal_gcd_controls':formal,'integer_response_compilers':compilers,'ordered_powers':powers,'strong_bank_recovery_hostile':[str(alpha),str(beta)]},sort_keys=True,indent=2))
if __name__=='__main__':run()
