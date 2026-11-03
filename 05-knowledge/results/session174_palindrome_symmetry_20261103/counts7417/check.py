#!/usr/bin/env python3
"""Elementary certificates and finite palindrome residue computations; stdlib only."""
import json, math
from pathlib import Path
P=7417

def check(ok,msg):
 if not ok: raise RuntimeError(msg)
def prime(n): return n>=2 and all(n%d for d in range(2,math.isqrt(n)+1))
def phi(n):return sum(math.gcd(a,n)==1 for a in range(n))
def orbit(a,m):
 seen={}; vals=[];x=1
 while x not in seen:
  seen[x]=len(vals);vals.append(x);x=x*a%m
 return {'preperiod':seen[x], 'period':len(vals)-seen[x], 'values_one_orbit':vals}
check(prime(P),'7417 prime')
cert={str(q):pow(10,2472//q,P) for q in [2,3,103]}
check(pow(10,2472,P)==1 and all(v!=1 for v in cert.values()),'order certificate')
check(2472**2>=9*P,'BHS sufficient inequality')
inv10=pow(10,-1,P)
w=[(pow(10,r,P)+pow(inv10,r,P))%P for r in range(2472)]
check([r for r,x in enumerate(w) if x==0]==[618,1854],'zero positions')
check(len(set(w))==1237,'trace image size')
check(len(set(w[:1237]))==1237,'fundamental domain injectivity')
check(all(w[(r+1236)%2472]==-w[r]%P for r in range(2472)),'half-period sign')
check(all(w[(r+1)%2472]==((10+inv10)*w[r]-w[(r-1)%2472])%P for r in range(2472)),'trace recurrence')
def census(m,q):
 check(m>=1,"census requires at least three digits")
 dp=[0]*q;dp[0]=1
 for j in range(m):
  weight=(pow(10,j,q)+pow(10,2*m-j,q))%q
  nxt=[0]*q
  for r,n in enumerate(dp):
   if n:
    for a in range(1 if j==0 else 0,10):nxt[(r+a*weight)%q]+=n
  dp=nxt
 center=pow(10,m,q);nxt=[0]*q
 for r,n in enumerate(dp):
  if n:
   for a in range(10):nxt[(r+a*center)%q]+=n
 check(sum(nxt)==9*10**m,'census total')
 return nxt
rows=[]
for m in range(1,7):
 c=census(m,P)
 rows.append({'digits':2*m+1,'total':sum(c),'occupied':sum(x>0 for x in c),'minimum':min(c),'maximum':max(c),'divisible_by_7417':c[0]})
m=618;A=10**(2*m)+1+10**m;B=2*10**(2*m)+2+10**m
check(str(A)==str(A)[::-1] and str(B)==str(B)[::-1],'palindrome collision')
check(A!=B and (B-A)%P==0,'residue collision')
mods={str(n):{'factorization':f,'unit_group_order':phi(n),'decimal_power_orbit':orbit(10,n)} for n,f in [(240,'2^4*3*5'),(1001,'7*11*13'),(2222,'2*11*101')]}
res={'prime':P,'trial_division_bound':math.isqrt(P),'p_minus_1_factorization':'2^3*3^2*103','decimal_order':2472,'order_witnesses':cert,'bhs_squared_inequality':[2472**2,9*P],'groups':{'Fp_units':7416,'decimal_subgroup':2472,'index':3,'trace_image':1237,'dihedral_order':4944,'affine_group_order':P*(P-1)},'trace_zero_distances':[618,1854],'trace_first_return':2472,'palindrome_censuses':rows,'clock_moduli':mods,'explicit_collision':{'digits':1237,'A':'10^1236+10^618+1','B':'2*10^1236+10^618+2','difference':'10^1236+1','common_residue':A%P}}
Path(__file__).with_name('certificate.json').write_text(json.dumps(res,indent=2)+'\n')
print(json.dumps(res,indent=2))
