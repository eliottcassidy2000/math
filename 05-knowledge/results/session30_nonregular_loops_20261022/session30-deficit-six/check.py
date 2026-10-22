"""Exact algebra controls; not a tournament census or closure audit."""
import json
import sympy as s
z,t=s.symbols('z t'); p=1-z-z**3; u=p.subs(z,z**2)
R=z**5+z**6+z**7
F9=(z**5+z**7)/u+z**10*(1+z+z**3)/u**2
F15=z**6/u
L=z**5*(1/p**2-1/u)/2
G=L+F9+F15
D=s.expand(p*p*u*u);N=s.cancel(G*D).expand()
def require(b):
    if not b: raise RuntimeError('exact algebra check failed')
require(s.degree(N,z)==18 and s.degree(D,z)==18)
require(s.gcd(N,D)==1)
q=1-t-t**3
require(s.cancel((F9+F9.subs(z,-z))/2-z**10/u**2)==0)
require(s.cancel((F9-F9.subs(z,-z))/2-((z**5+z**7)/u+(z**11+z**13)/u**2))==0)
cs=s.symbols('c15 c18 c20 c21 c24 c28'); orders=[15,18,20,21,24,28]
C=sum(c*z**n for c,n in zip(cs,orders))
require(s.cancel(u*((G+C/u)-G)-C)==0)
# Coefficients by polynomial recurrence, avoiding formal-series floats.
d=[int(D.coeff(z,i)) for i in range(19)]; nums=[int(N.coeff(z,i)) for i in range(19)]
a=[]
for n in range(81):
    a.append((nums[n] if n<len(nums) else 0)-sum(d[j]*a[n-j] for j in range(1,min(n,18)+1)))
require(all(x>=0 for x in a))
require(a[5:15]==[1,2,3,4,7,11,18,29,47,71])
# Exact normalized n=6 source, score test and Hamilton subset DP.
rows=[36,17,10,19,37,14]; n=6
require(sorted(x.bit_count() for x in rows)==[2]*3+[3]*3)
dp={(1<<i,i):1 for i in range(n)}
for mask in range(1,1<<n):
 for end in range(n):
  v=dp.get((mask,end),0)
  for j in range(n):
   if not(mask>>j&1) and rows[end]>>j&1:
    dp[(mask|1<<j,j)]=dp.get((mask|1<<j,j),0)+v
require(sum(dp.get(((1<<n)-1,j),0) for j in range(n))==45)
print(json.dumps({'status':'conditional global classification; exact algebra and known-family controls','numerator':str(N),'denominator':str(D),'gcd':1,'minimal_recurrence_order':18,'conditional_sharp_onset':19,'coefficients_0_to_40':a[:41],'remaining_source_orders':orders,'known_n6_H':45},indent=2))
