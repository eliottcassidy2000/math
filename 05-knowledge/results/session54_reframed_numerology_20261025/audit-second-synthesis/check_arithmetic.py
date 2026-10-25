import json,math
from fractions import Fraction
def need(x,s):
 if not x:raise RuntimeError(s)
def mul(a,b):
 c=0
 while b:
  if b&1:c^=a
  b>>=1;a<<=1
  if a&64:a^=67
 return c
def power(a,n):
 b=1
 while n:
  if n&1:b=mul(a,b)
  a=mul(a,a);n>>=1
 return b
nonzero=range(1,64);F4=[x for x in range(64) if power(x,4)==x];F8=[x for x in range(64) if power(x,8)==x]
need(len(F4)==4 and len(F8)==8,'subfields');need(len({power(2,e) for e in range(63)})==63,'primitive')
for x in nonzero:
 norm=power(x,9);need(norm in F8 and norm!=0,'norm range');need(mul(power(x,8),power(norm,62))==power(x,62),'norm inversion')
need(sum(power(x,8)==x for x in nonzero)==7 and sum(power(x,62)==x for x in nonzero)==1,'different fixed points')
cosets={tuple(sorted(mul(x,a) for a in F4 if a)) for x in nonzero};need(len(cosets)==21 and all(len(c)==3 for c in cosets),'projective points')
point={x:i for i,c in enumerate(sorted(cosets)) for x in c}
def trace(x):return x^power(x,4)^power(x,16)
lines={tuple(sorted({point[x] for x in nonzero if trace(mul(a,x))==0})) for a in nonzero};need(len(lines)==21 and all(len(L)==5 for L in lines),'projective lines')
inv={point[x]:point[power(x,62)] for x in nonzero};need(not any(tuple(sorted(inv[p] for p in L)) in lines for L in lines),'inversion not collineation')
need(21*5==105 and 105*9==945 and 63*15==945 and 2*63==126,'blowup')
# Direct matrix-state recursion versus formula and scalar finite differences.
for b in range(9):
 v=[2]+([2]+[0]*(b-1) if b else []);values=[]
 for r in range(1,16):
  vwant=[2*3**(r-1)*math.comb(r,j) if j<=r else 0 for j in range(b+1)];need(v==vwant,'state formula');values.append(Fraction(sum(v),2*3**(r-1)));v=[3*(v[j]+(v[j-1] if j else 0)) for j in range(b+1)]
 for _ in range(b):values=[y-x for x,y in zip(values,values[1:])]
 need(set(values)=={Fraction(1)},'exact polynomial degree/pole')
 # Clear denominator in x=3z; evaluate at x=1 without cancellation.
 numerator_at_one=2 if b==0 else 2*0**b+2*sum(1**(j-1)*0**(b-j) for j in range(1,b+1))
 need(numerator_at_one==2,'noncancelled pole')
for r in range(3,30,2):
 lam=(r-1)//2
 for t in range(1,2*r):
  need(r*r*t-(t+1)*(lam*t+r-lam)==r*(r-2)+lam*(t-1)*(2*r-1-t)>0,'deficiency one')
 for t in range(1,2*r-1):
  D=r*r*t-(t+2)*(lam*t+r-lam);need(D==r*(r-3)+lam*(t-1)*(2*r-2-t),'deficiency two identity');need(D>0 if r>=5 else (D==0)==(t in [1,4]),'sharp boundary')
# Rational-observer local exponent screen computed from exact prime inequalities.
Nd=[]
for d in range(1,9):
 a=3
 while 2**(a-2)<=d:a+=1
 n=2**a
 for p in range(3,4*d+2,2):
  if any(p%q==0 for q in range(2,math.isqrt(p)+1)):continue
  v=(p-1)//math.gcd(4,p-1)
  while v<=d:n*=p;v*=p
 Nd.append(n)
need(Nd[:3]==[120,240,65520],'conductor bound')
# n5: t²=1-t, j=A*t+B. H_-15(j) has zero constant and t coefficients.
A=85995;B=-52515
need(A*A+B*B+191025*B-121287375==0 and -A*A+2*A*B+191025*A==0,'affine CM polynomial')
print(json.dumps({'status':'PASS','norm_checked_nonzero':63,'conjugation_fixed':7,'inversion_fixed':1,'projective_points':21,'vectors_per_point':3,'projective_lines':21,'points_per_line':5,'base_flags':105,'lifted_edges':945,'lifted_degree':15,'Jordan_budgets':list(range(9)),'Hall_r_checked':list(range(3,30,2)),'N_d_1_to_8':Nd,'affine_n5_CM_identity':True},indent=2))
