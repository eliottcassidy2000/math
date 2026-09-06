"""Complete carried certificates for endpoints39,45,51,57,63.

Inherits the all-h normalization/degree interface, then proves positivity
for exactly h=6,...,10 by full rational polynomial certificates. No numerical
roots or repository producer imports. The general all-h claim remains open.
"""
import argparse,json,hashlib,sys
from pathlib import Path
from math import factorial,comb,gcd
from fractions import Fraction as F
import sympy as S
sys.stdout.reconfigure(newline='\n')
x=S.Symbol('x')
ZERO=S.Poly(0,x,domain=S.QQ);ONE=S.Poly(1,x,domain=S.QQ)
COUNT=0
def gate(ok,label):
 global COUNT
 COUNT+=1
 if not ok:raise RuntimeError(label)
def fall(a,n):
 r=ONE
 for j in range(n):r*=S.Poly(a-j,x,domain=S.QQ)
 return r
def rem(a,p):
 a=list(a);h=len(p)-1
 while len(a)>=len(p):
  z=a.pop()
  for j,pj in enumerate(p[:-1]):a[len(a)-len(p)+1+j]-=z*pj
 return a+[ZERO]*(h-len(a))
def rows(h):
 p=[fall(x+h,h-j).mul_ground(S.Rational(factorial(2*h+1),factorial(3*h-3*j)*factorial(1+2*j))) for j in range(h+1)]
 q=[fall(2*x+2*h,2*h-e).mul_ground(S.Rational(1,factorial(6*h-3*e)*factorial(2+2*e))) for e in range(-1,2*h+1)]
 carry=fall(x,1).mul_ground(S.Rational(2**(h+1)*factorial(3*h),factorial(6*h+3)*factorial(2*h+1)))
 for j in range(h):carry*=S.Poly(2*x+2*j+1,x,domain=S.QQ)
 gate(p[-1]==ONE,'monic first row')
 gate(q[0]==carry*p[0],'exact inverse-carry cancellation')
 r=rem(q[1:],p)
 r=[r[j]-carry*p[j+1] for j in range(h)]
 for j,a in enumerate(r):gate(a.degree()<=2*h-j,'response weight bound')
 difference=list(q)
 for j,a in enumerate(r):difference[j+1]-=a
 gate(all(a.is_zero for a in rem(difference,p)),'full t*q=t*R in original quotient')
 return p,q,r
def mm(A,B):
 h=len(A)
 return [[sum((A[i][j]*B[j][k] for j in range(h)),ZERO) for k in range(h)] for i in range(h)]
def char(p,r):
 h=len(r);cols=[rem([ZERO]*j+r,p) for j in range(h)]
 M=[[cols[j][i] for j in range(h)] for i in range(h)]
 for i in range(h):
  for j in range(h):gate(M[i][j].degree()<=2*h+j-i,'multiplication matrix weight')
 # Faddeev-LeVerrier computes the characteristic polynomial and yields
 # the complete Cayley-Hamilton matrix identity as its terminal state.
 B=[[ONE if i==j else ZERO for j in range(h)] for i in range(h)]
 cs=[ONE]
 for k in range(1,h+1):
  B=mm(M,B)
  ck=sum((B[i][i] for i in range(h)),ZERO).mul_ground(S.Rational(-1,k))
  cs.append(ck)
  for i in range(h):B[i][i]+=ck
  gate(ck.degree()==2*h*k,'exact characteristic degree')
 gate(all(a.is_zero for row in B for a in row),'full symbolic Cayley-Hamilton')
 return cs
def literal(h,p,q):
 primitive=0;fibre_count=0
 for xx in (1,2,3,7,16):
  g=xx+3*h+1;endpoint=6*h+3
  expected_first={j:F(factorial(g),factorial(xx+j)*factorial(3*h-3*j)*factorial(1+2*j)) for j in range(h+1)}
  expected_second={e:F(factorial(2*g),factorial(2*xx+e)*factorial(6*h-3*e)*factorial(2+2*e)) for e in range(-1,2*h+1)}
  gate(list(expected_first.values())==[comb(g,2*h+1)*F(c.eval(xx)) for c in p],'literal first normalization')
  scale=factorial(2*g)//factorial(2*g-4*h-2)
  gate(list(expected_second.values())==[scale*F(c.eval(xx)) for c in q],'literal doubled normalization')
  for mass,expected,gamma_offset in ((g,expected_first,1),(2*g,expected_second,2)):
   found={}
   for na in range(mass+1):
    for nb in range(mass-na+1):
     nc=mass-na-nb
     if -endpoint*na+(2*g-endpoint)*nb+(3*g-endpoint)*nc:continue
     gate((nc-gamma_offset)%2==0,'complete charge-fibre parity')
     j=(nc-gamma_offset)//2
     found[j]=F(factorial(mass),factorial(na)*factorial(nb)*factorial(nc))
   gate(found==expected,'exhaustive nonnegative charge fibre')
   fibre_count+=len(found)
  gate(expected_second[-1]>0,'strictly present lower carry')
  primitive+=gcd(g,endpoint)==1
 # Uniform gcd hostile: g=3h+3, gcd(g,6h+3)=3; first support mass h+1.
 gh=3*h+3;a=6*h+3
 gate(gcd(gh,a)==3,'gcd-hostile content')
 counts=(1,h-1,1)
 gate(sum(counts)==gh//3 and -a+(2*gh-a)*(h-1)+(3*gh-a)==0,'gcd-hostile earlier return')
 return primitive,fibre_count
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--certificate');args=ap.parse_args()
 certificate={'schema':1,'scope':'exactly h=6,...,10; x=g-3h-1>=1; char coefficients in y=x-1',
              'degree_bound':'deg c_k<=2hk','families':[]}
 total=0
 for h in range(6,11):
  p,q,r=rows(h);cs=char(p,r);cert=[];n=0
  for k in range(1,h+1):
   shifted=S.Poly(cs[k].as_expr().subs(x,x+1),x,domain=S.QQ)
   content,primitive=shifted.primitive();co=list(reversed(primitive.all_coeffs()))
   gate(content>0,'positive polynomial content')
   for c in co:gate(c.is_Integer and c>0,'positive integer shifted coefficient')
   n+=len(co)
   cert.append({'k':k,'degree':2*h*k,'content':str(content),'coefficients_ascending':[int(c) for c in co]})
  prim,fib=literal(h,p,q)
  total+=n
  family={'h':h,'endpoint':6*h+3,'g_min':3*h+2,'characteristic':cert}
  certificate['families'].append(family)
  print(f'FAMILY h={h} endpoint={6*h+3} g>={3*h+2}: {n} strictly positive coefficients; {prim}/5 primitive literal controls; {fib} complete fibres',flush=True)
 data=(json.dumps(certificate,sort_keys=True,indent=2)+'\n').encode()
 gate(total==3170,'entire coefficient certificate cardinality')
 if args.certificate:Path(args.certificate).write_bytes(data)
 print('CERTIFICATE coefficients=3170 bytes='+str(len(data))+' sha256='+hashlib.sha256(data).hexdigest())
 print('PROOF: full polynomial identities and shifted positivity hold for every real x>=1; every real first root has q_x(root)<0.')
 print('RELATIVE TO THM-4436: integral x>=1 has h simple negative first roots; gcd(g,6h+3)=1 gives first nonzero mass g or 2g.')
 print('SCOPE: five new fixed-endpoint, unbounded-g families; all-h positivity and general trinomial separation remain OPEN.')
 print('PASS',COUNT,'always-active exact gates; no numerical roots and no repository producer imports.')
if __name__=='__main__':main()
