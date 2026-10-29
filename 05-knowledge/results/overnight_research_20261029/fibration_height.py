"""Exact degree, coprimality and specialization census for the ICARM4239 family."""
from fractions import Fraction as F
from math import gcd
from pathlib import Path
from collections import Counter
import json

def trim(a):
 a=list(a)
 while len(a)>1 and a[-1]==0:a.pop()
 return a
def plus(*args):
 a=[0]*max(map(len,args))
 for b in args:
  for i,x in enumerate(b):a[i]+=x
 return trim(a)
def scale(a,c):return trim([c*x for x in a])
def mul(a,b):
 c=[0]*(len(a)+len(b)-1)
 for i,x in enumerate(a):
  for j,y in enumerate(b):c[i+j]+=x*y
 return trim(c)
def power(a,n):
 b=[1]
 for _ in range(n):b=mul(b,a)
 return b
def evaluate(a,t):
 y=0
 for x in a[::-1]:y=y*t+x
 return y
def rem(a,b,p):
 a=trim([x%p for x in a]);b=trim([x%p for x in b])
 while len(a)>=len(b) and a!=[0]:
  q=a[-1]*pow(b[-1],-1,p)%p;k=len(a)-len(b)
  for j,x in enumerate(b):a[j+k]=(a[j+k]-q*x)%p
  a=trim(a)
 return a
def gcdmod(a,b,p):
 a=trim([x%p for x in a]);b=trim([x%p for x in b])
 while b!=[0]:a,b=b,rem(a,b,p)
 return [(x*pow(a[-1],-1,p))%p for x in a]

L=[32416,-828066,466548]
p=[103008,-1076788,-478380]
q=[81622,-180205,82612]
D=[-53374213248,239872553168,11566653808,-91335649060,-26666192028]
E=[-24372268764,76235905664,-104259187849,118525002574,-47423161625]
B=plus(mul(p,E),mul(q,D),mul(mul(p,q),plus(L,scale(p,-1),scale(q,-1))))
b2=plus(power(L,2),scale(plus(D,E),4));b4=plus(scale(mul(D,E),2),mul(L,B));b6=power(B,2)
c4=plus(power(b2,2),scale(b4,-24));c6=plus(scale(power(b2,3),-1),scale(mul(b2,b4),36),scale(b6,-216))
jnum=scale(power(c4,3),1728);jden=plus(power(c4,3),scale(power(c6,2),-1))

def main():
 if (len(c4)-1,len(c6)-1,len(jnum)-1,len(jden)-1)!=(8,12,24,24):raise RuntimeError('degrees')
 certprime=next(r for r in (101,103,107,109) if c4[-1]%r and c6[-1]%r and gcdmod(c4,c6,r)==[1])
 x=scale(mul(p,q),-1);y=mul(p,plus(mul(p,q),scale(E,-1)))
 identity=plus(power(y,2),mul(plus(mul(L,x),B),y),scale(mul(mul(x,plus(x,D)),plus(x,E)),-1))
 if identity!=[0]:raise RuntimeError('marked rational section')
 out=[]
 for height in (5,10,20,40):
  values=Counter();singular=0;parameters=0
  for den in range(1,height+1):
   for num in range(-height,height+1):
    if gcd(num,den)!=1:continue
    parameters+=1;t=F(num,den);bottom=evaluate(jden,t)
    if bottom==0:singular+=1;continue
    values[evaluate(jnum,t)/bottom]+=1
  if max(values.values())>24:raise RuntimeError('fiber degree')
  out.append({'parameter_height':height,'primitive_parameters':parameters,'singular':singular,
              'distinct_j':len(values),'max_observed_fiber':max(values.values())})
 result={'source':'https://elliptic-rank.icarm.cloud/curve/4239',
         'polynomials_ascending':{'L':L,'p':p,'q':q,'D':D,'E':E,'c4':c4,'c6':c6},
         'degrees':[8,12,24,24],'coprimality_mod_prime':certprime,'gcd_mod_prime':[1],
         'leading_coefficients_nonzero_mod_prime':True,'section_identity':True,
         'finite_parameter_census':out,'proved_height_count_exponent':'1/12',
         'scope':'Geometric degree and standard height-count consequence. No generic-rank computation or rank-jump tail law.'}
 Path(__file__).with_suffix('.json').write_text(json.dumps(result,indent=2)+'\n')
 print(json.dumps({k:v for k,v in result.items() if k!='polynomials_ascending'},indent=2))
if __name__=='__main__':main()
