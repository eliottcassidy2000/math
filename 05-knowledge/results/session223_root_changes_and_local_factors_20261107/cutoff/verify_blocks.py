from itertools import product
from verify_drift import move,run,F

def pars(w):
 b=[int(x=='Z') for x in w];n=len(w)
 return (b[0],sum(x*2**j for j,x in enumerate(b)),sum(b[i]*b[j]*2**(i-1-j) for i in range(n) for j in range(i))+sum(b[i]*b[i-1] for i in range(1,n)),sum(x*2**(n-1-j) for j,x in enumerate(b)),2**n,b[-1])
def coeff(w):
 if not w:return 1,0,2
 i,a,c,q,L,e=pars(w);h=a+L;r=c+q+e
 return 1-i+h-r,i-1,1+h
words=['']+[''.join(x) for n in range(1,5) for x in product('OZ',repeat=n)]
checks=0
for A in words:
 s=run((1,0,2,1),A+'J');D,E,U,d=s
 for B in words[1:]:
  i,a,c,q,L,e=pars(B)
  if run(s,B)!=(D-i*E+a*U-c*d,e*d,L*U-q*d,d):raise RuntimeError(('block',A,B))
  for C in words:
   x,y,z=coeff(C)
   got=x*D-x*i*E+(x*a+z*L)*U-(x*c+z*q-y*e)*d
   if got!=F(run(s,B+'J'+C)):raise RuntimeError(('two J',A,B,C))
   checks+=1
print({'two_J_equations':checks,'status':'PASS'})
