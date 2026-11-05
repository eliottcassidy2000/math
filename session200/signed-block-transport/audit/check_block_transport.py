"""Independent exact affine block and physical fiber transport checks."""
from math import gcd
import json

def need(q,msg):
 if not q:raise RuntimeError(msg)
def apply(s,word):
 D,E,U,d=s
 for a in word:
  if a=='O':D,E,U,d=D,0,2*U,d
  elif a=='Z':D,E,U,d=D-E+U,d,2*U-d,d
  else:raise ValueError(a)
 return D,E,U,d

def g(n):
 acc=0
 while n>1:
  if n%2==0:n//=2
  else:n=(n+1)//2;acc+=n-n%2
 return 1+acc

def select(a,b):
 j=1
 while 2*(b-2**j*(b-a))>=b:j+=1
 t=b-2**j*(b-a);ell=1
 while 2**(ell+1)*t<=b:ell+=1
 return ell,j,2**ell*t,2**j*(b-a)+a-2*b

slopes=0;odd_slopes=0;largest=0
basis=[tuple(int(i==j) for i in range(4)) for j in range(4)]
for b in range(3,402,2):
 for a in range(b//2+1,b):
  if gcd(a,b)!=1:continue
  ell,j,aa,c=select(a,b);word='O'*ell+'Z'*j
  need(ell>=1 and j>=1 and ell+j<=2*b.bit_length(),'cost')
  need(b<2*aa<2*b and aa%2==0 and gcd(aa,b)==1,'effective slope')
  need(c%2==a%2,'shift parity')
  for s in basis:
   D,E,U,d=s;DD,EE,UU,dd=apply(s,word)
   need(a*UU-b*DD==aa*U-b*D+c*d,('affine identity',a,b,s))
   need(EE==d and dd==d,'endpoint')
  slopes+=1;odd_slopes+=a%2;largest=max(largest,ell+j)
# Complete bounded two-way arithmetic fiber equivalence, not used as proof.
bounded=0;positive=[]
for n in range(3,200002,2):
 left=(5*n-9*g(n)==24)
 q=(n+1)//16
 right=((n+1)%16==0 and q>=3 and q%2==1 and 8*q-9*g(q)==29)
 need(left==right,('two-way fiber',n))
 if left:positive.append(n)
 bounded+=1
q=675385;D=g(q)
need(D==600339 and 8*q-9*D==29,'positive predecessor')
out=apply((D,1,q,1),'OOOZ')
need(out==(6003419,1,10806159,1),'positive image')
need(g(out[2])==out[0] and 5*out[2]-9*out[0]==24,'positive physical membership')
j=(q-7)//18
need(j==37521 and q==18*j+7 and D==16*j+3,'index')
need(79232*j+30133==2972894005,'positive target')
need(g(7)!=3,'excluded j0')
print(json.dumps({'status':'PASS','maximum_odd_denominator':401,'coprime_slopes':slopes,'odd_numerator_slopes':odd_slopes,'largest_selected_block_length':largest,'fiber_odd_n_checked':bounded,'bounded_positive_fiber_inputs':positive,'positive_predecessor':(D,1,q,1),'positive_image':out,'positive_index_j':j,'C2_target':2972894005},sort_keys=True,indent=2))
