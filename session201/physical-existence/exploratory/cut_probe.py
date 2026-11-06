from itertools import product
from fractions import Fraction
from probe import T,B,state,step
sig={'O':'Z','Z':'O','J':'J'}
def scan(L,seed):
 best=(Fraction(0),None)
 for le in range(L+1):
  for letters in product('OZJ',repeat=le):
   w=''.join(letters);left=[T]
   for c in w:left.append(step(left[-1],c))
   right=[seed]
   for c in w[::-1]:right.append(step(right[-1],sig[c]))
   right=right[::-1]
   s=left[-1];n=s[0]-s[1]+s[2]if seed==T else s[0]+2*s[2]
   H=min(max(x[2],y[2])for x,y in zip(left,right));ratio=Fraction(H*H,n-1)
   if ratio>best[0]:best=(ratio,dict(word=w,n=n,H=H,heights=[(x[2],y[2])for x,y in zip(left,right)]))
 print('seed',seed,'best',best,flush=True)
if __name__=='__main__':scan(10,T);scan(10,B)
