from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'ans_inverse'))
from fixtures import state,forward,H

def K(w):
 s=0
 while w:
  v=w//2
  if w&1:s+=v+(v&1)
  w=v
 return s

def bits(w):return [(w>>i)&1 for i in range(w.bit_length())]
def rev(w,r):return int(f'{w:0{r}b}'[::-1],2)
L=(20,1,2,7,14);R=(16,1,3,7,14)
assert state(L,3)==(9,1,23,5)
assert state(R,3)==(9,1,25,5)
examples=[]
for r in range(2,43,4):
 q=(2**(r+1)+2)//5
 assert 5*q==2**(r+1)+2 and q%2==0
 u=rev(q,r)
 j=(r-2)//4
 assert u==(7*16**j-2)//5
 assert K(u)==(14*(16**j-1)-10*j)//25
 C=23*2**r-25*u+5*K(u)+(u&1)
 assert C==(299*16**j+36-10*j)//5+(j==0)
 assert C>0 and C%2==0 and u<2**(r-1)
 for V in [0,1,2,3,C//2,10**20+1]:
  wl=2**r+2**(r+2)*V;wr=u+2**(r+2)*V
  hl=9+23*wl-5*K(wl)-(wl&1)
  hr=9+25*wr-5*K(wr)-(wr&1)
  assert K(wl)==K(V)+2*V
  assert K(wr)==K(u)+K(V)+4*q*V
  assert hl-hr==C-2*V
  if V==C//2:
   assert hl==hr
   examples.append((r,u,V,hl))
  if (r==2 and V<=34) or (r==6 and V<=1):
   assert H(forward(L,3,bits(wl)))==hl
   assert H(forward(R,3,bits(wr)))==hr
print({'status':'PASS','literal_source_profiles':[(9,1,23,5),(9,1,25,5)],'first_recontacts':examples[:4], 'largest_r':42})
