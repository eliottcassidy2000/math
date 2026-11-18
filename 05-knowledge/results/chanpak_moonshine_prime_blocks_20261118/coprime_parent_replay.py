from unit_gain_replay import K,L,R,state,H,forward,bits
from math import gcd

def algebra(p,w):
 h,e,t,d=p
 for b in bits(w):
  if b:h,e,t,d=h+t-e,d,2*t-d,d
  else:h,e,t,d=h,0,2*t,d
 return h,e,t,d

new=[]
for A,w in ((L,548),(R,545)):
 G=forward(A,3,bits(w));root=len(G)-1
 p=state(G,root);new.append(p)
 assert H(tuple(a&((1<<root)-1) for a in G[:-1]))==p[3]
print('literal reroots',new)
assert new==[(12233,817,12233,817),(12233,753,12233,753)]
assert gcd(new[0][3],new[1][3])==1
for (A,w),p,expected in zip(((L,548),(R,545)),new,(23649,23713)):
 G=forward(A,3,bits(w));root=len(G)-1
 assert H(forward(G,root,[1]))==expected
for j in range(10):
 r=4*j+2;B=2**r;u=(7*16**j-2)//5
 C=(299*16**j+36-10*j)//5+(j==0);V=C//2;s=V.bit_length()
 wl=B+4*B*V;wr=u+4*B*V
 # Strip last chronological bit, i.e. highest set bit.
 k=wl.bit_length();assert k==wr.bit_length()
 ps=[]
 for p,w in (((9,1,23,5),wl),((9,1,25,5),wr)):
  h,e,t,d=algebra(p,w-(1<<(k-1)))
  # Pad skipped chronological zeros through length k-1.
  used=(w-(1<<(k-1))).bit_length()
  for _ in range(k-1-used):h,e,t,d=h,0,2*t,d
  ps.append((h+t-e,h-e,h+t,h))
 a,b=ps
 assert a[0]==b[0] and a[2]==b[2]
 delta=a[3]-b[3]
 assert delta==2**s and a[1]-b[1]==delta and gcd(a[3],b[3])==1
 assert a[3]<a[0] and b[3]<b[0]
 for w in (0,1,2,3,7,29,1001):
  assert algebra(a,w)[0]-algebra(b,w)[0]==-delta*(K(w)+(w&1))
print('coprime physical parents and continuation identities: PASS')
