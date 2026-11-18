from pathlib import Path
import sys
H=Path(__file__).resolve().parent
sys.path[:0]=[str(H),str(H.parent/'ans_inverse'),str(H.parent/'finite_edit_recontact')]
from fixtures import state
from unit_gain_replay import K,L,R,rev
def restrict(T,I):
 return tuple(sum(((T[i]>>j)&1)<<k for k,j in enumerate(I)) for i in I)
I=[i for i in range(5) if i!=4]
assert state(restrict(R,I),3)==(1,0,8,1)
for j in range(1,31):
 m=16**j;b=4*m;u=(7*m-2)//5;V=(299*m+36-10*j)//10
 assert 28*m<=V<32*m
 W=u+4*b*V;s=W.bit_length()-1;P=W-2**s
 assert P>=3*2**s//4
 A=1+8*W-K(W)
 d=9+25*P-5*K(P)-(P&1)
 assert d-A>=7+2**s
 print(j,'separation verified',d>A)
print('PASS')
