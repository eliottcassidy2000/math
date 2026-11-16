"""Exact rank and returned image for all odd mod-2048 classes on one seed."""
import json
from check_rank4_steering import *

def main():
 z=json.loads((directory/'RETURN_MOD2048.json').read_text())
 seed=(0,0,3,5,2,4); R=(3,)+seed[1:];P=matword(seed[1:]);K=O*P*B
 need(z['refresh']==list(R) and z['modulus']==2048 and rank(K)==4)
 sel=piv(K);need(len(sel)==4)
 C0=P*B[:,list(sel)];M=matword(R)
 # All full image outputs of R are contained in E (after applying the
 # observable quotient). Also verify the inclusion directly.
 E=O*C0;ER=O*M*B;need(rank(E.row_join(ER))==4)
 p=10007; q=2048
 t=np.array(M,dtype=np.int64)
 C=np.array(C0,dtype=np.int64)%p
 orow=np.array(O,dtype=np.int64)@t%p
 ft=h@t%q
 state=np.array(P*Matrix(v0),dtype=np.int64).reshape(-1)%q
 need(int(h@state)==25 and len(z['rows'])==q//2)
 for r,w in zip(range(1,q,2),z['rows']):
  need(len(w)<=5 and all(1<=a<=7 for a in w))
  v=state.copy();Y=C.copy()
  for a in w:v=A[a]@v%q;Y=A[a]@Y%p
  need(int(ft@v)%q==r and has_rank_four(orow@Y%p,p))
 print('1024 odd mod-2048 targets; each certified return acts invertibly on the SAME rank-four image.')
if __name__=='__main__': main()
