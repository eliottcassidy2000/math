"""Reproduce the exact ranks of the first-level and two-letter residue fans."""
from pathlib import Path
import sys
HERE=Path(__file__).resolve().parent
for folder in (HERE,HERE.parent/'target_ribbon_root',
               HERE.parent/'chanpak_target_ribbon_audit_20261116'):
    if (folder/'verify_rank_alias.py').exists():
        sys.path.insert(0,str(folder))
        break
from verify_rank_alias import A,O,B,h,v0,rank,Matrix
import numpy as np
def need(c):
    if not c: raise RuntimeError('certificate failure')
U=(0,1,3,3,2,5);V=(0,0,3,5,2,4)
W=[(5,7,257),(7,6,163),(6,2,229),(7,5,295),(3,3,137),(3,7,107),
 (6,4,237),(6,3,303),(7,0,81),(6,0,115),(6,6,277),(7,4,183),
 (7,1,249),(7,7,219),(4,6,317),(5,0,95)]
S=(0,3,5,2,4)
RET=[(0,4,7,16129),(1,5,5,14371),(0,4,5,13637),(1,6,5,13799),
 (0,4,6,13801),(1,6,7,19019),(1,7,3,13229),(1,4,7,19311),
 (1,4,6,15761),(0,6,6,13939),(1,7,5,12949),(5,4,3,29527),
 (0,6,7,18009),(2,6,6,16827),(0,6,5,14845),(1,7,7,14975)]
L=np.array(O,dtype=np.int64);R0=np.array(B,dtype=np.int64)
def walk(xs):
 P=R0.copy();q=np.array(v0,dtype=np.int64)
 for a in xs[1:]:
  M=A[a]
  # Check nonnegative integer matrix multiplications are in int64 range.
  need(int(max(M.sum(axis=1)))*max(int(P.max()),int(q.max()))<2**63)
  P=M@P;q=M@q
 need(int(max(L.sum(axis=1)))*int(P.max())<2**63)
 K=Matrix(L@P)
 return int(h@q),rank(K),tuple(K.T.to_DM().rref()[0].to_Matrix())
def main():
 need(walk(U)[:2]==walk(V)[:2]==(25,4))
 U1=[walk(U+(a,))[:2] for a in range(8)]
 V1=[walk(V+(a,))[:2] for a in range(8)]
 need(U1==list(zip((25,61,43,53,103,115,95,81),(3,3,3,3,4,4,4,4))))
 need(V1==list(zip((25,43,61,53,103,95,115,81),(4,)*8)))
 image=set(); k0=walk(V)[2]; ret=0
 for i,(a,b,H) in enumerate(W):
  n,r,k=walk(V+(a,b))
  need((n,r)==(H,4) and n%32==2*i+1)
  image.add(k);ret+=(k==k0)
 need(len(image)==16 and ret==0)
 for i,(a,b,c,H) in enumerate(RET):
  n,r,k=walk(V+(a,b,c)+S)
  need((n,r,k)==(H,4,k0) and n%32==2*i+1)
 need(walk(V+(3,)+S)==(2707,4,k0))
 q=np.array(v0,dtype=np.int64)
 for a in (V+RET[0][:3]+S)[1:]:q=(A[a]@q)%32
 z=[]
 for a,b,c,_ in RET:
  qq=q.copy()
  for aa in (a,b,c)+S:qq=(A[aa]@qq)%32
  z.append(int(h@qq)%32)
 need(set(z)=={1,3,5,9,11,13,17,21,27,31})
 print('PASS: all 16 odd mod-32 residues at the original image, using common refresh suffix',S)
 print('Selected 16 reset blocks after residue-1 branch cover only',sorted(set(z)))
 print('rank-four endpoint images of the unrefreshed fan:',len(image),'exact returns:',ret)
 print('PASS: all one-letter masks and rank-four witnesses for 16 odd mod-32 residues.')
if __name__=='__main__':main()
