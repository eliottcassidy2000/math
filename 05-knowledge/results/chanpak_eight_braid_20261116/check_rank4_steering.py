"""Replay exact semantic upper rank and 512 positive modular rank minors."""
import sys,json
from pathlib import Path
import numpy as np
directory=Path(__file__).resolve().parent
for prev in ('target_ribbon_root','chanpak_target_ribbon_audit_20261116'):
 sys.path.insert(0,str(directory.parent/prev))
from verify_rank_alias import A,B,O,rank,piv,Matrix,v0,h
from cap_check import need
from check_sharp_eight import dp

def matword(w):
 P=Matrix.eye(len(v0))
 for a in w:P=Matrix(A[a])*P
 return P

def adjacency(w):
 v=[0]*len(w)
 for j,s in enumerate(w):
  for i in range(j):
   if j-i<=3 and s>>(j-i-1)&1:v[j]|=1<<i
   else:v[i]|=1<<j
 return v

def strong(w):
 r=adjacency(w)
 for i in range(len(w)):r[i]|=1<<i
 for j in range(len(w)):
  for i in range(len(w)):
   if r[i]>>j&1:r[i]|=r[j]
 return all(x==(1<<len(w))-1 for x in r)

def has_rank_four(Q,p):
 for j in range(4):
  ix=np.flatnonzero(Q[j:,j])
  if not len(ix):return False
  i=j+ix[0];Q[[j,i]]=Q[[i,j]]
  row=Q[j]*pow(int(Q[j,j]),-1,p)%p
  Q[j+1:]=(Q[j+1:]-Q[j+1:,j,None]*row[None,:])%p
 return True

def main():
 a=json.loads((directory/'RANK4_MOD1024.json').read_text())
 word=tuple(a['prefix']);need(word==(0,0,3,5,2,4) and a['modulus']==1024)
 P=matword(word[1:]);K=O*P*B;sel=piv(K)
 need(rank(K)==4 and len(sel)==4 and strong(word))
 p=10007
 C=np.array(P*B[:,list(sel)],dtype=np.int64)%p
 Om=np.array(O,dtype=np.int64)%p
 old=np.array(P*Matrix(v0),dtype=np.int64).reshape(-1)%1024
 need(int(h@old)==25 and len(a['rows'])==512)
 for r,ww in zip(range(1,1024,2),a['rows']):
  need(len(ww)==5 and all(1<=x<=7 for x in ww))
  v=old.copy();Y=C.copy()
  for x in ww:v=A[x]@v%1024;Y=A[x]@Y%p
  need(int(h@v)%1024==r and has_rank_four(Om@Y%p,p))
  # Independent subset-path controls mixed through the residue range.
  if r%64==1 or r==1023:
   need(strong(word+tuple(ww)))
   need(sum(dp(adjacency(word+tuple(ww)))[-1])%1024==r)
 print("All 512 odd classes mod 1024 realized in five strong append steps with exact rank four throughout.")
if __name__=='__main__':main()
