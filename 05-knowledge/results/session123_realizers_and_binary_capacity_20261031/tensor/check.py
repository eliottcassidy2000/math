from math import comb
from fractions import Fraction as Q
import json

def need(x,m):
 if not x:raise RuntimeError(m)
def mul(a,b):return [[sum(x*y for x,y in zip(row,col)) for col in zip(*b)] for row in a]
def tr(a):return list(map(list,zip(*a)))
def eye(n):return [[int(i==j) for j in range(n)] for i in range(n)]
def det(a):
 a=[list(map(Q,r)) for r in a];v=Q(1)
 for i in range(len(a)):
  p=next((j for j in range(i,len(a)) if a[j][i]),None)
  if p is None:return 0
  if p!=i:a[p],a[i]=a[i],a[p];v=-v
  z=a[i][i];v*=z;a[i]=[t/z for t in a[i]]
  for j in range(i+1,len(a)):
   z=a[j][i];a[j]=[x-z*y for x,y in zip(a[j],a[i])]
 return v

def main():
 superchecks=strict=gram=parity=0
 ns=[]
 for s in range(128):
  p=[0]
  for i in range(256):p.append(p[-1]+int(i&s==s))
  ns.append(p)
  for x in range(128):
   for y in range(128):need(p[x+y]>=p[x]+p[y],'subset superadditivity');superchecks+=1
 for n in range(1,25):
  z=[[int(i&s==s) for s in range(n)] for i in range(n)]
  zi=[[(-1)**(i.bit_count()-s.bit_count()) if i&s==s else 0 for s in range(n)] for i in range(n)]
  t=[[int(j<=i) for j in range(n)] for i in range(n)]
  ti=[[int(i==j)-int(i==j+1) for j in range(n)] for i in range(n)]
  need(mul(z,zi)==eye(n),'Boolean inverse');need(mul(t,ti)==eye(n),'difference inverse')
  k=[[2**((i&j).bit_count()) for j in range(n)] for i in range(n)]
  need(mul(z,tr(z))==k,'density Gram')
  f=mul(mul(t,k),tr(t));b=mul(t,z)
  need(f==mul(b,tr(b)) and det(f)==1,'capacity unimodular Gram')
  bi=mul(zi,ti);need(mul(f,mul(tr(bi),bi))==eye(n),'integer capacity inverse');gram+=1
  for i in range(n):
   for j in range(n):need(z[i][j]==comb(i,j)%2 if j<=i else z[i][j]==0,'Pascal parity');parity+=1
 for A in range(1,65):
  L=1<<(A-1).bit_length();vals=[sum(ns[s][A]*ns[s][B] for s in range(A)) for B in range(2*L+1)]
  for B in range(2*L+1):
   if B+L<len(vals):need(vals[B+L]==vals[B]+vals[L],'origin translation')
  for B in range(1,2*L+1):
   irreducible=all(vals[B]>vals[x]+vals[B-x] for x in range(1,B))
   need(irreducible==(B<=L),'complete generator classification');strict+=max(0,B-1)
 tensor=[[1]]
 for d in range(7):
  need(tensor==[[2**((i&j).bit_count()) for j in range(1<<d)] for i in range(1<<d)],'tensor identity')
  m=[[1,1],[1,2]];tensor=[[tensor[i//2][j//2]*m[i%2][j%2] for j in range(2*len(tensor))] for i in range(2*len(tensor))]
 return dict(status='PASS',subset_superadditivity_controls=superchecks,generator_split_controls=strict,gram_sizes_checked=gram,binomial_parity_controls=parity,tensor_depths_checked=7,scope='finite exact controls for analytic all-parameter tensor/Gram/generator proofs')
if __name__=='__main__':print(json.dumps(main(),sort_keys=True,indent=2))
