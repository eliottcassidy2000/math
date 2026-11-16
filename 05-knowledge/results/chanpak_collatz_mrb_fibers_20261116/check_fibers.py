import json,sys
from pathlib import Path
from sympy import Matrix,diag
from sympy.matrices.normalforms import hermite_normal_form

def need(c):
 if not c: raise AssertionError
def load():
 here=Path(__file__).resolve().parent
 for p in [here/'return_automaton.json',here.parent/'eight_braid_root'/'return_automaton.json',here.parent/'chanpak_eight_braid_20261116'/'return_automaton.json']:
  if p.exists():return json.loads(p.read_text())
 raise FileNotFoundError('run previous lattice.py')
def A(M,v):return tuple(sum(a*b for a,b in zip(M[4*i:4*i+4],v)) for i in range(4))
def hp(w):
 n=len(w);a=[0]*n
 for j,m in enumerate(w):
  for i in range(j):
   if j-i<=3 and m>>(j-i-1)&1:a[j]|=1<<i
   else:a[i]|=1<<j
 d=[[0]*n for _ in range(1<<n)]
 for i in range(n):d[1<<i][i]=1
 for S,row in enumerate(d):
  for i,c in enumerate(row):
   if not c:continue
   r=a[i]&~S
   while r:
    b=r&-r;r-=b;d[S|b][b.bit_length()-1]+=c
 return sum(d[-1])
def main(stream=False):
 d=load();f=d['readout'];v0=d['start'];Ms=d['matrices'];b=[4,2,2,1]
 # Matrix count and exact HNF
 need(len(Ms)==168 and len(set(map(tuple,Ms)))==136)
 H=lambda v:sum(x*y for x,y in zip(f,v))
 v,w=(A(Ms[d['words'].index(list(t))],v0) for t in [(7,6,7),(7,1,7)])
 row=lambda M:tuple(sum(f[k]*M[4*k+i] for k in range(4)) for i in range(4))
 need(hermite_normal_form(Matrix(list(map(row,Ms))).T)==diag(*b))
 need(all(f[j]%b[j]==0 for j in range(4)))
 need(all(M[4*i+j]*b[i]%b[j]==0 for M in Ms for i in range(4) for j in range(4)))
 if stream:
  print(len(Ms),*v0,*f)
  for M in Ms:print(*M)
  return
 seed=(0,0,3,5,2,4);S=(0,3,5,2,4)
 need(H(v)==H(w)==hp(seed+(7,6,7)+S)==hp(seed+(7,1,7)+S)==33423)
 U={H(A(M,v)) for M in Ms};V={H(A(M,w))for M in Ms}
 need(len(U)==len(V)==136 and not U&V)
 need(min((abs(a-b),a,b)for a in U for b in V)==(12,18060541,18060529))
 need({a%32 for a in U}=={b%32 for b in V}==set(range(1,32,2)))
 X={A(M,v0) for M in Ms};need(len(X)==136 and len(set(map(H,X)))==135)
 print('Exact 33423 fiber: 136 disjoint successors on each side; both full mod 32. Observable lattice diag(4,2,2,1).')
if __name__=='__main__':main('--stream' in sys.argv)
