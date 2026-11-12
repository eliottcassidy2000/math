from pathlib import Path
import sys,itertools,json
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'session275_exact_selection_20261111/switches'))
from carrier import hp,profile

def cyclic(n):return tuple(sum(1<<((i+j)%n) for j in range(1,(n+1)//2)) for i in range(n))
def apex(C,A):
 t=[sum(1<<(i+1) for i in range(len(C)) if i not in A)]
 t += [(r<<1)|(i in A) for i,r in enumerate(C)]
 return tuple(t)
def covers(T):
 h=hp(T); f=[None]*len(h);f[0]=(1,)
 for S in range(1,len(h)):
  r=[0]*(S.bit_count()+1); b=S&-S;sub=S
  while sub:
   if sub&b:
    for k,x in enumerate(f[S^sub]):r[k+1]+=h[sub]*x
   sub=(sub-1)&S
  f[S]=tuple(r)
 return f[-1]
for n in (3,5,7,9):
 C=cyclic(n);groups={}
 for A in itertools.combinations(range(n),2):
  T=apex(C,A);p=profile(T);c=covers(T);groups.setdefault((p['H'],c),[]).append(A)
 print(n,[(x,len(y)) for x,y in groups.items()],flush=True)
