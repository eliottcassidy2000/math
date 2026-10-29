import sys,json
sys.path.insert(0,str(__import__('pathlib').Path(__file__).resolve().parents[1]/'singularity'))
from actual_grammar import build
from native_kernel import observer,sectors
from itertools import combinations

def stars(w):
 A,c,h=build(w);p,M=observer(A);S=sectors(p,M)[1];pairs=list(combinations(range(len(A)),2));idx={p:i for i,p in enumerate(pairs)};cols=list(zip(*S));out=[]
 def col(u,k):
  return cols[idx[u,k]] if u<k else tuple(-x for x in cols[idx[k,u]])
 for u,v in combinations(range(len(A)),2):
  for k in range(len(A)):
   if k not in (u,v) and col(u,k)==col(v,k):out.append((u,v,k))
 return out
if __name__=='__main__':
 W=['ee','cbe']
 for n in range(3):
  for w in W:
   s=stars(w);A,c,h=build(w)
   print(w,'parent',h[-1]['parent'],'stars',s)
  W=[w+x for w in W for x in ('ab' if w.endswith(('a','c')) else 'cde')]
