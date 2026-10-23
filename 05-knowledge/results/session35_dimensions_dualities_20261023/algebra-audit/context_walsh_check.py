from itertools import product
import sympy as s
import json

def require(c,label):
 if not c:raise RuntimeError(label)
def binary_rank(values):
 basis={}
 for a in values:
  while a:
   p=a.bit_length()-1
   if p in basis:a^=basis[p]
   else:basis[p]=a;break
 return len(basis)
r=3;N=1<<r
for mask in range(1<<N):
 support=[a for a in range(N) if mask>>a&1]
 f=[sum((a+1)*(-1)**((a&x).bit_count()%2) for a in support) for x in range(N)]
 observed=[tuple(f[x^u] for u in range(N)) for x in range(N)]
 classes=len(set(observed))
 require(classes==2**binary_rank(support),('context quotient',mask))
 require(s.Matrix(observed).rank()==len(support),('translate rank',mask))
 for x in range(N):
  for y in range(N):
   require((observed[x]==observed[y])==all(((a&(x^y)).bit_count()%2)==0 for a in support),('annihilator',mask,x,y))
for v in range(3,10):
 pairs=[(1<<i)|(1<<j) for i in range(v) for j in range(i+1,v)]
 require(binary_rank(pairs)==v-1,('even weight span',v))
print(json.dumps({'all_supports_F2_3':256,'annihilator_pair_checks':256*64,'exact_linear_rank_checks':256,'even_weight_rank_checks':'v=3..9','q7_context_states':2**56,'q7_linear_span':1+57*56//2},indent=2))
