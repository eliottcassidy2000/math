"""Additional independent source-free width, GF(2), and mesh controls."""
from kernel import *
from functools import reduce
import json,hashlib,pathlib
q=tournament(5,10)
a=substitute(q,[(0,),tournament(2,1),(0,),(0,),(0,)])
require(a==(40,5,1,6,15,30) and h_literal(a)==19,'width5 H19')
C=(2,4,1);L=(2,0)
models=[tournament(4,63),ordinal(C,(0,)),ordinal((0,),C),substitute(C,[L,(0,),(0,)])]
counts=[0]*4
for code in range(64):
 q=tournament(4,code);found=[]
 for k,m in enumerate(models):
  if any(induced(q,p)==m for p in permutations(range(4))):found.append(k)
 require(len(found)==1,('Q4 iso classes',code));counts[found[0]]+=1
require(counts==[24,8,8,24],'Q4 multiplicities')

def rank(columns):
 span={0}
 for col in columns:span|={v^col for v in list(span)}
 return len(span).bit_length()-1

def rank_powers(a):
 n=len(a)
 columns=tuple(sum((a[i]>>j&1)<<i for i in range(n))for j in range(n))
 vectors=tuple(1<<i for i in range(n));ranks=[]
 for _ in range(n+1):
  ranks.append(rank(vectors))
  vectors=tuple(reduce(int.__xor__,[columns[i]for i in range(n)if v>>i&1],0)for v in vectors)
 return ranks

b=tournament(6,280);ranks=rank_powers(b)
require(ranks==[6,5,4,3,2,1,0] and h_literal(b)==41,'root nilpotent hostile')
trans=tournament(6,32767)
require(rank_powers(trans)==ranks and h_literal(trans)==1,'transitive lens control')
require(67**137<65**138 and 67**138>=65**139,'m0 example')
residuals=[(3,3),(5,4),(9,5),(17,6),(35,7),(63,8)]
for r,s in residuals:require(80267**(s-1)<r**22,('compression residual',r,s))

j={'status':'PASS','width5_H19_rows':a,'Q4_isomorphism_multiplicities':counts,'root_nilpotent_H41_rows':b,'root_GF2_rank_powers':ranks,'cofinality_65_67_m0':138,'compression_residual_inequalities':residuals,'script_sha256':hashlib.sha256(pathlib.Path(__file__).read_bytes()).hexdigest()}
print(json.dumps(j,indent=2,sort_keys=True))
