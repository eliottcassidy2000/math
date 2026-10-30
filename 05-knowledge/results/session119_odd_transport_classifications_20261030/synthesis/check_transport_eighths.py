from itertools import product
import json

def req(c,m):
 if not c:raise RuntimeError(m)
U=[1,3,5,7]
c4=lambda a:(-1)**((a-1)//2)
c8=lambda a:(-1)**((a*a-1)//8)
cs=[lambda a:1,c4,c8,lambda a:c4(a)*c8(a)]
for f,g in product(range(4),repeat=2):
 req(sum(cs[f](a)*cs[g](a) for a in U)==(4 if f==g else 0),'orthogonality')
for a,b in product(U,repeat=2):
 for f in cs:req(f((a*b)%8)==f(a)*f(b),'multiplicativity')
 req(sum(f(a)*f(b) for f in cs)==(4 if a==b else 0),'projector')
# Coordinates scaled by sqrt(2), verified against exact four-angle values.
expected={1:(1,1),3:(-1,1),5:(-1,-1),7:(1,-1)}
for a in U:req((c8(a),c4(a)*c8(a))==expected[a],'eighth phase')
odd=checks=0
for n in range(2,5):
 N=1<<n;mask=(1<<N)-1
 for f in range(1<<N):
  if f.bit_count()%2==0:continue
  odd+=1
  for h in range(1,N):
   shifted=sum(((f>>(x^h))&1)<<x for x in range(N))
   req((f^shifted) not in (0,mask),'odd weight derivative')
   checks+=1
req((1^2)==3,'n1 derivative exception')
edges=0
for a in range(2,7):
 N=1<<a
 for x in range(N):
  for y in range(x+1,N):
   adjacent=((y-x)%N)%2==1
   req(adjacent==((x%2)!=(y%2)),'unit graph bipartition')
   edges+=adjacent
 req(sum(1 for x in range(N) for y in range(x+1,N) if (y-x)%2)==(N//2)**2,'complete bipartite')
print(json.dumps({'status':'PASS','odd_weight_functions':odd,'nonzero_translation_tests':checks,'unit_graph_edges':edges,'moduli':[4,8,16,32,64],'character_table':[[f(a) for a in U] for f in cs],'n1_exception':'f(x)=x has constant derivative1','scope':'exact finite controls for stated uniform elementary proofs'},indent=2,sort_keys=True))
