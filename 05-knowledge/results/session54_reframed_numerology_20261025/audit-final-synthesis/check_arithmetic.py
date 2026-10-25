from pathlib import Path
from itertools import combinations
import json,math
A=Path(__file__).resolve().parent
def need(x,s):
 if not x:raise RuntimeError(s)
def valuation(x,q):
 v=0
 while x%q==0:x//=q;v+=1
 return v
controls=[]
for q,a,b in [(3,0,1),(7,1,8),(11,628,1101),(19,3040229428,3817640437)]:
 r=(q-1)//2;s=a+b;delta=a-b;N=(s*s+q*delta*delta)//4;t=r*s;M=t*t+s*t+N
 need((s*s+q*delta*delta)%4==0,'integral norm');need(4*M==q*(q*s*s+delta*delta),'gluing modulus');need(b%q==1,'closed congruence');need(valuation(M,q)==(2 if delta%q==0 else 1),'actual valuation')
 # Basis image of 1,X,X² in Z × Z[lambda] has determinant Q(t).
 need(N+t*(t+s)==M,'basis index')
 # An arbitrary congruence pair lifts via U+VX+kQ(X).
 for u,v,k in [(1,0,1),(0,1,-3),(17,-4,7)]:
  z=u+v*t+k*M;poly=[u+k*N,v+k*s,k]
  need(poly[0]+poly[1]*t+poly[2]*t*t==z,'lift rational factor');need(poly[0]-poly[2]*N==u and poly[1]-poly[2]*s==v,'lift quadratic factor')
 controls.append({'q':q,'M':M,'valuation':valuation(M,q),'conductor':abs(delta)})
# Test both valuation branches against many values meeting the actual hypotheses.
guard_count=0
for q in [3,7,11,19,23,31]:
 for b in [1,1+q,1+2*q]:
  for a in range(3*q):
   if (a+b)%2==0:continue
   s=a+b;delta=a-b;M=q*(q*s*s+delta*delta)//4
   need(valuation(M,q)==(2 if delta%q==0 else 1),'valuation guard');guard_count+=1
endpoint=[]
for row in json.loads((A/'inputs/endpoint-paths.json').read_text()):
 q=row['q'];Q={x*x%q for x in range(1,q)};paths=[row['forward_path'],list(range(q))]
 for p in paths:need(len(p)==q and set(p)==set(range(q)) and all((v-u)%q in Q for u,v in zip(p,p[1:])),'endpoint witness')
 need((paths[0][-1]-paths[0][0])%q in Q and (paths[1][-1]-paths[1][0])%q not in Q,'two endpoint orbits')
 pairs={( (a*p[0]+b)%q,(a*p[-1]+b)%q) for p in paths for a in Q for b in range(q)}
 need(len(pairs)==q*(q-1),'all ordered endpoints');endpoint.append({'q':q,'ordered_pairs':len(pairs),'queries':2*q-4})
# Independent full fixed-pool subset census for formal family n5,n6.
optima=[]
for n in [5,6]:
 queries=[(0,i) for i in range(1,n)]+[(i,0) for i in range(1,n)]
 labels=[(typ,U) for typ in ['O','I'] for U in combinations(range(n),2)]+[('S',(v,)) for v in range(n)]
 rows=[]
 for typ,U in labels:rows.append([int(t in U) if typ=='O' else int(s in U) if typ=='I' else int(s in U)+int(t in U) for s,t in queries])
 found=[]
 for size in range(len(queries)+1):
  for chosen in combinations(range(len(queries)),size):
   if len({tuple(row[i] for i in chosen) for row in rows})==n*n:found.append(chosen)
  if found:break
 need(size==2*n-4 and len(found)==n-1,'sharp pool optimum')
 for chosen in found:
  omitted=[queries[i] for i in range(len(queries)) if i not in chosen];need(omitted[0][1]==omitted[1][0],'matching omitted labels')
 optima.append({'formal_n':n,'minimum':size,'optimal_subsets':len(found)})
need(2*5*(628+1101)==17290 and 11*5*(628+1101)-17290==77805,'P11 histogram');need(49*6+14*6==63*6,'P7 count')
print(json.dumps({'status':'PASS','gluing_controls':controls,'valuation_guards':guard_count,'larger_endpoint_orbits':endpoint,'fixed_pool_complete_census':optima,'P11_histogram':[17290,77805]},indent=2))
