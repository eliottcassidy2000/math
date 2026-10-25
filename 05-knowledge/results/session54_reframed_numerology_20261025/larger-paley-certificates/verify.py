#!/usr/bin/env python3
import pathlib,json
D=pathlib.Path(__file__).resolve().parent
P=101
def ck(p,m):
 if not p:raise RuntimeError(m)
def det_mod(A):
 A=[[x%P for x in row] for row in A];n=len(A);det=1
 for k in range(n):
  p=next((i for i in range(k,n) if A[i][k]),None)
  if p is None:return 0
  if p!=k:A[p],A[k]=A[k],A[p];det=-det
  pivot=A[k][k];det=det*pivot%P;inv=pow(pivot,-1,P)
  for i in range(k+1,n):
   f=A[i][k]*inv%P
   if f:
    for j in range(k+1,n):A[i][j]=(A[i][j]-f*A[k][j])%P
    A[i][k]=0
 return det%P

def rank_mod(A):
 basis={}
 for row in A:
  w=[x%P for x in row]
  for p,v in sorted(basis.items()):
   f=w[p]
   if f:w=[(x-f*y)%P for x,y in zip(w,v)]
  p=next((i for i,x in enumerate(w) if x),None)
  if p is not None:basis[p]=[x*pow(w[p],-1,P)%P for x in w]
 return len(basis)
results=[]
for n in (11,19,23,31):
 cert=json.loads((D/('inputs/P11-rank.json' if n==11 else f'certificate-{n}.json')).read_text());r=(n-1)//2;lam=(n-3)//4;V=set(range(n));Q={x*x%n for x in range(1,n)};E=[(u,v) for u in range(n) for v in range(n) if (v-u)%n in Q];EI={e:i for i,e in enumerate(E)};target=len(E)-2*n+1
 ck(len(Q)==r and Q.isdisjoint({(-x)%n for x in Q}),'Paley orientation')
 IN={v:{u for u,w in E if w==v} for v in V}
 ck(all(len(IN[v])==r for v in V),'regular degree')
 ck(all(len(IN[u]&IN[v])==lam for u in V for v in V if u!=v),'double regularity')
 ck(all(a*b%n in Q for a in Q for b in Q),'multiplier closure')
 ck({(b,(a+b)%n) for a in Q for b in V}==set(E),'arc transitivity')
 graph={i:[] for i in range(2*n)}
 for u,v in E:graph[u].append(n+v);graph[n+v].append(u)
 seen={0};todo=[0];tree=[]
 while todo:
  u=todo.pop(0)
  for v in graph[u]:
   if v not in seen:seen.add(v);todo.append(v);tree.append([u,v])
 ck(len(seen)==2*n,'connected support')
 def cyclevec(C):
  ck(len(C)==n and set(C)==V,'cycle vertex set');edges=list(zip(C,C[1:]+C[:1]));ck(all(e in EI for e in edges),'cycle arcs')
  return [int(e in edges) for e in E]
 if n!=11:
  ck(len(cert['seed_cycles'])==1,'first seed orbit claim')
  seed=cert['seed_cycles'][0]
  for C in cert['independent_cycles']:
   a=(C[1]-C[0])*pow((seed[1]-seed[0])%n,-1,n)%n;b=(C[0]-a*seed[0])%n
   ck(a in Q and C==[(a*x+b)%n for x in seed],'affine seed orbit membership')
 base=cyclevec(cert['base_cycle']);cycle_vectors=[cyclevec(C) for C in cert['independent_cycles']];A=[[x-y for x,y in zip(C,base)] for C in cycle_vectors];piv=cert['pivot_columns']
 ck(len(A)==target and len(piv)==target and len(set(piv))==target and all(0<=p<len(E) for p in piv),'rank dimensions')
 determinant=det_mod([[row[j] for j in piv] for row in A]);ck(determinant!=0,'full rank minor modulo101')
 endpoint=cert['endpoint_paths'][0] if n==11 else cert['endpoint_path']
 ck(len(endpoint)==n and set(endpoint)==V and endpoint[0]==0 and endpoint[-1]==1 and all(e in EI for e in zip(endpoint,endpoint[1:])),'forward endpoint HP')
 # Constant-difference cycles are a hostile, NOT sufficient rank witnesses.
 constant=[cyclevec([k*d%n for k in range(n)]) for d in sorted(Q)]
 constant_rank=rank_mod([[x-y for x,y in zip(C,constant[0])] for C in constant[1:]])
 ck(constant_rank==r-1 and constant_rank<target,'constant-cycle rank hostile')
 # Check general cut inequalities numerically for the entire relevant finite range too.
 ck(all(r*(r-2)+lam*(t-1)*(2*r-1-t)>0 for t in range(1,2*r)),'deficiency one expansion')
 ck(all(r*(r-3)+lam*(t-1)*(2*r-2-t)>0 for t in range(1,2*r-1)),'deficiency two expansion')
 results.append({'n':n,'arcs':len(E),'r':r,'lambda':lam,'cycle_difference_rank':target,'minor_determinant_mod101':determinant,'listed_Hamilton_cycles':len(A)+1,'support_tree':tree,'forward_endpoint_path':endpoint,'constant_difference_cycle_orbit_rank':constant_rank,'criterion_conclusion':{'beta_HP':n-1,'minimum_blocker_count':n*n,'two_source':n*(n-1)//2,'two_sink':n*(n-1)//2,'stars':n},'scope':'Exact application of sufficient criterion, not blocker enumeration or all-order inference'})
print(json.dumps({'status':'PASS','modulus':101,'results':results},indent=2,sort_keys=True))
