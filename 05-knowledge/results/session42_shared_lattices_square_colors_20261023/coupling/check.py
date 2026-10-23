import itertools,math,json
from fractions import Fraction as F
from collections import Counter

def need(x,m):
 if not x:raise RuntimeError(m)
def psd_ldl(M):
 A=[row[:] for row in M];positive=zero=0
 for k in range(len(A)):
  pivot=A[k][k]
  need(pivot>=0,'negative PSD pivot')
  if not pivot:
   need(all(A[k][j]==0 for j in range(k+1,len(A))),'zero pivot with nonzero row')
   zero+=1;continue
  positive+=1
  for i in range(k+1,len(A)):
   for j in range(i,len(A)):
    v=A[i][j]-A[i][k]*A[k][j]/pivot
    A[i][j]=A[j][i]=v
 return positive,zero

def run(n):
 edges=list(itertools.combinations(range(n),2));idx={e:k for k,e in enumerate(edges)};N=len(edges);size=1<<N;full=size-1
 perms=list(itertools.permutations(range(n)));maps=[]
 for p in perms:
  pos={v:i for i,v in enumerate(p)}
  mp=[idx[tuple(sorted((p[a],p[b])))] for a,b in edges]
  inversion=sum(1<<k for k,(a,b) in enumerate(edges) if pos[a]>pos[b])
  maps.append((mp,inversion))
 gc={};tc={}
 for x in range(size):
  plain=[sum(1<<mp[k] for k in range(N) if x>>k&1) for mp,_ in maps]
  gc[x]=min(plain);tc[x]=min(v^s for v,(_,s) in zip(plain,maps))
 # Exact signed inversion averages, one representative per ordinary graph shape.
 mus={h:F(sum((-1)**((h & mask).bit_count()) for _,mask in maps),math.factorial(n)) for h in set(gc.values())}
 max_nonconstant=max((abs(mu) for h,mu in mus.items() if h),default=F(0))
 need(max_nonconstant<=F(1,3),'uniform third bound')
 if n>=3:need(max_nonconstant==F(1,3),'finite sharp value')
 cubic_star_control=None
 if n==4:
  def cubic_star(x):
   score=[0]*n
   for bit,(a,b) in enumerate(edges):score[b if x>>bit&1 else a]+=1
   return sum((-1)**(n-1-d) for d in score)
  for g in set(gc.values()):need(sum(cubic_star(x) for x in range(size) if gc[x]==g)==0,'star observable not erased')
  for t in set(tc.values()):need(len({cubic_star(x) for x in range(size) if tc[x]==t})==1,'star not tournament invariant')
  cubic_star_control=sorted(set(cubic_star(x) for x in range(size)))
 reports=[]
 for paired in (False,True):
  G=lambda x:min(gc[x],gc[full^x]) if paired else gc[x]
  T=lambda x:min(tc[x],tc[full^x]) if paired else tc[x]
  cells=Counter((G(x),T(x)) for x in range(size));gs=sorted({a for a,b in cells});ts=sorted({b for a,b in cells})
  weights=[[cells[a,b] for b in ts] for a in gs]
  gm=[sum(row) for row in weights];tm=[sum(row[j] for row in weights) for j in range(len(ts))]
  K=[[sum((F(weights[i][j]*weights[k][j],gm[i]*tm[j]) for j in range(len(ts))),F(0)) for k in range(len(gs))] for i in range(len(gs))]
  need(all(sum(row)==1 for row in K),'not stochastic')
  need(all(K[i][i]>0 for i in range(len(gs))),'not positive diagonal')
  need(all(gm[i]*K[i][k]==gm[k]*K[k][i] for i in range(len(gs)) for k in range(len(gs))),'not reversible')
  todo=[0];seen={0}
  while todo:
   i=todo.pop()
   for k,p in enumerate(K[i]):
    if p and k not in seen:seen.add(k);todo.append(k)
  need(len(seen)==len(gs),'not irreducible')
  # Exact rational positivity of D_g K and spectral-gap certificate.
  positive_kernel=[[gm[i]*K[i][k] for k in range(len(gs))] for i in range(len(gs))]
  rank,zeros=psd_ldl(positive_kernel)
  gap=F(1) if n==2 else F(8,9)
  certificate=[[F(gm[i]*(i==k))-gm[i]*K[i][k]-gap*(F(gm[i]*(i==k))-F(gm[i]*gm[k],size)) for k in range(len(gs))] for i in range(len(gs))]
  pivots=psd_ldl(certificate)
  # Walsh predicted positive spectrum, with multiplicity, and trace-power checks.
  predicted=[mus[h]**2 for h in sorted(mus) if not paired or h.bit_count()%2==0]
  need(len(predicted)==len(gs),'Walsh dimension')
  for power in (1,2):
   trace=sum(K[i][i] for i in range(len(gs))) if power==1 else sum(K[i][j]*K[j][i] for i in range(len(gs)) for j in range(len(gs)))
   need(trace==sum(v**power for v in predicted),'Walsh trace identity')
  positive_spectrum=Counter(str(x) for x in predicted if x)
  # Independently verify every predicted Walsh eigenvector on graph quotient states.
  eigenvector_checks=0
  for h,mu in mus.items():
   if paired and h.bit_count()%2:continue
   orbit=[x for x in range(size) if gc[x]==h]
   values=[sum((-1)**((f&x).bit_count()) for f in orbit) for x in gs]
   for i in range(len(gs)):
    need(sum(K[i][j]*values[j] for j in range(len(gs)))==mu*mu*values[i],'exact Walsh eigenvector')
    eigenvector_checks+=1
  reports.append({'n':n,'complement_converse_paired':paired,'graph_states':len(gs),'tournament_states':len(ts),'incidence_edges':len(cells),'exact_gap_lower_bound':str(gap),'factorial_gap_lower_bound':str(F(4,n*math.factorial(n))),'exact_max_nonconstant_signed_average':str(max_nonconstant),'exact_positive_spectrum':dict(sorted(positive_spectrum.items())),'exact_positive_kernel_rank':rank,'invisible_tournament_dimensions':len(ts)-rank,'invisible_graph_dimensions':len(gs)-rank,'order_four_cubic_star_values':cubic_star_control,'gap_certificate_positive_pivots':pivots[0],'gap_certificate_zero_pivots':pivots[1],'exact_eigenvector_equations':eigenvector_checks})
 return reports
print(json.dumps({'status':'PASS','controls':[r for n in range(2,6) for r in run(n)],'scope':'exact finite stochastic, reversibility, connectivity, positive-semidefinite and gap checks; all-size bound proved separately'},indent=2,sort_keys=True))
