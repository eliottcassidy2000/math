"""Native regular hosts, exact coset bounds and integer transfer clocks."""
import itertools as it, hashlib,json

def need(ok,msg):
 if not ok:raise RuntimeError(msg)

def host(p,t,E):
 D={(e+1)%p for e in E};n=p*t;a=[[0]*n for _ in range(n)]
 for i in range(t):
  for u in range(p):
   for v in range(p):a[i*p+u][i*p+v]=int((v-u)%p in range(1,(p+1)//2))
 for i in range(t):
  for j in range(i+1,t):
   outer=(j-i)%t in range(1,(t+1)//2)
   for u in range(p):
    for v in range(p):
     bit=int((v-u)%p in D) if outer else int((u-v)%p not in D)
     a[i*p+u][j*p+v]=bit;a[j*p+v][i*p+u]=1-bit
 need(all(sum(row)==(n-1)//2 for row in a),'regular host')
 accepted=0;fixed=[0]*p;free=[0]*p;digest=hashlib.sha256();witness=None
 for phase in it.product(range(p),repeat=t):
  path=tuple(i*p+(phase[i]+u)%p for i in range(t) for u in range(p))
  valid=all(a[x][y] for x,y in zip(path,path[1:]))
  predicted=all((phase[i+1]-phase[i])%p in E for i in range(t-1))
  need(valid==predicted,'literal native path predicate')
  if valid:
   accepted+=1;delta=(phase[-1]-phase[0])%p;free[delta]+=1
   if phase[0]==0:fixed[delta]+=1
   digest.update(bytes(path));witness=path
 need(accepted==p*len(E)**(t-1),'path count')
 need(free==[p*x for x in fixed],'free versus pinned initial owner')
 dp=[1]+[0]*(p-1)
 for _ in range(t-1):dp=[sum(dp[(j-e)%p] for e in E) for j in range(p)]
 need(dp==fixed,'shared-state transfer')
 if p==3 and E=={0,1}:
  r=t-1;clock=[2,1,-1,-2,-1,1]
  need(fixed==[(2**r+clock[(r-2*d)%6])//3 for d in range(3)],'period-six formula')
 return {'p':p,'cells':t,'vertices':n,'outdegree':(n-1)//2,'E':sorted(E),'family_paths':accepted,'minimum_affine_cosets':len(E)**(t-1),'fixed_start_endpoint_counts':fixed,'free_start_difference_counts':free,'literal_path_sha256':digest.hexdigest(),'last_witness':witness}

def exhaustive_cosets(p,E):
 # Construct all additive subspaces of F_p^3 directly, no assumed classification.
 vs=list(it.product(range(p),repeat=3));zero=(0,0,0);subs={frozenset([zero])};todo=list(subs)
 while todo:
  H=todo.pop()
  for v in vs:
   if v in H:continue
   K=frozenset(tuple((h[j]+c*v[j])%p for j in range(3)) for h in H for c in range(p))
   if K not in subs:subs.add(K);todo.append(K)
 S={v for v in vs if (v[1]-v[0])%p in E and (v[2]-v[1])%p in E}
 cosets=set();contained=[]
 for H in subs:
  for v in vs:
   C=frozenset(tuple((h[j]+v[j])%p for j in range(3)) for h in H)
   cosets.add(C)
 for C in cosets:
  if C<=S:contained.append(C)
 max_size=max(map(len,contained))
 need(max_size==p,'maximum allowed affine-coset size')
 fibres={}
 for v in S:fibres.setdefault(((v[1]-v[0])%p,(v[2]-v[1])%p),set()).add(v)
 need(len(fibres)==len(E)**2 and all(frozenset(F) in cosets and len(F)==p for F in fibres.values()),'optimal diagonal cover')
 return {'p':p,'dimension':3,'E':sorted(E),'all_subspaces':len(subs),'all_distinct_affine_cosets':len(cosets),'contained_cosets':len(contained),'maximum_contained_size':max_size,'minimum_cover':len(fibres)}

# The formula applies to every length, including partial forests before odd-cell closure.
clock=[2,1,-1,-2,-1,1];row=[1,0,0];clock_rows=[]
for r in range(31):
 predicted=[(2**r+clock[(r-2*d)%6])//3 for d in range(3)]
 need(row==predicted and all((2**r+clock[(r-2*d)%6])%3==0 for d in range(3)),'integer clock')
 clock_rows.append({'r':r,'fixed_start_counts':row,'correction':[clock[(r-2*d)%6] for d in range(3)]})
 row=[row[d]+row[(d-1)%3] for d in range(3)]
result={'status':'PASS','scope':'specified contiguous cyclic-cut path families in actual regular tournaments; not all Hamilton paths','hosts':[host(3,t,{0,1}) for t in (3,5,7,9)]+[host(5,3,E) for E in ({0,1},{0,1,3})],'coset_audits':[exhaustive_cosets(3,{0,1}),exhaustive_cosets(5,{0,1,3})],'period_six_clock':clock_rows}
print(json.dumps(result,sort_keys=True,indent=2))
