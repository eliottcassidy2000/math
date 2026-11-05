import json,math,itertools,functools

def need(x,s):
 if not x:raise RuntimeError(s)
def graph(w):
 owners=[]
 for i,x in enumerate(w):owners.extend([i]*x)
 n=len(owners);a=[0]*n
 for i in range(n):
  for j in range(i+1,n):
   forward=owners[i]==owners[j] or (owners[j]-owners[i])%3==1
   u,v=(i,j) if forward else (j,i);a[u]|=1<<v
 return a
@functools.lru_cache(None)
def actual_sc(w):
 a=graph(w);n=len(a);degrees=[r.bit_count() for r in a];groups={d:[i for i,x in enumerate(degrees) if x==d] for d in set(degrees)}
 for d,vs in groups.items():
  if len(vs)!=len(groups.get(n-1-d,[])):return False
 keys=sorted(groups);choices=[list(itertools.permutations(groups[n-1-d])) for d in keys]
 for assignment in itertools.product(*choices):
  p=[None]*n
  for d,images in zip(keys,assignment):
   for v,img in zip(groups[d],images):p[v]=img
  if all(((a[i]>>j)&1)!=((a[p[i]]>>p[j])&1) for i in range(n) for j in range(i)):return True
 return False
def comps(d):
 for a in range(d+1):
  for b in range(d-a+1):yield a,b,d-a-b
def main():
 limit=200;inventory=[{} for _ in range(limit+1)];classes=0
 for n in range(3,limit+1):
  for x in range(1,n//3+1):
   for y in range(x,(n-x)//2+1):
    z=n-x-y
    if z<y:continue
    d=min(y-x,z-y);multiplicity=1 if d==0 else 2
    inventory[n][d]=inventory[n].get(d,0)+multiplicity;classes+=multiplicity
  U=(math.comb(n-1,2)+2*(n%3==0))//3
  need(sum(inventory[n].values())==U,'ordinary cyclic orbit total')
  need(max(inventory[n])==(n-3)//3,'exact maximum repair cost')
  for d in range((n-3)//3+2):
   formula=(n-1)//2 if d==0 else (2*((n-1-3*d)//2) if n>=3+3*d else 0)
   need(inventory[n].get(d,0)==formula,'exact sorted-partition cost coefficient')
 # Independent formal-series multiplication for the bivariate rational function.
 A={(0,0):1}
 for size,mark in ((1,0),(2,0),(3,1)):
  B=dict(A)
  for n in range(size,limit-2):
   for d in range((limit-3)//3+1):
    B[n,d]=B.get((n,d),0)+B.get((n-size,d-mark),0)
  A=B
 gf=0
 for n in range(3,limit+1):
  for d in range((n-3)//3+2):
   v=A.get((n-3,d),0)+A.get((n-6,d-1),0)
   need(v==inventory[n].get(d,0),'bivariate rational GF');gf+=1
 # Tail coefficients use a separate univariate coin-product recurrence.
 tail=[0]*(limit+1);tail[0]=1
 for step in (1,2,3):
  for i in range(step,limit+1):tail[i]+=tail[i-step]
 tails=0
 for n in range(3,limit+1):
  for d in range(1,(n-3)//3+2):
   actual=sum(v for j,v in inventory[n].items() if j>=d)
   N=n-3-3*d;expected=2*tail[N] if N>=0 else 0
   need(actual==expected,'sorted-gap tail shift');tails+=1
 moments=[0]*10001;moments[0]=1
 for step in (1,2,3,3):
  for i in range(step,len(moments)):moments[i]+=moments[i-step]
 for n in range(3,limit+1):
  need(sum(d*v for d,v in inventory[n].items())==(2*moments[n-6] if n>=6 else 0),'exact first moment GF')
 physical=searches=0
 for x in range(1,6):
  for y in range(x,6):
   for z in range(y,6):
    w=(x,y,z);expected=min(y-x,z-y);found=None
    for cost in range(expected+1):
     for inc in comps(cost):
      v=tuple(a+b for a,b in zip(w,inc));searches+=1
      if actual_sc(v):found=cost;break
     if found is not None:break
    need(found==expected,'literal physical minimum without weight-symmetry shortcut');physical+=1
 ratios=[]
 for n in (200,1000,10000):
  U=(math.comb(n-1,2)+2*(n%3==0))//3;moment=2*moments[n-6]
  ratios.append({'n':n,'classes':U,'sum_cost':moment,'mean_over_order':f'{moment}/{n*U}'})
 print(json.dumps({'status':'PASS','ordinary_classes':classes,'max_order':limit,'bivariate_coefficient_checks':gf,'tail_checks':tails,'physical_initial_graphs':physical,'physical_candidate_searches':searches,'distinct_actual_graphs_checked':actual_sc.cache_info().currsize,'exact_mean_ratios':ratios},sort_keys=True,indent=2))
if __name__=='__main__':main()
