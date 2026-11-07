import json
from gap import rows,coeff

def matching(a,v):
 A=[i for i in range(len(a)) if a[v]>>i&1]
 B=[i for i in range(len(a)) if a[i]>>v&1]
 def go(k,used):
  if k==len(A):return 0
  ans=go(k+1,used)
  for b in B:
   if not(used>>b&1) and a[A[k]]>>b&1:ans=max(ans,1+go(k+1,used|(1<<b)))
  return ans
 return go(0,0)
checks=0;hist={}
for n in range(1,6):
 for bits in range(1<<(n*(n-1)//2)):
  a=rows(n,bits)
  for v in range(n):
   nu=matching(a,v);c=coeff(a,v)
   if any(x<=0 for x in c[:nu+2]) or any(c[nu+2:]):raise RuntimeError((a,v,nu,c))
   checks+=1;hist[nu]=hist.get(nu,0)+1
print(json.dumps({'status':'PASS','rooted_labeled_checks':checks,'matching_histogram':hist},indent=2))
