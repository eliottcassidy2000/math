"""Exhaustive finite regression for the direct-readiness classification."""
import json
M=50000
cur={(1,2,0,1)};ready=[set(),set()];n=2;states=0
while cur:
 states+=len(cur)
 for h,r,e,d in cur:
  C=r-h+e-d
  if not (r>=h>=d+e and d>=e>=0 and d%2 and h%2 and C>=0 and (r-e)%2==0):raise RuntimeError('invariant')
  if r-e==h-1:ready[0].add(h)
  if r-e==h+1:ready[1].add(h)
 nxt=set()
 for h,r,e,d in cur:
  if h+r-e<=M:
   nxt.add((h+r-e,2*r-d,d,d));nxt.add((h+r-e,h+r,h-e,h))
  if h+2*r<=M:nxt.add((h,2*r,0,d))
 cur=nxt;n+=1
limit=M//2
pred=[{2**m+1 for m in range(1,30)}|{3*2**m+3 for m in range(1,30)}, {1}|{3*2**t-1 for t in range(1,30)}|{5*2**t-1 for t in range(1,30)}|{3*(2**m+1)*2**t-1 for m in range(1,30) for t in range(1,30)}]
for i in range(2):
 a={h for h in ready[i] if h<limit};b={h for h in pred[i] if h<limit}
 if a!=b:raise RuntimeError((i,sorted(a-b),sorted(b-a)))
print(json.dumps({'status':'PASS','bound':M,'ready_parent_exclusive_bound':limit,'final_empty_order':n,'states_by_order_total':states,'minus_ready':sorted(h for h in ready[0] if h<limit),'plus_ready':sorted(h for h in ready[1] if h<limit)},indent=2))
