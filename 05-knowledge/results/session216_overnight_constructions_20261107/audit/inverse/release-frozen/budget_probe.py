from probe import make
import json,time
init,step,val=make(2,[[0,1],[0,0]])
seen=[set() for _ in range(6)]; costs={}; witnesses={}
def dfs(v,w,a,b):
 h=val(v)
 if h>10000:return
 cost=min(a,b);old=costs.get(h,99)
 if cost<old:costs[h]=cost;witnesses[h]=w
 if len(w)==10:return
 for m in range(4):dfs(step(v,m),w+(m,),a+(m&1),b+((m>>1)&1))
dfs(init,(),0,0)
print(json.dumps({'missing_by_exit_budget':{q:[h for h in range(23,10001,2) if costs.get(h,99)>q] for q in range(1,4)},'max_min_cost':max(costs.values())},indent=2))
