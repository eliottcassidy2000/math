from probe import make,hp
import json,time
init,step,val=make(2,[[0,1],[0,0]])
seen={}; counts=[0]*11

def dfs(v,w):
 h=val(v);counts[len(w)]+=1
 if h<=10000:seen.setdefault(h,w)
 if len(w)==10 or h>10000:return
 for m in range(4):dfs(step(v,m),w+(m,))
t=time.time();dfs(init,())
for h,w in sorted(seen.items()):
 if h%101==0 or h in (7477,27371,9999):assert hp(w)==h
with open('witnesses.json','w') as f:json.dump({str(h):w for h,w in sorted(seen.items())},f)
print(json.dumps({'counts':counts,'seconds':time.time()-t,'distinct':len(seen),'missing':[h for h in range(1,10001,2) if h not in seen],'controls':{h:seen.get(h) for h in [7477,9999]}},indent=2))
