#!/usr/bin/env python3
import json

def closure(coins,N):
    f=[0]+[-10**9]*N
    for n in range(1,N+1):
        f[n]=max(f[n-l]+p for l,p in coins if l<=n)
    return f

def guard(ok):
    if not ok: raise RuntimeError('profile obstruction check failed')
f=closure([(1,0),(2,19),(5,50)],500)
g=closure([(1,0),(4,39),(5,50)],500)
h=closure([(1,0),(2,19),(4,39),(5,50)],500)
def bar(x,n): return max(x[n+5*t]-50*t for t in range((len(x)-1-n)//5+1))
guard([bar(f,0),bar(g,0),bar(f,1),bar(g,1),bar(h,1)]==[0,0,7,6,8])
for n in range(101):
    guard(h[n]==max(f[i]+g[n-i] for i in range(n+1)))
print(json.dumps({'status':'PASS','bar_f_at_1':7,'bar_g_at_1':6,'tensor_of_bars_at_1':7,'bar_of_tensor_at_1':8,'finite_controls_only':True},sort_keys=True))
