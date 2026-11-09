def ends(a):
 n=len(a); d=[[0]*n for _ in range(1<<n)]
 for v in range(n): d[1<<v][v]=1
 for s in range(1,1<<n):
  for v in range(n):
   if s>>v&1:
    t=s^(1<<v)
    if t: d[s][v]=sum(d[t][u] for u in range(n) if t>>u&1 and a[u]>>v&1)
 return d[-1]
