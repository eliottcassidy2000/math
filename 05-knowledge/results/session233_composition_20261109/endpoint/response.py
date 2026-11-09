import sys

def require(condition, message):
 if not condition:
  raise RuntimeError(message)

if "--fail" in sys.argv:
 require(False, "Deliberate verification failure sentinel")

from probe import ends
def base(n):
 a=[sum(1<<j for j in range(i+1,n)) for i in range(n)]
 a[0]^=1<<(n-1); a[-1]|=1
 return a

def extend(n,N):
 a=base(n)+[0]
 for i in range(n):
  if 1<=i<=n-2 and N>>(i-1)&1:a[i]|=1<<n
  else:a[n]|=1<<i
 return a

def stats(a):
 n=len(a);full=(1<<n)-1;d=[[0]*n for _ in range(1<<n)]
 for v in range(n):d[1<<v][v]=1
 for s in range(1,1<<n):
  for v in range(n):
   if s>>v&1:
    t=s^(1<<v)
    if t:d[s][v]=sum(d[t][u] for u in range(n) if t>>u&1 and a[v]>>u&1)
 h=[sum(row) for row in d];h[0]=1
 B=[sum(d[s][v]*h[full^s] for s in range(1,1<<n)) for v in range(n)]
 return h[-1],B
for n in range(4,10):
 a=base(n);h,B=stats(a);print('base',n,'h',h,'B',B)
 for N in range(1,1<<(n-2)):
  u=extend(n,N);hu,Bu=stats(u)
  expected=h+B[0]+B[-1]+sum(B[i] for i in range(1,n-1) if not N>>(i-1)&1)
  require(Bu[-1]==expected, f"Response identity failed: n={n}, N={N}")
  if n==6:print('selector',N,hu,Bu[-1])
