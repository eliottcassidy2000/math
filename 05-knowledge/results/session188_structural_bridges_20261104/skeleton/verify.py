"""Exact skeleton coefficients and direct Hamilton-path DP; standard library only."""
import itertools,math,json
SEEDS=((16,5,33,39,46,3),(16,5,49,39,42,3),(18,4,49,39,42,3))
def need(ok,msg):
 if not ok:raise RuntimeError(msg)
def skeleton(a,root):
 ext=[v for v in range(len(a)) if v!=root]; c=[0]*(len(ext)+2)
 for p in itertools.permutations(ext):
  good=[];bad=[]
  for g in range(len(p)+1):
   u=None if g==0 else p[g-1];v=None if g==len(p) else p[g]
   if (u is None or a[u]>>root&1) and (v is None or a[root]>>v&1):good.append(g)
   if u is not None and v is not None and not(a[u]>>v&1):bad.append(g)
  if set(bad)<=set(good):
   for k in range(len(good)+1):
    if 0<=k-len(bad)<=len(good)-len(bad):c[k]+=math.comb(len(good)-len(bad),k-len(bad))
 while c and not c[-1]:c.pop()
 return c
def inflate(a,r):
 t=[0]*(5+r)
 for u in range(5):
  t[u]=a[u]&31
  for v in range(5,5+r):
   if a[u]>>5&1:t[u]|=1<<v
   else:t[v]|=1<<u
 for u,v in itertools.combinations(range(5,5+r),2):t[u]|=1<<v
 return t
def hp(a):
 n=len(a);d=[[0]*n for _ in range(1<<n)]
 for i in range(n):d[1<<i][i]=1
 for m in range(1,1<<n):
  for v in range(n):
   if m>>v&1:
    z=m^(1<<v)
    if z:d[m][v]=sum(d[z][u] for u in range(n) if z>>u&1 and a[u]>>v&1)
 return sum(d[-1])
stir=[[0]*8 for _ in range(8)];stir[0][0]=1
for n in range(1,8):
 for k in range(1,n+1):stir[n][k]=stir[n-1][k-1]+k*stir[n-1][k]
rows=[]
for j,a in enumerate(SEEDS):
 b=skeleton(a,5); need(b==[9+2*j,27+2*j,21,4],('skeleton',j,b)); vals=[]
 for r in range(7):
  actual=hp(inflate(a,r));pred=sum(b[k]*math.factorial(k)*stir[r][k] for k in range(min(r+1,len(b))))
  need(actual==pred,('hp',j,r,actual,pred));vals.append(actual)
  if r:need(actual==4*3**r+9*2**r-3+2*j,('formula',j,r))
 rows.append({'state':j,'skeleton':b,'H_r0_to_6':vals})
print(json.dumps({'status':'PASS','rows':rows},indent=2))
