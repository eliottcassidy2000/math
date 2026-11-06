import itertools,json

def require(p,msg):
 if not p:raise AssertionError(msg)
def rmat(n):return [[int(j==i+1 or i>=j+2) for j in range(n)]for i in range(n)]
def target(a,b):
 n=a+b+3;t=rmat(n);t[a][a+2]=1;t[a+2][a]=0;return t
def hp(t):
 n=len(t);dp=[[0]*n for _ in range(1<<n)]
 for i in range(n):dp[1<<i][i]=1
 for s in range(1,1<<n):
  for v in range(n):
   if s>>v&1:dp[s][v]+=sum(dp[s^(1<<v)][u] for u in range(n) if s>>u&1 and t[u][v])
 return sum(dp[-1])
def allcycles(t):
 n=len(t);out=set()
 for i in range(n):
  for k in range(2,n-i):
   for tail in itertools.permutations(range(i+1,n),k):
    cyc=(i,)+tail
    if all(t[cyc[z]][cyc[(z+1)%len(cyc)]] for z in range(len(cyc))):out.add(cyc)
 return out
seed=[[1,5,9],[5,13,23],[9,23,41]]
g={}
for a in range(21):
 for b in range(21):
  if a<3 and b<3:g[a,b]=seed[a][b]
  elif a>=3:g[a,b]=sum(g[a-k,b]for k in (1,2,3))
  else:g[a,b]=sum(g[a,b-k]for k in (1,2,3))
  if b>=3:require(g[a,b]==sum(g[a,b-k]for k in (1,2,3)),'commuting recurrences')
h=[1,1,1]
for n in range(3,45):h.append(sum(h[-3:]))
rows=[];cycle_cases=0
for n in range(3,14):
 for a in range(n-2):
  b=n-a-3;t=target(a,b);actual=hp(t);require(actual==g[a,b],('DP',a,b))
  conv=h[n]-2*h[a]*h[b]
  conv+=2*sum(h[a-p]*h[b-q]for p in range(a+1)for q in range(b+1)if (p+q)%2)
  conv+=2*sum(h[a-p]*h[b-q]for p in range(1,a+1)for q in range(1,b+1)if not(p+q)%2)
  require(conv==actual,('cycle convolution',a,b))
  cells={}
  for v in range(n):cells.setdefault(sum(t[v]),[]).append(v)
  for cell in cells.values():
   require(all(not(t[i][j]and t[j][k]and t[k][i])for i,j,k in itertools.permutations(cell,3)),('nontransitive score cell',a,b,cell))
  if n<=8:
   old=allcycles(rmat(n));new=allcycles(t)
   removed={(a,a+1,a+2)};added=set()
   for i in range(a+1):
    for j in range(a+2,n):
     if j-i>=3:added.add(tuple(range(i,a+1))+tuple(range(a+2,j+1)))
   for i in range(a):
    for j in range(a+3,n):added.add(tuple(range(i,a+1))+tuple(range(a+2,j+1))+(a+1,))
   require(old-new==removed,('removed',a,b));require(new-old==added,('added',a,b));cycle_cases+=1
  rows.append([n,a,b,actual])
# Verify rational numerator by clearing both recurrence denominators.
N={(0,0):1,(1,0):4,(0,1):4,(2,0):3,(0,2):3,(1,1):4,(2,1):2,(1,2):2,(2,2):1}
d=[1,-1,-1,-1]
for a in range(21):
 for b in range(21):
  coeff=sum(d[i]*d[j]*g[a-i,b-j]for i in range(min(3,a)+1)for j in range(min(3,b)+1))
  require(coeff==N.get((a,b),0),('rational numerator',a,b))
print(json.dumps({'status':'PASS','dp_cases':len(rows),'direct_cycle_inventory_cases':cycle_cases,'recurrence_grid':'0<=a,b<=20','rigidity_control':'every score cell transitive in all DP cases','rows':rows},indent=2))
