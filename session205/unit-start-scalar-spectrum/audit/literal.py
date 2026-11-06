import itertools,json

def counts(out,root):
 n=len(out);full=(1<<n)-1
 dp=[[0]*n for _ in range(1<<n)]
 for v in range(n):dp[1<<v][v]=1
 for mask in range(1,1<<n):
  for v in range(n):
   val=dp[mask][v]
   if val:
    for w in range(n):
     if out[v]>>w&1 and not(mask>>w&1):dp[mask|1<<w][w]+=val
 H=[sum(x) for x in dp];H[0]=1
 st=[[0]*n for _ in range(1<<n)];st[1<<root][root]=1
 for mask in range(1,1<<n):
  for v in range(n):
   if st[mask][v]:
    for w in range(n):
     if out[v]>>w&1 and not(mask>>w&1):st[mask|1<<w][w]+=st[mask][v]
 R=sum(sum(st[m])*H[full^m] for m in range(1,1<<n))
 return (H[full],R,dp[full][root],H[full^(1<<root)]),sum(st[full])
checks=0
for l in range(6):
 for word in itertools.product('OZJ',repeat=l):
  out=[2,0];root=0;v=(1,2,0,1)
  for c in word:
   h,r,e,d=v;n=len(out);out.append(0)
   for i in range(n):
    if c=='O' or i!=root:out[i]|=1<<n
    else:out[n]|=1<<i
   if c=='O':v=(h,2*r,0,d)
   elif c=='Z':v=(h+r-e,2*r-d,d,d)
   else:v=(h+r-e,h+r,h-e,h);root=n
  actual,s=counts(out,root)
  if actual!=v or s!=1:raise RuntimeError((word,v,actual,s))
  checks+=1
print(json.dumps({'status':'PASS','literal_recipe_checks':checks,'max_vertices':7},indent=2))
