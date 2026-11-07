"""Independent endpoint DP; no imported producer code or external data."""
import json,itertools

def need(x,msg):
 if not x:raise RuntimeError(msg)
def move(g,c):
 n=len(g)
 if c=='Q':return tuple([g[0]|(1<<n),*g[1:],(1<<n)-2])
 return tuple([g[0],*[r|(1<<n) for r in g[1:]],1])
def profile(g):
 n=len(g);M=1<<n;end=[[0]*n for _ in range(M)];start=[[0]*n for _ in range(M)]
 for v in range(n):end[1<<v][v]=start[1<<v][v]=1
 for mask in range(1,M):
  for v in range(n):
   if not(mask>>v&1):continue
   for w in range(n):
    if mask>>w&1:continue
    if g[v]>>w&1:end[mask|1<<w][w]+=end[mask][v]
    if g[w]>>v&1:start[mask|1<<w][w]+=start[mask][v]
 h=[sum(x) for x in end];h[0]=1;allmask=M-1
 L=sum(end[m][0]*h[allmask^m] for m in range(M))
 R=sum(start[m][0]*h[allmask^m] for m in range(M))
 D=allmask^1;b=sum(h[m]*h[D^m] for m in range(M) if not(m&1))
 return h[-1],start[-1][0],end[-1][0],L,R,h[D],b

def transfer(p,c):
 h,s,e,L,R,d,b=p
 if c=='Q':return h+L-s,d,e,2*L-d,R+b,d,2*b
 return h+R-e,s,d,L+b,2*R-d,d,2*b

def closed(p,q,z):
 h,s,e,L,R,d,b=p;need(s==e==d,'stable');a=L-d;r=R-d;u=2**q;v=2**z
 return(h+(u-1)*a+(v-1)*r+(u-1)*(v-1)*b,d,d,u*(a+(v-1)*b)+d,v*(r+(u-1)*b)+d,d,u*v*b)
controls=literal=0
for n in range(2,5):
 for code in range(1<<(n*(n-1)//2)):
  g=[0]*n;k=0
  for i in range(n):
   for j in range(i+1,n):
    if code>>k&1:g[i]|=1<<j
    else:g[j]|=1<<i
    k+=1
  p=profile(g)
  need(transfer(transfer(p,'Q'),'Z')==transfer(transfer(p,'Z'),'Q'),'commuting algebra')
  for c in 'QZ':need(profile(move(g,c))==transfer(p,c),'literal transfer');literal+=1
  stable=move(move(g,'Q'),'Z');ps=profile(stable)
  need(ps[1]==ps[2]==ps[5],'stable reset')
  for q in range(4):
   for z in range(4):
    v=ps
    for c in 'Q'*q+'Z'*z:v=transfer(v,c)
    need(v==closed(ps,q,z),'closed formula');controls+=1
    if q+z<=2:
     graph=stable
     for c in 'Z'*z+'Q'*q:graph=move(graph,c)
     need(profile(graph)==v,'independent literal closed formula');literal+=1
  # Full seven-state commutation survives nonstable starting coordinates.
  for word in itertools.product('QZ',repeat=5):
   v=p
   for c in word:v=transfer(v,c)
   w=p
   for c in 'Q'*word.count('Q')+'Z'*word.count('Z'):w=transfer(w,c)
   need(v==w,'all-word abelianization');controls+=1
print(json.dumps({'status':'PASS','algebra_controls':controls,'literal_DP_controls':literal,'source_universe':'All labeled rooted tournaments of orders 2..4, root 0','scope':'Finite checks of all-length formulas proved in COMMUTING.md'},indent=2))
