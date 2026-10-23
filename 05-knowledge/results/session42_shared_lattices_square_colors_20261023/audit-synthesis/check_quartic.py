import itertools,json

def need(x,msg):
 if not x:raise RuntimeError(msg)
checks=0
for n in range(2,6):
 E=list(itertools.combinations(range(n),2))
 for mask in range(1<<len(E)):
  X=[[0]*n for _ in range(n)];deg=[0]*n
  for k,(a,b) in enumerate(E):
   X[a][b]=X[b][a]=-1 if mask>>k&1 else 1
   if mask>>k&1:deg[a]+=1;deg[b]+=1
  wedge=sum(X[v][a]*X[v][b] for v in range(n) for a,b in itertools.combinations([u for u in range(n) if u!=v],2))
  need(2*wedge==sum((n-1-2*d)**2-(n-1) for d in deg),'quadratic wedge')
  X2=[[sum(X[a][k]*X[k][b] for k in range(n)) for b in range(n)] for a in range(n)]
  trace4=sum(X2[a][b]*X2[b][a] for a in range(n) for b in range(n))
  cycles=0
  for a,b,c,d in itertools.combinations(range(n),4):
   for p in ((a,b,c,d),(a,b,d,c),(a,c,b,d)):
    val=1
    for k in range(4):val*=X[p[k]][p[(k+1)%4]]
    cycles+=val
  need(trace4==n*(n-1)*(2*n-3)+8*cycles,'quartic C4')
  checks+=2
print(json.dumps({'status':'PASS','exact_quadratic_quartic_checks':checks,'orders':[2,3,4,5],'matrix_scope':'zero diagonal, every offdiagonal entry plus or minus one'},indent=2))
