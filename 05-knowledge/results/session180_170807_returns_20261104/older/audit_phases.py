import itertools,json,pathlib
P=pathlib.Path(__file__).parent.parent/'phases'
def need(x,msg):
 if not x:raise RuntimeError(msg)
def perm(rows,w):
 h=f=0;n=len(rows)
 for p in itertools.permutations(range(n)):
  bad=[j for j in range(n-1) if not rows[p[j]]>>p[j+1]&1]
  if not bad:h+=1
  k=p.index(w)+1
  if not bad or bad==[k-1]:f+=1
 return h,f
def birth(rows,ins):
 n=len(rows);return [a|(((ins>>j)&1)<<n) for j,a in enumerate(rows)]+[((1<<n)-1)^ins]
def strong(r):
 n=len(r)
 for s in range(n):
  reached=1<<s
  while True:
   new=reached
   for v in range(n):
    if reached>>v&1:new|=r[v]
   if new==reached:break
   reached=new
  if reached!=(1<<n)-1:return False
 return True
def endpoint(rows,w):
 n=len(rows);N=1<<n;dp=[{} for _ in range(N)];h=[0]*N;h[0]=1
 for s in range(1,N):
  for v in range(n):
   if s>>v&1:
    t=s^(1<<v);dp[s][v]=sum(a for u,a in dp[t].items() if rows[u]>>v&1) if t else 1
  h[s]=sum(dp[s].values())
 return h[-1],sum(dp[s].get(w,0)*h[(N-1)^s] for s in range(N))
bank=json.loads((P/'small_suffix_search.json').read_text())
for z in bank['bank']:need(perm(z['rows'],z['port'])==(z['h'],z['f']),'bank literal mismatch')
pairs={(z['h'],z['f']) for z in bank['bank']}
need(len(pairs)==180,'bank cardinality')
need({str(d):sum((3*h+d,3*f) in pairs for h,f in pairs) for d in [2,4]}=={'2':6,'4':7},'affine pair count')
B=[8,1,3,6];need(perm(B,0)==(5,7),'base pair');receipts=[]
for s in range(1,5):
 root=4+s;source=B[:]+[0]*(s+1)
 for i in range(s):source[4+i]=15|sum(1<<(4+j) for j in range(i+1,s))
 source[0]|=1<<root;source[root]=((1<<root)-1)^1
 n=len(source);X=birth(source,((1<<n)-1)^(1<<1))
 for d in [2,4]:
  ins=((1<<n)-1) if d==2 else (1|sum(1<<j for j in range(4,n)))
  Y=birth(X,ins)
  need(strong(source) and strong(X) and strong(Y),'strongness')
  need([a&((1<<n)-1) for a in Y[:n]]==source,'old induced labels')
  need([j for j in range(len(Y)) if Y[j]>>root&1]==[0],'root unique inport')
  oldh,oldf=endpoint(source,root);h,f=endpoint(Y,root)
  need(oldh==oldf==5+7*2**s and h==f==3*oldh+d,'fixed label all-clock sample')
  if s<=2:need(perm(Y,root)==(h,f),'literal full permutation')
  receipts.append({'s':s,'d':d,'root':root,'newborns':[n,n+1],'h':h,'f':f})
# The stated failure of recipe recurrence on the actual displayed C4.
C4=[56,1,19,22,2,30];fail=[]
for a in range(6):
 if C4[a]&1:
  X=birth(C4,63^(1<<a))
  for d in [2,4]:
   Y=birth(X,63 if d==2 else 1);h,f=endpoint(Y,0)
   need((h,f)!=(57+d,63),'unexpected recurrence hit');fail.append({'a':a,'d':d,'h':h,'f':f})
out={'status':'PASS','bank_witnesses_permutation_checked':180,'named_full_birth_checks':receipts,'actual_C4_recipe_failures':fail,'scope':'finite bank witnesses and mathematical all-clock proof checked; original complete labelled census replay separate'}
(pathlib.Path(__file__).parent/'phase_audit.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
