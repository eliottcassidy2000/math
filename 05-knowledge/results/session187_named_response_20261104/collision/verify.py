import itertools,json,math
from collections import defaultdict

def need(x,msg):
 if not x:raise RuntimeError(msg)
def hp(rows):
 n=len(rows);dp=[[0]*n for _ in range(1<<n)]
 for u in range(n):dp[1<<u][u]=1
 for s in range(1,1<<n):
  for u in range(n):
   if not dp[s][u]:continue
   b=rows[u]&~s
   while b:
    v=b&-b;b-=v;dp[s|v][v.bit_length()-1]+=dp[s][u]
 return sum(dp[-1])
def tournament(q,code):
 rows=[0]*q
 for k,(u,v) in enumerate(itertools.combinations(range(q),2)):
  if code>>k&1:rows[u]|=1<<v
  else:rows[v]|=1<<u
 return rows
def chain(order):
 rows=[0]*len(order)
 for i,u in enumerate(order):
  for v in order[i+1:]:rows[u]|=1<<v
 return rows
def extend(rows,a,c):
 n=len(rows);out=list(rows)+[1<<(n+1),0]
 for v in range(n):
  for k,m in enumerate((a,c)):
   if m>>v&1:out[v]|=1<<(n+k)
   else:out[n+k]|=1<<v
 return out
def substitute(q,module):
 r=len(module);out=list(module)+[0]*(len(q)-1)
 for u in range(1,len(q)):
  for v in range(1,len(q)):
   if q[u]>>v&1:out[r+u-1]|=1<<(r+v-1)
  for v in range(r):
   if q[u]&1:out[r+u-1]|=1<<v
   else:out[v]|=1<<(r+u-1)
 return out
def profile(rows):
 n=len(rows);h=[1]
 for s in range(1,1<<n):
  ids=[u for u in range(n) if s>>u&1]
  h.append(hp([sum(1<<j for j,v in enumerate(ids) if rows[u]>>v&1) for u in ids]))
 dp=[0]*(1<<n);dp[0]=1;out=[]
 for k in range(1,n+1):
  nxt=[0]*(1<<n)
  for s in range(1,1<<n):
   b=s
   while b:nxt[s]+=h[b]*dp[s^b];b=(b-1)&s
  dp=nxt;out.append(dp[-1])
 return out
def P(r,i,t):return t**(r+2)+t**(r-i+1)*((t+1)**i-(t-1)**i)
def polynomial(rows):
 n=len(rows);out=defaultdict(int)
 for word in itertools.permutations(range(n+2)):
  terms={0:1}
  for u,v in zip(word,word[1:]):
   if u<n and v<n:
    if not rows[u]>>v&1:break
   elif u>=n and v>=n:
    if u!=n:break
   else:
    bit=1<<((v-n)*n+u if u<n else (u-n)*n+v);nxt=defaultdict(int)
    for m,z in terms.items():
     if u<n:nxt[m|bit]+=z
     else:nxt[m]+=z;nxt[m|bit]-=z
    terms=nxt
  else:
   for m,z in terms.items():out[m]+=z
 return {m:z for m,z in out.items() if z}
def diff(a,b):return {m:a.get(m,0)-b.get(m,0) for m in a.keys()|b.keys() if a.get(m,0)!=b.get(m,0)}
def pairs(rows):return [list((u,v)) for u,v in itertools.combinations(range(len(rows)),2) if (rows[u]^rows[v])&~((1<<u)|(1<<v))==0]
def main():
 checks=0
 for r in range(1,7):
  last=None
  for i in range(r):
   p=profile(extend(chain(range(r)),0,1<<i))
   for t in range(r+4):need(sum(math.comb(t,k)*v for k,v in enumerate(p,1) if k<=t)==P(r,i,t),'color polynomial')
   if last is not None:need(p[0]>last[0] and all(a>=b for a,b in zip(p,last)),'coefficient monotonicity')
   last=p;checks+=1
 ambient=0
 for qn in range(1,5):
  for code in range(1<<(qn*(qn-1)//2)):
   q=tournament(qn,code)
   for r in range(1,5):
    vals=[hp(substitute(q,extend(chain(range(r)),0,1<<i))) for i in range(r)]
    need(all(x<y for x,y in zip(vals,vals[1:])),'arbitrary quotient strictness');ambient+=r
 q=[2,4,1];T=substitute(q,chain([0,1,2,3]));U=substitute(q,chain([1,3,0,2]));n=len(T)
 need(not set(map(tuple,pairs(T)))&set(map(tuple,pairs(U))),'no shared twin pair')
 O=[((1<<n)-1)^(1<<u)^T[u] for u in range(n)]
 ps=[polynomial(a) for a in (T,U,O)];du=diff(ps[0],ps[1]);do=diff(ps[0],ps[2])
 need(max(m.bit_count() for m in du)==4,'arbitrary full quartic');need(max(m.bit_count() for m in do)<=3,'global cubic')
 masks=[s for s in range(1<<n) if s.bit_count()<=2];queries=0
 for rows,p in zip((T,U,O),ps):
  for a in masks:
   for c in masks:
    m=a|(c<<n);need(hp(extend(rows,a,c))==sum(z for s,z in p.items() if s&m==s),'independent full basis query');queries+=1
 gate=[]
 seeds=[(16,5,33,39,46,3),(16,5,49,39,42,3),(18,4,49,39,42,3)]
 for j,s in enumerate(seeds):
  # Move root 5 to root 0; other names retain their relative order.
  ids=[5,0,1,2,3,4];q=[sum(1<<k for k,v in enumerate(ids) if s[u]>>v&1) for u in ids]
  for r in range(1,6):
   vals=[]
   for i in range(r):
    val=hp(substitute(q,extend(chain(range(r)),0,1<<i)))
    need(val==(2*j-3)*P(r,i,1)+9*P(r,i,2)+4*P(r,i,3),'gate exact law');vals.append(val)
   gate.append({'state':j,'r':r,'rank_queries':vals})
 result={'status':'PASS','profile_controls':checks,'ambient_queries':ambient,'full_basis_queries':queries,'T':T,'U':U,'pairs_T':pairs(T),'pairs_U':pairs(U),'arbitrary_difference':{str(m):v for m,v in sorted(du.items())},'global_converse_difference':{str(m):v for m,v in sorted(do.items())},'gate':gate}
 print(json.dumps(result,indent=2,sort_keys=True))
if __name__=='__main__':main()
