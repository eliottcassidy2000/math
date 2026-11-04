import itertools,json,argparse
from collections import defaultdict

def need(b,s):
 if not b:raise RuntimeError(s)
def ham(r):
 n=len(r);dp=[[0]*n for _ in range(1<<n)]
 for u in range(n):dp[1<<u][u]=1
 for mask in range(1,1<<n):
  for u in range(n):
   val=dp[mask][u];bits=r[u]&~mask
   while bits:
    b=bits&-bits;bits-=b;dp[mask|b][b.bit_length()-1]+=val
 return sum(dp[-1])
def query(r,mask):
 n=len(r); t=list(r)+[1<<(n+1),0]
 for port in range(2):
  for u in range(n):
   if mask>>(port*n+u)&1:t[u]|=1<<(n+port)
   else:t[n+port]|=1<<u
 return ham(t)
def polynomial(r):
 n=len(r);result=defaultdict(int)
 for word in itertools.permutations(range(n+2)):
  terms={0:1}
  for u,v in zip(word,word[1:]):
   if u<n and v<n:
    if not r[u]>>v&1:break
   elif u>=n and v>=n:
    if u!=n:break
   else:
    bit=1<<((v-n)*n+u if u<n else (u-n)*n+v)
    nxt=defaultdict(int)
    for s,c in terms.items():
     if u<n:nxt[s|bit]+=c
     else:nxt[s]+=c;nxt[s|bit]-=c
    terms=nxt
  else:
   for s,c in terms.items():result[s]+=c
 return {s:c for s,c in result.items() if c}
def expand(r):
 q=len(r);t=list(r)+[r[0]];t[0]|=1<<q
 for v in range(1,q):
  if r[v]&1:t[v]|=1<<q
 return t
def strong(r):
 n=len(r)
 for u in range(n):
  s=1<<u
  while True:
   t=s
   for v in range(n):
    if s>>v&1:t|=r[v]
   if s==t:break
   s=t
  if s!=(1<<n)-1:return False
 return True
def contact(r,mask):
 n=len(r);q=n-1;u=0;v=q;d=((mask>>u)&1)-((mask>>v)&1);e=((mask>>(n+u))&1)-((mask>>(n+v))&1)
 a=(mask>>u)&1;b=(mask>>v)&1;c=(mask>>(n+u))&1;f=(mask>>(n+v))&1
 quotient=[r[i]&((1<<q)-1) for i in range(q)]
 counts=[0,0,0,0]
 for word in itertools.permutations(range(q+2)):
  good=True
  for x,y in zip(word,word[1:]):
   if (x==0 and y>=q) or (y==0 and x>=q):continue
   if x<q and y<q:ok=quotient[x]>>y&1
   elif x>=q and y>=q:ok=x==q
   elif x<q:ok=mask>>((y-q)*n+x)&1
   else:ok=1-(mask>>((x-q)*n+y)&1)
   if not ok:good=False;break
  if not good:continue
  k=word.index(0);left=word[k-1] if k else -1;right=word[k+1] if k+1<len(word) else -1
  x= q in (left,right);y=q+1 in (left,right)
  if x and not y:counts[0]+=1
  elif y and not x:counts[1]+=1
  elif x and y:counts[2 if left==q else 3]+=1
 return -d*counts[0]-e*counts[1]+((1-a)*f-(1-b)*c)*counts[2]+((1-c)*b-(1-f)*a)*counts[3]

def cover_skew(r,removed,z,j):
 vertices=[v for v in range(len(r)) if v!=removed]
 result=0
 for p in itertools.permutations(vertices):
  for k in range(1,len(p)):
   if (p[k-1],p[k]) not in ((z,j),(j,z)):continue
   if all(r[p[t]]>>p[t+1]&1 for t in range(len(p)-1) if t!=k-1):
    result+=1 if (p[k-1],p[k])==(z,j) else -1
 return result

def main():
 counts={}
 for q in range(1,5):
  hist=defaultdict(int)
  for code in range(1<<(q*(q-1)//2)):
   r=[0]*q
   for k,(u,v) in enumerate(itertools.combinations(range(q),2)):
    if code>>k&1:r[u]|=1<<v
    else:r[v]|=1<<u
   t=expand(r);tt=t[:];tt[0]^=1<<q;tt[q]^=1
   p=polynomial(t);pp=polynomial(tt)
   d={s:p.get(s,0)-pp.get(s,0) for s in p.keys()|pp.keys() if p.get(s,0)!=pp.get(s,0)}
   for i,j in itertools.permutations(range(1,q),2):
    mon=(1<<i)|(1<<j)|(1<<(q+1+i))|(1<<(q+1+q))
    need(d.get(mon,0)==cover_skew(r,i,0,j),'deletion skew coefficient')
   degree=max(s.bit_count() for s in d)
   need(degree<=4,'degree')
   if q<=3:need(degree<=3,'small quartic')
   hist[str(degree)]+=1
  counts[str(q)]={'quotients':1<<(q*(q-1)//2),'degrees':dict(hist)}
 T=[24,17,19,6,8];N=[8,17,19,6,9];need(strong(T) and strong(N),'strong hostile')
 mon=(1<<1)|(1<<3)|(1<<6)|(1<<9)
 controls=[];coefficient=0
 sub=mon
 while True:
  a=query(T,sub);b=query(N,sub);sgn=(-1)**(4-sub.bit_count());coefficient+=sgn*(a-b)
  controls.append({'mask':sub,'F_T':a,'F_N':b,'sign':sgn})
  if sub==0:break
  sub=(sub-1)&mon
 need(coefficient==2,'finite difference')
 p=polynomial(T);pp=polynomial(N);need(p.get(mon,0)-pp.get(mon,0)==2,'word coefficient')
 # Independent full truth tables, then Boolean interpolation.
 values=[query(T,s)-query(N,s) for s in range(1<<10)]
 need(all(contact(T,s)==values[s] for s in range(1<<10)), 'contact factorization')
 mob=values[:]
 for k in range(10):
  for s in range(1<<10):
   if s>>k&1:mob[s]-=mob[s^(1<<k)]
 need(all(mob[s]==p.get(s,0)-pp.get(s,0) for s in range(1<<10)),'full coefficients')
 print(json.dumps({'status':'PASS','quotient_census':counts,'strong_hostile':{'T':T,'N':N,'monomial':'a_1 a_3 c_1 c_4','coefficient':coefficient,'queries':controls},'full_truth_queries':2048,'degree_counts':{str(k):sum(v!=0 and s.bit_count()==k for s,v in enumerate(mob)) for k in range(5)}},indent=2,sort_keys=True))
if __name__=='__main__':main()
