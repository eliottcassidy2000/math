"""Independent exact controls for cyclic-run endpoint separation. Stdlib only."""
from itertools import combinations,product
from functools import lru_cache
from pathlib import Path
import hashlib,json,math


def need(ok,msg):
 if not ok:raise RuntimeError(msg)
def rows(n,code):
 T=[0]*n
 for k,(u,v) in enumerate(combinations(range(n),2)):
  if code>>k&1:T[u]|=1<<v
  else:T[v]|=1<<u
 return T
def extend(T,a,c):
 n=len(T);R=list(T)+[1<<(n+1),0]
 for u in range(n):
  if a>>u&1:R[u]|=1<<n
  else:R[n]|=1<<u
  if c>>u&1:R[u]|=1<<(n+1)
  else:R[n+1]|=1<<u
 return R
def starts(T):
 n=len(T);N=(1<<n)-1
 @lru_cache(None)
 def rec(S,u):
  if S==N:return 1
  return sum(rec(S|1<<v,v) for v in range(n) if not S>>v&1 and T[u]>>v&1)
 return [rec(1<<u,u) for u in range(n)]
def hp(T):return sum(starts(T))
def induced(T,vs):return [sum(1<<j for j,v in enumerate(vs) if T[u]>>v&1) for u in vs]
def cycle(blocks,op=False):
 offsets=[0,len(blocks[0]),len(blocks[0])+len(blocks[1])];N=sum(map(len,blocks));T=[0]*N
 for i,A in enumerate(blocks):
  for u,r in enumerate(A):T[offsets[i]+u]|=r<<offsets[i]
  j=(i+(-1 if op else 1))%3
  for u in range(len(A)):T[offsets[i]+u]|=((1<<len(blocks[j]))-1)<<offsets[j]
 return T,offsets
def profile(T):
 n=len(T);N=(1<<n)-1;H=[1]+[hp(induced(T,[u for u in range(n) if S>>u&1])) for S in range(1,N+1)]
 D=[[0]*(n+1) for _ in range(N+1)];D[0][0]=1
 for S in range(1,N+1):
  bit=S&-S;part=S
  while part:
   if part&bit:
    for k in range(1,n+1):D[S][k]+=H[part]*D[S^part][k-1]
   part=(part-1)&S
 return [D[N][k]*math.factorial(k) for k in range(n+1)]
def stirling(m):
 s=[1]
 for n in range(1,m+1):s=[0]+[k*(s[k] if k<len(s) else 0)+s[k-1] for k in range(1,n+1)]
 return s
def skew(ps,i):
 n=max(map(len,ps));p=[x+[0]*(n-len(x)) for x in ps];j=(i+1)%3;l=(i+2)%3
 return sum(p[i][k]*(p[j][k]*p[l][k-1]-p[j][k-1]*p[l][k]) for k in range(1,n))

def main():
 tomography=0
 for n in range(1,5):
  for code in range(1<<(n*(n-1)//2)):
   T=rows(n,code);ss=starts(T);H=sum(ss)
   for u in range(n):
    f=hp(extend(T,1<<u,0));g=hp(extend(T,0,1<<u))
    need(f-2*g+H==2*ss[u],'start tomography');tomography+=1
 profiles=[stirling(m) for m in range(1,31)];minors=0
 for b in range(2,31):
  for c in range(1,b):
   P=profiles[b-1];R=profiles[c-1]+[0]*(b-c)
   for k in range(1,b+1):
    need(P[k]*R[k-1]>=P[k-1]*R[k],'Stirling ratio minor');minors+=1
 # Arbitrary small children validate the exact run formula independently of Stirling positivity.
 children=[[0],[2,0],[0,1],[6,4,0],[2,4,1],[4,1,2]]
 run_cases=0
 for blocks in product(children,repeat=3):
  if sum(map(len,blocks))>8:continue
  T,offsets=cycle(blocks);U,_=cycle(blocks,True);S=starts(T);R=starts(U);ps=list(map(profile,blocks))
  for i,A in enumerate(blocks):
   need(sum(S[offsets[i]:offsets[i]+len(A)])-sum(R[offsets[i]:offsets[i]+len(A)])==skew(ps,i),'cyclic run identity');run_cases+=1
 examples=[]
 for sizes in product(range(2,5),repeat=3):
  if sum(sizes)>10 or len(set(sizes))==1:continue
  blocks=[rows(m,(1<<(m*(m-1)//2))-1) for m in sizes];T,offsets=cycle(blocks);U,_=cycle(blocks,True)
  ps=[[math.factorial(k)*v for k,v in enumerate(stirling(m))] for m in sizes]
  i=next(i for i in range(3) if sizes[(i+1)%3]!=sizes[(i+2)%3]);want=skew(ps,i)
  need(want!=0,'strict nonzero skew');need((want>0)==(sizes[(i+1)%3]>sizes[(i+2)%3]),'skew sign')
  data=[]
  for u in range(offsets[i],offsets[i]+sizes[i]):
   for a,c in ((1<<u,0),(0,1<<u)):data.append([a,c,hp(extend(T,a,c)),hp(extend(U,a,c))])
  need(sum((1 if c==0 else -2)*(f-g) for a,c,f,g in data)==2*want,'literal signed battery')
  need(any(f!=g for a,c,f,g in data),'no actual separator')
  examples.append({'sizes':sizes,'block':i,'start_skew':want,'queries':data})
 bad=[8,1,3,6];p=profile(bad);r=profile(induced(bad,[1,2,3]));need(p==[0,5,22,36,24] and r==[0,1,6,6],'arbitrary-child hostile');need(p[2]*r[1]-p[1]*r[2]==-8,'hostile minor')
 root=Path(__file__).resolve().parent
 print(json.dumps({'status':'PASS','tomography_checks':tomography,'stirling_minors_through_order_30':minors,'arbitrary_child_run_checks':run_cases,'literal_unequal_batteries':examples,'arbitrary_child_hostile':{'rows':bad,'deleted':0,'profile':p,'deleted_profile':r,'minor':-8},'source_sha256':{name:hashlib.sha256((root/name).read_bytes()).hexdigest() for name in ('PROOF.md','verify.py')}},indent=2,sort_keys=True))
if __name__=='__main__':main()
