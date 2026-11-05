import itertools,json
from pathlib import Path
HERE=Path(__file__).resolve().parent
def need(ok,msg):
 if not ok:raise RuntimeError(msg)
def op(a,r,c):
 n=len(a);a=a[:]+[0]
 if c=='D':
  for i in range(n):a[i]|=1<<n
 else:
  for i in range(n):
   if i==r:a[n]|=1<<i
   else:a[i]|=1<<n
 if c=='A':r=n
 return a,r

def run(w,a=None,r=0):
 if a is None:a=[0]
 for c in w:a,r=op(a,r,c)
 return a,r

def counts(a,r):
 n=len(a);N=1<<n;dp=[[0]*n for _ in range(N)];dq=[[0]*n for _ in range(N)]
 for i in range(n):dp[1<<i][i]=1
 dq[1<<r][r]=1
 for S in range(1,N):
  for i in range(n):
   if not(S>>i&1):continue
   for j in range(n):
    if not(S>>j&1) and a[i]>>j&1:
     dp[S|1<<j][j]+=dp[S][i];dq[S|1<<j][j]+=dq[S][i]
 return sum(dp[-1]),sum(dq[-1])

def converse(a):return [sum(1<<j for j in range(len(a)) if j!=i and not(a[i]>>j&1)) for i in range(len(a))]
def induced(a,V):return [sum(1<<j for j,y in enumerate(V) if a[x]>>y&1) for x in V]
def endpoint(a,u,w):
 n=len(a);m=len(w);v=n+m;out=a[:]+[0]*(m+1)
 for i in range(n):
  if i==u:out[i]|=1<<v
  else:out[v]|=1<<i
 for j,b in enumerate(w):
  z=n+j;out[z]=(1<<n)-1
  for k in range(j+1,m):out[z]|=1<<(n+k)
  if b=='1':out[v]|=1<<z
  else:out[z]|=1<<v
 return out,v

def main():
 wordchecks=0;initialK=0
 for k in range(8):
  for letters in itertools.product('AKD',repeat=k):
   w=''.join(letters);a,r=run(w);t=len(w)-len(w.lstrip('K'));b,s=run(w[t:]);V=[0]+list(range(t+1,len(a)))
   need(induced(a,V)==b and V[s]==r,('strip_graph',w))
   need(all(a[i]>>j&1 for i in range(1,t+1) for j in V),('permanent_source',w))
   need(all(a[i]>>j&1 for i in range(1,t+1) for j in range(i+1,t+1)),('prefix_order',w))
   h,start=counts(a,r);h2,start2=counts(b,s)
   need(h==h2 and start2==1,('counts',w,h,h2,start2))
   need(start==(0 if t else 1),('initial_exception',w,start))
   wordchecks+=1;initialK+=bool(t)
 cores=[([0],0),([2,4,1],0),([2,12,9,17,7],2)]
 macrochecks=0
 for a,u in cores:
  n=len(a)
  for m in range(7):
   for bits in itertools.product('01',repeat=m):
    w=''.join(bits);target,v=endpoint(a,u,w)
    c='A'+''.join('K' if b=='1' else 'D' for b in w[::-1]);got,r=run(c,converse(a),u)
    order=list(range(n))+list(range(n+m,n,-1))+[n]
    need(induced(got,order)==converse(target) and order[v]==r,('macro_graph',a,u,w,c))
    macrochecks+=1
 receipt={'status':'PASS','all_elementary_words_through_length':7,'elementary_word_count':wordchecks,'words_with_nonempty_initial_K_prefix':initialK,'graph_macro_compilations':macrochecks,'tested_macro_lengths':'0..6, all bit words; includes nonstrong boundary cases','scope':'Graph-level compiler and source stripping support the proved all-order scalar-support equality; no unrestricted spectrum claim.'}
 (HERE/'receipt.json').write_text(json.dumps(receipt,sort_keys=True,indent=2)+'\n');print(json.dumps(receipt,sort_keys=True))
if __name__=='__main__':main()
