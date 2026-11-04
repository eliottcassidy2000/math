import itertools,json,hashlib
from pathlib import Path

def need(b,m):
 if not b:raise RuntimeError(m)
def signature(rows,W):
 n=len(rows);N=(1<<n)-1
 need(all(not rows[i]>>i&1 for i in range(n)),'loop')
 need(all(((rows[i]>>j)&1)+((rows[j]>>i)&1)==1 for i in range(n) for j in range(i+1,n)),'tournament')
 E=[[0]*n for _ in range(N+1)];S=[[0]*n for _ in range(N+1)];h=[0]*(N+1);h[0]=1
 for A in range(1,N+1):
  for v in range(n):
   if A>>v&1:
    B=A^(1<<v)
    if not B:E[A][v]=S[A][v]=1
    else:
     E[A][v]=sum(E[B][u] for u in range(n) if rows[u]>>v&1)
     S[A][v]=sum(S[B][u] for u in range(n) if rows[v]>>u&1)
  h[A]=sum(E[A]);need(h[A]==sum(S[A]),'startend')
 e=[sum(row[v] for v in range(n) if W>>v&1) for row in E]
 s=[sum(row[v] for v in range(n) if W>>v&1) for row in S]
 return(h[N],sum(h[A]*h[N^A] for A in range(N+1)),sum(e[A]*h[N^A] for A in range(N+1)),sum(h[A]*s[N^A] for A in range(N+1)),sum(e[A]*s[N^A] for A in range(N+1)),e[N],s[N])
def join(a,b):
 if a is None:return b
 if b is None:return a
 h,p,f,g,j,e,s=a;H,P,F,G,J,E,S=b
 return(h*H,p*P,f*H+p*F,g*P+h*G,j*H+h*J+e*S+g*F,h*E,s*H)
def graphjoin(parts):
 rows=[];W=0
 for r,w in parts:
  n=len(rows);k=len(r)
  rows=[x|(((1<<k)-1)<<n) for x in rows]+[x<<n for x in r];W|=w<<n
 return rows,W
cases=[((2,20,9,3,45,79,31),42,(35,676,252,267,70,13,14)),((2,4,33,23,71,155,47,95),52,(89,2144,504,627,118,17,24))]
contexts=[((),0,None)]
for n in range(1,4):
 edges=list(itertools.combinations(range(n),2))
 for bits in range(1<<len(edges)):
  r=[0]*n
  for i,(u,v) in enumerate(edges):
   if bits>>i&1:r[u]|=1<<v
   else:r[v]|=1<<u
  for W in range(1<<n):contexts.append((r,W,signature(r,W)))
count=0;direct=0
for r,W,z in cases:
 need(signature(r,W)==z,'signature')
 c=z[0]-z[4]
 for x,y in itertools.product(contexts,repeat=2):
  q=join(join(x[2],z),y[2]);C=q[0]-q[4]
  hx=x[2][0] if x[2] else 1;hy=y[2][0] if y[2] else 1
  need(C<=hx*hy*c,'general inequality')
  need(C not in (3*c+2,3*c+4),'barrier');count+=1
 # Direct whole-graph checks against assembly, including both empty contexts.
 for x,y in [(contexts[0],contexts[-1]),(contexts[-1],contexts[0]),(contexts[-1],contexts[-1]),(contexts[1],contexts[2])]:
  rr,ww=graphjoin([(x[0],x[1]),(r,W),(y[0],y[1])]);q=signature(rr,ww)
  need(q==join(join(x[2],z),y[2]),'literal contextual signature');direct+=1
result={'status':'PASS','independent_full_signatures':2,'context_pair_inequalities_and_exclusions':count,'literal_whole_context_graphs':direct,'normal_empty_context_identity':'Explicit None identity, no fictitious endpoint coordinates.'}
print(json.dumps(result,sort_keys=True))
