"""Independent permutation oracle for the two-hole uniform padding bound."""
import itertools,json

def H(g):
 return sum(all(g[a]>>b&1 for a,b in zip(p,p[1:])) for p in itertools.permutations(range(len(g))))
def graph(n,code):
 g=[0]*n
 for k,(a,b) in enumerate(itertools.combinations(range(n),2)):
  if code>>k&1:g[a]|=1<<b
  else:g[b]|=1<<a
 return g

def join_tail(g,m):
 n=len(g);return [r+(((1<<m)-1)<<n) for r in g]+[sum(1<<j for j in range(n+i+1,n+m)) for i in range(m)]
def blow(q,parts):
 offsets=[];n=0
 for p in parts:offsets.append(n);n+=len(p)
 out=[0]*n
 for i,p in enumerate(parts):
  for v,r in enumerate(p):
   out[offsets[i]+v]=r<<offsets[i]
   for j,t in enumerate(parts):
    if q[i]>>j&1:out[offsets[i]+v]|=((1<<len(t))-1)<<offsets[j]
 return out

def unsafe(q,i):
 return any(q[i]>>a&1 and q[a]>>b&1 and q[b]>>i&1 for a in range(len(q)) for b in range(len(q)))
inputs=([0],[2,0],[2,4,1]);checks=0
for code in range(8):
 q=graph(3,code)
 for a,b in itertools.product(inputs,repeat=2):
  hp=H(a)*H(b)
  for m in (0,1):
   if len(a)+len(b)+2*m+1>8:continue
   out=H(blow(q,[join_tail(a,m),join_tail(b,m),[0]]))
   if any(unsafe(q,i) for i in (0,1)):
    if out<(H(q)+2**(m+1)-2)*hp:raise RuntimeError(('lower',code,a,b,m))
   elif out!=H(q)*hp:raise RuntimeError(('safe',code,a,b,m))
   checks+=1
print(json.dumps({'status':'PASS','independent_permutation_controls':checks,'universe':'Every labeled order3 quotient; holes 0,1; singleton/TT2/C3 inputs; padding m=0,1; total order <=8'},indent=2))
