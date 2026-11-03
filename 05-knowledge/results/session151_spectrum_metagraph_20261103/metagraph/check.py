from itertools import permutations
from collections import Counter,defaultdict
from math import factorial

def run(n):
 p=[(i,j) for i in range(n) for j in range(i+2,n)]; m=len(p); pairs=[(i,j) for i in range(n) for j in range(i+1,n)]; idx={x:k for k,x in enumerate(pairs)}
 perms=list(permutations(range(n))); pmaps=[[idx[tuple(sorted((q[i],q[j])))] for i,j in pairs] for q in perms]
 def full(t):
  return sum((1<<idx[e]) for k,e in enumerate(p) if t>>k&1) # bit 1 reverse edge; path naturally forward
 def permute(b,q):
  return sum(((((b>>idx[tuple(sorted((q[i],q[j])))])&1) ^ (q[i]>q[j]))<<idx[i,j]) for i,j in pairs)
 def canon(b):return min(permute(b,q) for q in perms)
 C=[canon(full(t)) for t in range(1<<m)]; conv={c:canon(c^((1<<len(pairs))-1)) for c in set(C)}
 merged=[min(c,conv[c]) for c in C]; sig=[p.index((n-1-j,n-1-i)) for i,j in p]
 blue=lambda t:all((t>>i&1)==(t>>j&1) for i,j in enumerate(sig))
 def mirror(t):return sum(((t>>i)&1)<<j for i,j in enumerate(sig))
 if not (sum(i==j for i,j in enumerate(sig))==(n-1)//2): raise AssertionError("metagraph check failed")
 if not (all(C[mirror(t)]==conv[C[t]] for t in range(1<<m))): raise AssertionError("metagraph check failed")
 bundles=Counter()
 for t in range(1<<(m-1)):
  v=t^((1<<m)-1)
  if not blue(t): bundles[tuple(sorted((merged[t],merged[v])))]+=1
 if not (all(k%2==0 for k in bundles.values())): raise AssertionError("metagraph check failed")
 if not (all(C.count(c)%2==1 for c in set(C))): raise AssertionError("metagraph check failed")
 if not (all(sum(blue(t) for t in range(1<<m) if C[t]==c)%2==int(conv[c]==c) for c in set(C))): raise AssertionError("metagraph check failed")
 if n==5:
  neighbors=set()
  for t in range(1<<(m-1)):
   v=t^((1<<m)-1); a,b=merged[t],merged[v]
   if blue(t) and (a==76 or b==76): neighbors.add(b if a==76 else a)
  if not (len(neighbors)==3 and 76 not in neighbors): raise AssertionError("metagraph check failed")
  print('n=5 pure-blue code 76 has three distinct blue neighbors:',sorted(neighbors))
 print('CHECKS PASS: S converse, fixed coordinates, color fiber parity, all black bundles including loops even')
 def stats(ids):
  nodes=sorted(set(ids)); s={u:Counter() for u in nodes}
  for t in range(1<<(m-1)):
   v=t^((1<<m)-1); a,b=ids[t],ids[v]; col='B' if blue(t) else 'K';
   s[a][col+'d']+=1; s[b][col+'d']+=1
   if a==b:s[a][col+'loop']+=1
   else:s[a][col+'cross']+=1;s[b][col+'cross']+=1
  result=[]
  for u in nodes:
   c=next(C[t] for t,x in enumerate(ids) if x==u); a=s[u];mass=ids.count(u); aut=sum(permute(c,q)==c for q in perms)
   result.append(dict(code=u,SC=conv[c]==c,mass=mass,aut=aut,blue=a['Bd'],black=a['Kd'],blue_loop=a['Bloop'],black_loop=a['Kloop'],blue_cross=a['Bcross'],black_cross=a['Kcross']))
  return result
 print('\nn=',n,'m=',m,'classes=',len(set(C)),'merged=',len(set(merged)))
 for label,ids in [('original',C),('merged',merged)]:
  print(label)
  for x in stats(ids):print(x)
for n in (3,4,5):run(n)
