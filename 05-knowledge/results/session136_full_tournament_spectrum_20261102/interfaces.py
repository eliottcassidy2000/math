"""Independent small interfaces between Session136's mathematical lenses."""
import json

def h(a):
 n=len(a);dp=[[0]*n for _ in range(1<<n)]
 for i in range(n):dp[1<<i][i]=1
 for m in range(1,1<<n):
  for v,c in enumerate(dp[m]):
   if not c:continue
   for w in range(n):
    if a[v]>>w&1 and not m>>w&1:dp[m|1<<w][w]+=c
 return sum(dp[-1]) if n else 1

def word_graph(b):
 n=len(b);a=[sum(1<<j for j in range(i+1,n)) for i in range(n)]+[0]
 for i,x in enumerate(b):
  if x=='1':a[n]|=1<<i
  else:a[i]|=1<<n
 return a

def triangle_context(n):
 # Transitive n-vertex module A, A->b->c->A.
 a=[sum(1<<j for j in range(i+1,n))|(1<<n) for i in range(n)]
 return a+[1<<(n+1),(1<<n)-1]

def main():
 examples=[]
 for target,b,expected in [(7,'101010000010',2007),(13,'101011100010',2013),(21,'10001001000',1021)]:
  value=h(word_graph(b))
  if value!=expected or value%1000!=target or value==target:raise RuntimeError('local/exact boundary')
  examples.append({'target':target,'word':b,'order':len(b)+1,'actual_H':value,'modulus':1000})
 contexts=[]
 for n in range(1,9):
  value=h(triangle_context(n))
  if value!=2**n+1:raise RuntimeError('context bound equality')
  contexts.append({'module_order':n,'module_H':1,'context_H':value,'cap127':min(value,128)})
 return {'status':'PASS','scope':'Small literal subset-DP interfaces, not the universal proofs. The missing exact values and one-apex13 are inherited.','congruence_hostiles':examples,'contextual_padding_hostiles':contexts}
if __name__=='__main__':print(json.dumps(main(),indent=2,sort_keys=True))
