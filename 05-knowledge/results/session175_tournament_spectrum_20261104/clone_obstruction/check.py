#!/usr/bin/env python3
import itertools,json

def require(x,s):
 if not x:raise RuntimeError(s)
def labelled(n,code):
 out=[0]*n;k=0
 for i in range(n):
  for j in range(i+1,n):
   if code>>k&1:out[i]|=1<<j
   else:out[j]|=1<<i
   k+=1
 return out

def h(out):
 n=len(out);dp=[[0]*n for _ in range(1<<n)]
 for v in range(n):dp[1<<v][v]=1
 for s in range(1,1<<n):
  for v,val in enumerate(dp[s]):
   if not val:continue
   t=out[v]&~s
   while t:
    bit=t&-t;t-=bit;dp[s|bit][bit.bit_length()-1]+=val
 return sum(dp[-1])
def induce(g,vs):return [sum(1<<j for j,w in enumerate(vs) if g[v]>>w&1) for v in vs]
def strong(g):
 n=len(g)
 for v in range(n):
  seen=1<<v
  while True:
   nex=seen
   for w in range(n):
    if seen>>w&1:nex|=g[w]
   if nex==seen:break
   seen=nex
  if seen!=(1<<n)-1:return False
 return True

def inflate(g,rs):
 blocks=[];offset=0
 for r in rs:blocks.append(list(range(offset,offset+r)));offset+=r
 out=[0]*offset
 for i,vs in enumerate(blocks):
  for a,v in enumerate(vs):
   out[v]=sum(1<<w for w in vs[a+1:])
   for j,ws in enumerate(blocks):
    if g[i]>>j&1:out[v]|=sum(1<<w for w in ws)
 return out

def main():
 insertion=triangles=blowups=0
 for n in range(1,6):
  for code in range(1<<(n*(n-1)//2)):
   g=labelled(n,code);total=h(g)
   if n>1:
    for v in range(n):require(total>=h(induce(g,[w for w in range(n) if w!=v])),'monotonicity');insertion+=1
   if n>=3 and strong(g):
    for v in range(n):
     require(any(g[v]>>a&1 and g[a]>>b&1 and g[b]>>v&1 for a in range(n) for b in range(n)),'triangle');triangles+=1
    if n<=4:
     for rs in itertools.product(range(1,4),repeat=n):
      G=inflate(g,rs);H=h(G)
      require(strong(G),'strong inflation')
      require(all(H>=2**r+1 for r in rs),'parameter bound');blowups+=1
 for r in range(1,10):require(h(inflate([2,4,1],[r,1,1]))==2**r+1,'triangle formula')
 print(json.dumps({'status':'PASS','vertex_deletion_monotonicity_checks':insertion,'strong_vertex_triangle_checks':triangles,'strong_blowup_parameter_checks':blowups,'inflated_triangle_sizes':list(range(1,10))},indent=2))
if __name__=='__main__':main()
