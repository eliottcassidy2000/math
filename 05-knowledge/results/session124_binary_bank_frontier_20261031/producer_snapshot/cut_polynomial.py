from full_binary import h,extend

def cut_counts(out):
 n=len(out);N=1<<n;rev=tuple(sum(1<<j for j in range(n) if out[j]>>i&1) for i in range(n))
 def ends(edges):
  d=[[0]*n for _ in range(N)]
  for i in range(n):d[1<<i][i]=1
  for m in range(1,N):
   for i in range(n):
    if d[m][i]:
     todo=edges[i]&~m
     while todo:
      b=todo&-todo;j=b.bit_length()-1;d[m|b][j]+=d[m][i];todo-=b
  return d
 E=ends(out);R=ends(rev);H=sum(E[-1]);W=[[0]*n for _ in range(n)]
 for m in range(1,N-1):
  q=N-1-m
  for u in range(n):
   if E[m][u]:
    for v in range(n):W[u][v]+=E[m][u]*R[q][v]
 L=[R[-1][v]-E[-1][v]+sum(W[u][v] for u in range(n)) for v in range(n)]
 F=[H]*N
 for c in range(1,N):
  bit=c&-c;v=bit.bit_length()-1;b=c-bit
  F[c]=F[b]+L[v]-sum(W[u][v]+W[v][u] for u in range(n) if b>>u&1)
 return F

if __name__=='__main__':
 import random,json
 rng=random.Random(124);checks=0
 for n in range(1,8):
  for _ in range(20):
   out=[0]*n
   for i in range(n):
    for j in range(i+1,n):
     if rng.randrange(2):out[i]|=1<<j
     else:out[j]|=1<<i
   out=tuple(out);F=cut_counts(out)
   for cut,f in enumerate(F):
    actual=h(extend(out,cut))
    if actual!=f:raise ValueError((out,cut,f,actual))
    checks+=1
 print(json.dumps({'independent_cut_dp_checks':checks}))
