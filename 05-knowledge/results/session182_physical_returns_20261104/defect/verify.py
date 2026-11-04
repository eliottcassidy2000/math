"""Exact named-state transport and clone-spectrum checks; standard library only."""
import itertools,json,math,pathlib
ROOT=pathlib.Path(__file__).resolve().parent

def need(x,msg):
 if not x: raise RuntimeError(msg)
def add(*ps):
 d={}
 for p in ps:
  for k,v in p.items():d[k]=d.get(k,0)+v
 return {k:v for k,v in d.items() if v}
def mul(*ps):
 d={0:1}
 for p in ps:
  q={}
  for a,x in d.items():
   for b,y in p.items():
    if not a&b:q[a|b]=q.get(a|b,0)+x*y
  d=q
 return d
def matrix(rows):
 n=len(rows);g=[[{} for _ in rows] for _ in rows]
 for k in range(1,n+1):
  for p in itertools.permutations(range(n),k):
   if all(rows[i]>>j&1 for i,j in zip(p,p[1:])):
    mask=sum(1<<i for i in p);a=g[p[0]][p[-1]];a[mask]=a.get(mask,0)+1
 return g
def birth(rows,cut):
 n=len(rows)
 return [r|(((cut>>i)&1)<<n) for i,r in enumerate(rows)]+[((1<<n)-1)^cut]
def strong(rows):
 for root in range(len(rows)):
  seen=1<<root
  while True:
   new=seen
   for v,r in enumerate(rows):
    if seen>>v&1:new|=r
   if new==seen:break
   seen=new
  if seen!=(1<<len(rows))-1:return False
 return True
def direct(rows):
 n=len(rows);h=0;f=[0]*n;j=[[0]*n for _ in rows]
 for p in itertools.permutations(range(n)):
  edges=[bool(rows[a]>>b&1) for a,b in zip(p,p[1:])]
  if all(edges):h+=1
  for k in range(1,n+1):
   if all(edges[:k-1]) and all(edges[k:]):
    f[p[k-1]]+=1
    if k<n:j[p[k-1]][p[k]]+=1
 return h,f,j
def tournaments(n):
 es=list(itertools.combinations(range(n),2))
 for code in range(1<<len(es)):
  a=[0]*n
  for k,(i,j) in enumerate(es):a[i if code>>k&1 else j]|=1<<(j if code>>k&1 else i)
  yield a
def h_dp(rows):
 n=len(rows);dp=[[0]*n for _ in range(1<<n)]
 for s in range(1,1<<n):
  for v in range(n):
   if s>>v&1:
    q=s^(1<<v);dp[s][v]=1 if not q else sum(dp[q][u] for u in range(n) if rows[u]>>v&1)
 return sum(dp[-1])
def clone(rows,x,t):
 old=[v for v in range(len(rows)) if v!=x];n=len(old);a=[0]*(n+t)
 for i,u in enumerate(old):
  for j,v in enumerate(old):
   if rows[u]>>v&1:a[i]|=1<<j
  for k in range(t):
   if rows[u]>>x&1:a[i]|=1<<(n+k)
   else:a[n+k]|=1<<i
 for k in range(t):
  for l in range(k+1,t):a[n+k]|=1<<(n+l)
 return a
def coefficients(rows,x):
 old=[v for v in range(len(rows)) if v!=x];u=[sum(1<<j for j,v in enumerate(old) if rows[w]>>v&1) for w in old];g=matrix(u);n=len(u)
 I=[i for i,v in enumerate(old) if rows[v]>>x&1];O=[i for i,v in enumerate(old) if rows[x]>>v&1]
 L=add({0:1},*(g[i][j] for i in range(n) for j in I));R=add({0:1},*(g[i][j] for i in O for j in range(n)));M=add(*(g[i][j] for i in O for j in I));p=mul(L,R);A=[]
 for r in range(1,n//2+2):A.append(p.get((1<<n)-1,0));p=mul(p,M)
 need(not p,'nilpotence bound')
 return A
def surj(t,r):return sum((-1)**(r-k)*math.comb(r,k)*k**t for k in range(r+1))
def main():
 cases=matrices=guards=clones=0;one={0:1}
 for n in range(1,5):
  full=(1<<n)-1
  for rows in tournaments(n):
   g=matrix(rows);e=[add(*(g[i][j] for i in range(n))) for j in range(n)];s=[add(*r) for r in g];P=add(*e)
   for cut in range(1<<n):
    born=birth(rows,cut);gg=matrix(born);h,f,j=direct(born)
    X=add(*(e[i] for i in range(n) if cut>>i&1));Y=add(*(s[i] for i in range(n) if not cut>>i&1));A=add(one,X);B=add(one,Y)
    need(mul(A,B).get(full,0)==h,'h birth')
    if strong(rows):need(strong(born)==(cut not in (0,full)),'strong guard');guards+=1
    # Full endpoint matrix, including both new physical endpoints.
    for i in range(n+1):
     for k in range(n+1):
      if i==n and k==n:pred={1<<n:1}
      elif i==n:
       p=add(*(g[v][k] for v in range(n) if not cut>>v&1));pred={m|(1<<n):v for m,v in p.items()}
      elif k==n:
       p=add(*(g[i][v] for v in range(n) if cut>>v&1));pred={m|(1<<n):v for m,v in p.items()}
      else:
       p=mul(add(*(g[i][v] for v in range(n) if cut>>v&1)),add(*(g[v][k] for v in range(n) if not cut>>v&1)))
       pred=add(g[i][k],{m|(1<<n):v for m,v in p.items()})
      need(pred==gg[i][k],'named endpoint birth')
    for w in range(n+1):
     for v in range(n+1):
      if w==n and v==n:pred=0
      elif w==n:pred=mul(A,s[v]).get(full,0)
      elif v==n:pred=mul(e[w],B).get(full,0)
      else:
       Cw=add(*(g[k][w] for k in range(n) if not cut>>k&1));Dv=add(*(g[v][k] for k in range(n) if cut>>k&1))
       pred=add(mul(A,Cw,s[v]),mul(e[w],Dv,B)).get(full,0)
      need(pred==j[w][v],'individual joint transport')
     need(f[w]-h==sum(j[w][v] for v in range(n+1) if born[v]>>w&1),'inner defect actual incoming ports')
    matrices+=1
    for wm in range(1<<n):
     W=[w for w in range(n) if wm>>w&1];E=add(*(e[w] for w in W));S=add(*(s[w] for w in W))
     C=add(*(g[v][w] for v in range(n) if not cut>>v&1 for w in W));D=add(*(g[w][v] for w in W for v in range(n) if cut>>v&1))
     for delta in (0,1):
      CP=add(C,one) if delta else C;DP=add(D,one) if delta else D;WW=W+([n] if delta else [])
      F=add(mul(A,CP,add(one,P)),mul(E,A,B)).get(full,0)
      J=add(mul(A,CP,S),mul(E,DP,B)).get(full,0)
      need(F==sum(f[w] for w in WW),'multiport F');need(J==sum(j[w][v] for w in WW for v in WW),'joint J')
      need(h-J==h-sum(j[w][v] for w in WW for v in WW),'alpha')
      cases+=1
 for n in range(1,5):
  for rows in tournaments(n):
   for x in range(n):
    A=coefficients(rows,x)
    for r,a in enumerate(A,1):need(a%math.factorial(r-1)==0,'ordered middle divisibility')
    for t in range(1,5):
     out=clone(rows,x,t);need(h_dp(out)==sum(a*surj(t,r) for r,a in enumerate(A,1)),'clone spectrum')
     if strong(rows) and n>1:need(strong(out),'clone strongness')
     clones+=1
 # Nontrivial ternary spectrum: demonstrates that the obstruction is not merely absence of base 3.
 sample=[20,1,2,7,14];AA=coefficients(sample,0);need(AA==[9,10,2],'three-run sample')
 for t in range(1,8):need(h_dp(clone(sample,0,t))==2*3**t+4*2**t-5,'mixed spectrum')
 pure=[26,24,67,116,36,71,19];PA=coefficients(pure,6)
 need(strong(pure) and PA==[145,150,50,0],'pure ternary coefficients')
 for t in range(1,7):need(h_dp(clone(pure,6,t))==50*3**t-5,'pure ternary +10 response')
 out={'status':'PASS','pure_ternary_control':{'rows':pure,'physical_clone_vertex':6,'A':PA,'H_t':'50*3^t-5','digit':10,'tested_t':[1,6]},'joint_birth_cases':cases,'full_matrix_birth_cases':matrices,'local_strong_guards':guards,'clone_response_checks':clones,'mixed_spectrum':{'rows':sample,'physical_clone_vertex':0,'A':AA,'H_t':'2*3^t+4*2^t-5','tested_t':[1,7]},'scope':'All labelled tournaments order1..4, all birth cuts, all old W and newborn port memberships; all cloned vertices t1..4. No ternary closure asserted.'}
 (ROOT/'verification.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
if __name__=='__main__':main()
