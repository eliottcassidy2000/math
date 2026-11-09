import itertools as it
B=[40,1,3,6,15,30]
def states(B):
 n=len(B);N=1<<n;h=[1]+[0]*(N-1);ep=[[0]*n for _ in range(N)];sp=[[0]*n for _ in range(N)]
 for s in range(1,N):
  for v in range(n):
   if s>>v&1:
    r=s^(1<<v)
    ep[s][v]=sum(ep[r][u] for u in range(n) if B[u]>>v&1) if r else 1
    sp[s][v]=sum(sp[r][u] for u in range(n) if B[v]>>u&1) if r else 1
  h[s]=sum(ep[s])
 ans=[]
 for m in range(N):
  es=[sum(ep[s][v] for v in range(n) if m>>v&1) for s in range(N)]
  ss=[sum(sp[s][v] for v in range(n) if m>>v&1) for s in range(N)]
  ans.append((h[-1],sum(h[s]*h[N-1-s] for s in range(N)),sum(ss[s]*h[N-1-s] for s in range(N)),sum(es[s]*h[N-1-s] for s in range(N)),sum(es[s]*ss[N-1-s] for s in range(N)),ss[-1],es[-1]))
 return ans
if __name__=='__main__':
 st=states(B)
 for m,k in [(17,40),(24,38),(19,42),(21,44),(25,39),(23,46)]:print(m,k,st[m],st[k],tuple(b-a for a,b in zip(st[m],st[k])))
def compose(x,y):
 h,q,A,G,T,S,E=x;j,k,B,J,U,V,W=y
 return (h*j,q*k,A*k+h*B,J*q+j*G,A*J+T*j+h*U+E*V,S*j,h*W)
if __name__=='__main__':
 cs=states([16,1,3,7,14]);x=compose(cs[2],cs[28]);y=compose(cs[4],cs[0]);print('carry',x,y,tuple(b-a for a,b in zip(x,y)))
