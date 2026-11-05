import json
from functools import lru_cache
@lru_cache(None)
def stats(out):
 n=len(out); size=1<<n; full=size-1
 end=[[0]*n for _ in range(size)]
 h=[0]*size;h[0]=1
 for S in range(1,size):
  for v in range(n):
   if S>>v&1:
    U=S^(1<<v)
    end[S][v]=1 if U==0 else sum(end[U][u] for u in range(n) if U>>u&1 and out[u]>>v&1)
  h[S]=sum(end[S])
 return h,end
@lru_cache(None)
def g(v):
 if v==0:return 1
 u=v//2;return g(u)+(u+(u&1) if v%2==0 else 0)
def apex(v):
 bits=[int(x) for x in bin(v)[2:]] if v else []
 n=len(bits);out=[0]*(n+1)
 for i,b in enumerate(bits):
  for j in range(i+1,n):out[i]|=1<<j
  if b:out[n]|=1<<i
  else:out[i]|=1<<n
 return out,n

def join(v,S,root):
 A,a=apex(v);n=len(A);m=len(S);out=A+[0]*m
 for i in range(n):out[i]|=((1<<m)-1)<<n
 out[a]^=1<<(n+root);out[n+root]|=1<<a
 for i in range(m):
  for j in range(m):
   if S[i]>>j&1:out[n+j]|=1<<(n+i)
 return tuple(out)
if __name__ == '__main__':
 rows=[list(map(int,s.split())) for s in open('source_bank7.graphs.txt')]
 w=[]
 for row in rows:
  n,r,h,e,R,*out=row
  for v in range(127//R):
   if h*g(v)-(v%2==0)*e+(v+1)*R==127:
    J=join(v,out,r);Hs,_=stats(J)
    if Hs[-1]!=127: raise RuntimeError('join H != 127')
    masks=[s for s,x in enumerate(Hs) if x==41]
    w.append(dict(order=len(J),bank_order=n,bank_root=r,bank_triple=[h,e,R],v=v,bank_out=out,join_out=J,H41_masks=masks))
 print(json.dumps({'F2_H127_witnesses':len(w),'max_order':max(x['order'] for x in w),'H41_containing':sum(bool(x['H41_masks']) for x in w),'witnesses':w},indent=2))
