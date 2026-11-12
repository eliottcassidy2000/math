from collections import deque,Counter
from itertools import combinations
import json,time

def coredata():
 u=[6,12,24,17,3];h=[0]*32;st=[[0]*5 for _ in h];en=[[0]*5 for _ in h];h[0]=1
 for mask in range(1,32):
  for i in range(5):
   if mask>>i&1:
    rest=mask^(1<<i)
    st[mask][i]=1 if not rest else sum(st[rest][j] for j in range(5) if u[i]>>j&1)
    en[mask][i]=1 if not rest else sum(en[rest][j] for j in range(5) if u[j]>>i&1)
  h[mask]=sum(st[mask])
 return h,st,en
CD=coredata()
def endpoints(mask,l,r):
 h,st,en=CD;L=mask&((1<<l)-1);C=(mask>>l)&31;R=mask>>(l+5)
 if not mask:return [],[]
 ss=([(L&-L).bit_length()-1,h[C]] if L else None)
 starts=[ss] if ss else [[l+i,v] for i,v in enumerate(st[C]) if v] if C else [[l+5+(R&-R).bit_length()-1,1]]
 ee=([l+5+R.bit_length()-1,h[C]] if R else None)
 ends=[ee] if ee else [[l+i,v] for i,v in enumerate(en[C]) if v] if C else [[L.bit_length()-1,1]]
 return starts,ends

def form(l,r):
 n=l+5+r;full=(1<<n)-1;C=[[0]*n for _ in range(n)]
 for mask in range(1,full):
  _,en=endpoints(mask,l,r);st,_=endpoints(full^mask,l,r)
  for i,a in en:
   for j,b in st:C[i][j]+=a*b
 st,en=endpoints(full,l,r);sd=dict(st);ed=dict(en)
 lin=[sd.get(j,0)-ed.get(j,0)+sum(C[i][j] for i in range(n)) for j in range(n)]
 pair=[[C[i][j]+C[j][i] for j in range(n)] for i in range(n)]
 hs=[0]*(1<<n);hs[0]=15
 for mask in range(1,1<<n):
  bit=mask&-mask;i=bit.bit_length()-1;rest=mask^bit
  v=hs[rest]+lin[i];rr=rest
  while rr:b=rr&-rr;rr^=b;v-=pair[i][b.bit_length()-1]
  hs[mask]=v
 return hs

def run(l,r):
 n=l+5+r;hs=form(l,r);initial=(28<<l)|(((1<<r)-1)<<(l+5)); assert hs[initial]==41
 # Relaxed root-star two-flip graph: shared-frame condition tested separately on witnesses.
 ds=[1<<i for i in range(n)]+[(1<<i)|(1<<j) for i in range(n) for j in range(i+1,n)]
 seen={initial};q=deque([initial]);parent={initial:None};best=initial
 while q:
  x=q.popleft()
  if hs[x]>hs[best]:best=x
  for d in ds:
   y=x^d
   if y not in seen and abs(hs[y]-hs[x])<=2:seen.add(y);q.append(y);parent[y]=x
 vals=sorted({hs[x] for x in seen});assert vals==list(range(vals[0],vals[-1]+1,2))
 path=[];p=best
 while p is not None:path.append(p);p=parent[p]
 return dict(left=l,right=r,order=n+1,total_states=len(hs),component_states=len(seen),count_range=[min(vals),max(vals)],path=path[::-1],path_counts=[hs[x] for x in path[::-1]])
if __name__=='__main__':
 import sys
 out=[]
 for total in range(int(sys.argv[1]) if len(sys.argv)>1 else 9):
  for l in range(total+1):
   z=run(l,total-l);out.append(z); print(json.dumps(z),flush=True)
