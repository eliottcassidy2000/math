"""Exact apex-cut menus. Literal subset DP, not fusion matrices."""
from itertools import combinations,product
from collections import defaultdict
import json
from pathlib import Path

def endpoints(T):
 n=len(T);E=[[0]*n for _ in range(1<<n)]
 for i in range(n):E[1<<i][i]=1
 for mask in range(1,1<<n):
  for i in range(n):
   if not E[mask][i]:continue
   for j in range(n):
    if not(mask>>j&1) and T[i]>>j&1:E[mask|1<<j][j]+=E[mask][i]
 return E

def data(T):
 n=len(T);full=(1<<n)-1;E=endpoints(T);O=tuple(sum(1<<j for j in range(n) if T[j]>>i&1) for i in range(n));S=endpoints(O)
 H=[sum(z) for z in E];H[0]=1;h=H[-1]
 W=[[0]*n for _ in range(n)]
 for mask in range(1,full):
  for i in range(n):
   for j in range(n):W[j][i]+=E[mask][j]*S[full^mask][i]
 a=[sum(S[mask][i]*H[full^mask] for mask in range(1,full+1)) for i in range(n)]
 b=[sum(E[mask][i]*H[full^mask] for mask in range(1,full+1)) for i in range(n)]
 c=[sum((-1)**mask.bit_count()*S[mask][i]*H[full^mask] for mask in range(1,full+1)) for i in range(n)]
 d=[sum((-1)**mask.bit_count()*E[mask][i]*H[full^mask] for mask in range(1,full+1)) for i in range(n)]
 Q=sum(H[mask]*H[full^mask] for mask in range(full+1))
 return dict(h=h,Q=Q,W=W,a=a,b=b,c=c,d=d,start=S[-1],end=E[-1])

def cut(D,mask):
 n=len(D['W']);x=[mask>>i&1 for i in range(n)]
 return sum(x[i]*D['start'][i]+(1-x[i])*D['end'][i] for i in range(n))+sum((1-x[j])*x[i]*D['W'][j][i] for i in range(n) for j in range(n))

def apex(T,mask):return (mask<<1,)+tuple((T[i]<<1)|(1-(mask>>i&1)) for i in range(len(T)))
def circulant(n,code):
 ds=[i+1 if code>>i&1 else n-i-1 for i in range((n-1)//2)]
 return tuple(sum(1<<((v+d)%n) for d in ds) for v in range(n))

def run():
 reports=[]
 for n in [3,5,7,9,11]:
  for code in range(1<<((n-1)//2)):
   T=circulant(n,code);D=data(T);menus=defaultdict(dict)
   for mask in range(1<<n):menus[mask.bit_count()].setdefault(cut(D,mask),mask)
   entries=[]
   for r,menu in menus.items():
    vals=sorted(menu);runs=[]
    for v in vals:
     if runs and v==runs[-1][-1]+2:runs[-1].append(v)
     else:runs.append([v])
    best=max(runs,key=len)
    entries.append(dict(r=r,values=vals,longest_run=best,witnesses={v:menu[v] for v in best}))
   reports.append(dict(n=n,code=code,T=T,h=D['h'],Q=D['Q'],a=D['a'],c=D['c'],menus=entries))
   print(n,code,D['h'],max((len(z['longest_run']),z['longest_run'][0],z['longest_run'][-1],z['r']) for z in entries),flush=True)
 Path(__file__).with_name('central_menus.json').write_text(json.dumps(reports,indent=2)+'\n')
if __name__=='__main__':run()
