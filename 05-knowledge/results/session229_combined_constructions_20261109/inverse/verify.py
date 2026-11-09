import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parent))
from inherited import *
def tournaments(n):
 for mask in range(1<<(n*(n-1)//2)):
  rows=[0]*n;k=0
  for i in range(n):
   for j in range(i+1,n):
    u,v=(i,j) if mask>>k&1 else (j,i);rows[u]|=1<<v;k+=1
  yield rows
def strong(a):
 for i in range(len(a)):
  seen=1<<i
  while True:
   nxt=seen
   for j in range(len(a)):
    if seen>>j&1:nxt|=a[j]
   if nxt==seen:break
   seen=nxt
  if seen!=(1<<len(a))-1:return False
 return True
def join(a,b):
 n=len(a);return [r|(((1<<len(b))-1)<<n) for r in a]+[r<<n for r in b]
def main():
 atoms={};starts={};gadgets={}
 for n in range(1,7):
  for a in tournaments(n):
   hp,_,_,_=counts(a,0);H=hp[-1];Q=sum(hp[s]*hp[(1<<n)-1-s] for s in range(1<<n))
   if n<=4:gadgets.setdefault((H,Q),a)
   if not strong(a):continue
   if Q<=31:atoms.setdefault((H,Q),a)
   if H>20:continue
   for r in range(n):
    x=tuple(profile(a,r))
    if x[2]<=31:starts.setdefault(x,(a,r))
 states=dict(starts);todo=list(starts)
 while todo:
  x=todo.pop();a,r=states[x];D,E,U,d=x
  for (h,q),b in atoms.items():
   y=(D*h,0,U*q,d*h)
   if y[0]>20 or y[2]>31:continue
   if y not in states:states[y]=(join(a,b),r);todo.append(y)
 for x,(a,r) in states.items():
  require(tuple(profile(a,r))==x,"literal source vs SCC transfer")
 out={}
 for x,(a,r) in states.items():
  if not a[r]:continue
  D,E,U,d=x
  for (H,Q),b in gadgets.items():
   N=H*(D-E+U)+(Q-2*H)*(U-d)
   if N>31:continue
   g,rr=graft(a,r,b);require(profile(g,rr)[0]==N,'literal graft')
   require(strong(g),'strong output')
   require(H+Q<=N and H<N and D<N and U<=N//H-1,'all bounds')
   require(N>=H*(3*D+1)//2+Q-2*H,'descent')
   out.setdefault(N,{'source':a,'root':r,'state':x,'gadget':b,'H_Q':[H,Q]})
 require(set(out)==set(range(3,32,2))-{7,21},'exact image')
 require(profile(graft([0],0,[4,1,2])[0],1)[0]==3,'sink descent hostile')
 print(json.dumps({'status':'FINITE-EXACT PASS','scope':'Complete reachable-root proper-graft image through31, using cited strong floor f(7)=25','source_D_cap':20,'source_U_cap':31,'labelled_initial_graphs':33867,'states':len(states),'image_through31':sorted(out),'witnesses':out},indent=2))
if __name__=='__main__':main()
