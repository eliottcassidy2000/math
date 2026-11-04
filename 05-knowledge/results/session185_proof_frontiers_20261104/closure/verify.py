"""Literal strong-suffix escape certificate; Python standard library only."""
import json,itertools
from pathlib import Path

def need(x,msg):
 if not x:raise RuntimeError(msg)
def tables(r):
 n=len(r);N=1<<n;h=[0]*N;h[0]=1;e=[[0]*n for _ in range(N)]
 for v in range(n):e[1<<v][v]=1
 for S in range(1,N):
  for v in range(n):
   if S>>v&1:
    T=S^(1<<v)
    if T:e[S][v]=sum(e[T][u] for u in range(n) if r[u]>>v&1)
  h[S]=sum(e[S])
 return h,e
def direct(r,W):
 n=len(r);h,e=tables(r);rev=[sum(1<<j for j in range(n) if r[j]>>i&1) for i in range(n)];_,s=tables(rev);N=(1<<n)-1
 ep=[sum(x[i] for i in W) for x in e];sp=[sum(x[i] for i in W) for x in s]
 return(h[N],sum(h[x]*h[N^x] for x in range(N+1)),sum(ep[x]*h[N^x] for x in range(N+1)),sum(h[x]*sp[N^x] for x in range(N+1)),sum(ep[x]*sp[N^x] for x in range(N+1)),ep[N],sp[N])
def strong(r):
 N=(1<<len(r))-1
 for v in range(len(r)):
  S=1<<v
  while True:
   T=S
   for i in range(len(r)):
    if S>>i&1:T|=r[i]
   if T==S:break
   S=T
  if S!=N:return False
 return True
def tournament(r):
 n=len(r)
 return all(not(r[i]>>i&1) and r[i]<(1<<n) for i in range(n)) and all(bool(r[i]>>j&1)!=bool(r[j]>>i&1) for i in range(n) for j in range(i))
def attach(r,W,t):
 n=len(r);root=t+n
 rows=[sum(1<<j for j in range(i+1,t+n)) for i in range(t)]
 rows += [(x<<t)|((1<<root) if i in W else 0) for i,x in enumerate(r)]
 rows += [((1<<t)-1)|sum(1<<(t+i) for i in range(n) if i not in W)]
 return rows

def run():
 ws=json.loads(Path(__file__).with_name('witnesses.json').read_text())
 # Original entry and siblings are literal ordered suffixes.
 states={'':{'rows':[6,4,0],'mask':7,'sig':[1,8,7,7,6,1,1]},'2':{'rows':[30,28,8,16,4],'mask':7,'sig':[3,48,21,39,16,3,1]},'4':{'rows':[30,28,24,16,0],'mask':21,'sig':[1,32,21,21,12,1,1]}}
 # Compute rather than trust nonessential source signature coordinates.
 for x in states.values():x['sig']=list(direct(x['rows'],{i for i in range(len(x['rows'])) if x['mask']>>i&1}))
 for x in ws:states[x['address']]=x
 report=[];shifts={'':0}
 for name,x in states.items():
  r=x['rows'];W={i for i in range(len(r)) if x['mask']>>i&1};z=tuple(x['sig'])
  need(tournament(r),'tournament validity');need(direct(r,W)==z,'literal seven signature coordinates')
  if name not in ('','2','4'):need(strong(r),'suffix strong')
  for t in (1,2):
   a=attach(r,W,t);need(tournament(a) and strong(a),'full attachment strong')
   need(tables(a)[0][-1]==z[0]-z[4]+2**t*z[2],'literal attachment Hamilton count')
  if name:
   p=states[name[:-1]]['sig'];d=int(name[-1]);q=0;lhs=z[2];rhs=3*p[2]
   while lhs<rhs:lhs*=2;q+=1
   while lhs>rhs:rhs*=2;q-=1
   need(lhs==rhs,'dyadic slope relation');need(z[0]-z[4]==3*(p[0]-p[4])+d,'digit constant relation')
   shifts[name]=shifts[name[:-1]]+q
   report.append({'source':name[:-1],'target':name,'digit':d,'local_clock_shift':q,'total_clock_shift':shifts[name]})
 # A simultaneous concrete witness at root clock4.
 counts={name:x['sig'][0]-x['sig'][4]+2**(4+shifts[name])*x['sig'][2] for name,x in states.items()}
 for name in states:
  if name:need(counts[name]==3*counts[name[:-1]]+int(name[-1]),'common-clock digit counts')
 # Unbounded normalized valuation with constants coprime to3 infinitely often.
 seed=(11,100,64,49,28,9,6)
 need(direct((24,1,3,6,14),{1,2,3})==seed,'valuation-family seed')
 def mul(x,y):
  h,p,f,g,j,e,s=x;H,P,F,G,J,E,S=y
  return(h*H,p*P,f*H+p*F,g*P+h*G,j*H+h*J+e*S+g*F,h*E,s*H)
 A=671088;block=None
 for i in range(20):
  a=(A>>i)&1;v=(1,2,a,a,0,a,a);block=v if block is None else mul(block,v)
 need(block==(1,1048576,671088,60357,38626,1,0),'twenty singleton block')
 z=seed
 for k in range(1,31):
  z=mul(z,block)
  need(z[2]==2**(6+20*k),'unbounded normalized valuation')
  need((z[0]-z[4])%3==(1+k)%3,'eligible constant residue')
 result={'status':'PASS','literal_states':len(states),'new_strong_suffixes':len(ws),'attachment_clocks_checked':[1,2],'edges':report,'root_clock4_counts':counts,'valuation_family':'F=2^(6+20k), c=1+k mod3; k=1..30 checked; all k>=1 proved'}
 Path(__file__).with_name('receipt.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,sort_keys=True))
if __name__=='__main__':run()
