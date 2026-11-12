"""Session280 receiver rearming: literal paths, finite seeds, exact infinite-family identities."""
import json
from pathlib import Path
from itertools import combinations
from functools import lru_cache

def require(p,msg):
 if not p:raise ValueError(msg)

def endpoints(T):
 n=len(T);N=1<<n;dp=[[0]*n for _ in range(N)];H=[0]*N;H[0]=1
 for v in range(n):dp[1<<v][v]=1
 for s in range(1,N):
  for v in range(n):
   x=dp[s][v]
   if not x:continue
   H[s]+=x;todo=T[v]&~s
   while todo:
    bit=todo&-todo;todo-=bit;dp[s|bit][bit.bit_length()-1]+=x
 return H,dp

def profile(T):
 n=len(T);N=1<<n;full=N-1;R=full^1
 H,dp=endpoints(T);V=e=0
 for s in range(1,N,2):
  val=dp[s][0]*H[full^s];V+=val;e+=val*(1 if s.bit_count()%2 else -1)
 Q=k=0;s=R
 while True:
  val=H[s]*H[R^s];Q+=val;k+=val*(-1 if s.bit_count()%2 else 1)
  if not s:break
  s=(s-1)&R
 rev=tuple(sum(1<<j for j in range(n) if T[j]>>i&1) for i in range(n))
 Hr,dr=endpoints(rev);U=d=0
 for s in range(1,N,2):
  val=dr[s][0]*Hr[full^s];U+=val;d+=val*(1 if s.bit_count()%2 else -1)
 return dict(H=H[full],h=H[R],Q=Q,k=k,U=U,V=V,d=d,e=e)

def context(word):
 n=len(word);a=sum((b=='1')<<i for i,b in enumerate(word))
 return tuple([a<<1]+[((((1<<n)-1)^((1<<(i+1))-1))<<1)|(1-(a>>i&1)) for i in range(n)])

def fuse(A,B):
 n=len(A);m=len(B);r=[0]*(n+m-1)
 for T,ids in ((A,list(range(n))),(B,[0]+list(range(n,n+m-1)))):
  for i in range(len(T)):
   for j in range(len(T)):
    if T[i]>>j&1:r[ids[i]]|=1<<ids[j]
 for i in range(1,n):
  for j in range(n,n+m-1):r[i]|=1<<j
 return tuple(r)

def fm(a,b):
 return dict(h=a['h']*b['h'],Q=a['Q']*b['Q'],k=a['k']*b['k'],
 H=a['h']*b['H']+(a['H']-a['h'])*b['h']+(a['U']-a['h'])*(b['V']-b['h'])+(a['d']-a['h'])*(b['e']-b['h']),
 U=a['h']*b['U']+(a['U']-a['h'])*b['Q'],
 d=a['h']*b['d']+(a['d']-a['h'])*b['k'],
 V=a['Q']*b['V']+(a['V']-a['Q'])*b['h'],
 e=a['k']*b['e']+(a['e']-a['k'])*b['h'])

@lru_cache(None)
def g(p):
 if p==0:return 1
 v=p//2
 return g(v) if p%2 else g(v)+v+v%2

def ctx_profile(p,n=None):
 if n is None:n=p.bit_length()
 if n==0:return dict(H=1,h=1,Q=1,k=1,U=1,V=1,d=1,e=1)
 rev=sum((p>>i&1)<<(n-1-i) for i in range(n))
 return dict(H=g(p),h=1,Q=1<<n,k=0,U=p+1,V=(1<<n)-rev,d=1-(p&1),e=(p>>(n-1))&1)

def tournament(T):
 for i in range(len(T)):
  require(not(T[i]>>i&1),'loop')
  for j in range(i):require(((T[i]>>j)&1)+((T[j]>>i)&1)==1,'pair')

def is_frame(a,b,f):
 return sorted(f)==list(range(len(a))) and all((a[i]&b[i])>>j&1 for i,j in zip(f,f[1:]))

def main():
 T41=(56,13,25,48,34,6);T43=(52,13,24,49,34,6);T45=(52,13,24,33,42,6)
 R=(386,508,489,497,485,193,257,320,32)
 bank=[fuse(context('100'),t) for t in (T41,T43)]+[R]
 profiles=list(map(profile,bank));require([p['H'] for p in profiles]==[269,271,273],'base counts')
 require(all((p['h'],p['V'],p['e'])==(15,441,15) for p in profiles),'shared receiver')
 require(profiles[2]==dict(H=273,h=15,Q=768,k=0,U=441,V=441,d=9,e=15),'repair packet')
 upper=fuse(context('100'),T45)
 diff=lambda a,b:[(i,j) for i in range(len(a)) for j in range(i+1,len(a)) if (a[i]^b[i])>>j&1]
 require(diff(upper,R)==[(0,5),(2,4),(4,7),(4,8)],'four edits')
 require(diff(bank[1],R)==[(0,5),(2,4),(4,7),(4,8),(6,7)],'five free edits')
 require(is_frame(bank[1],R,[0,1,2,3,4,6,8,5,7]),'five-edit frame')
 require(is_frame(upper,R,[0,1,2,3,4,5,7,6,8]),'zero-edit frame')
 require(is_frame(bank[0],bank[1],[1,2,3,0,7,4,5,6,8]),'first bank edge frame')
 for T in bank+[upper]:tournament(T)
 # Literal endpoint-cover inequalities and old obstruction across nonunit contexts.
 tri=(2,4,1)
 source=(14,4,8,2);sink=(0,5,9,3);mixed=(2,4,9,3)
 contexts=[(0,),(2,0),(0,1),tri,source,sink,mixed]
 for t in contexts:tournament(t)
 pcontexts=list(map(profile,contexts));pt43=profile(T43);pt45=profile(T45)
 sandwich=[]
 for li,L in enumerate(contexts):
  lp=pcontexts[li]
  require(lp['h']<=lp['V']<=lp['Q'],'universal bounds')
  require(lp['V']==lp['h'] or lp['V']>=2*lp['h'],'receiver gap')
  require(lp['V']==lp['Q'] or lp['Q']-lp['V']>=lp['h'],'complement gap')
  for ri,Rc in enumerate(contexts):
   rp=pcontexts[ri];x=profile(fuse(fuse(L,T43),Rc));y=profile(fuse(fuse(L,T45),Rc))
   delta=2*(lp['Q']-lp['V'])*rp['h']-16*lp['Q']*(rp['V']-rp['h'])
   require(y['V']-x['V']==delta,'sandwich literal defect')
   expected=(lp['Q']==lp['V'] and rp['V']==rp['h'])
   require((delta==0)==expected,'sandwich iff')
   if expected:require(y['H']-x['H']==2*lp['h']*rp['h'],'multiplication cost')
   sandwich.append(delta)
 # Actual unit contexts; inherited formula checked against graph DP, then new bank transport.
 unit=0
 for n in range(5):
  for p in range(1<<n):
   w=f'{p:0{n}b}' if n else '';C=context(w)
   cp=ctx_profile(p,n);require(profile(C)==cp,'unit formula')
   out=[fm(cp,t) for t in profiles]
   require([r['H']-out[0]['H'] for r in out]==[0,2,4],'arbitrary unit translate')
   require(len({(r['h'],r['V'],r['e']) for r in out})==1,'arbitrary unit rearm')
   if n<=2:
    for T,r in zip(bank,out):require(profile(fuse(C,T))==r,'literal bank fusion')
   unit+=1
 # Literal two repetitions, all three graphs, and exact all-scale algebraic consequences.
 C1=context('10');pc=profile(C1);cur=bank
 for k in range(3):
  for i,T in enumerate(cur):
   p=profile(T);require((p['H'],p['h'],p['V'],p['e'])==(269+2*i+872*(4**k-1)//3+10*k,15,436*4**k+5,15),'literal repeated rearm')
  cur=[fuse(C1,T) for T in cur]
 levels=[]
 for k in range(0,41,2):
  V=436*4**k+5;z=(V-21)//15;require((V-21)%15==0 and z>0 and z%2==0,'carry integer')
  P=2*z;H=269+872*(4**k-1)//3+10*k
  shift=lambda p:15*(g(p)-1)+(V-15)*p
  require(shift(P+1)-shift(P)==6,'carry six')
  values=[H+shift(P)+2*i for i in range(3)]+[H+shift(P+1)+2*i for i in range(3)]
  require(values==list(range(values[0],values[0]+12,2)),'six interval')
  levels.append(dict(k=k,V=str(V),P=str(P),N=str(values[0])))
 literal=[];receivers=[]
 for P in (56,57):
  C=context(bin(P)[2:])
  for T in bank:
   pp=profile(fuse(C,T));literal.append(pp['H']);receivers.append(pp['V'])
 require(literal==[24875,24877,24879,24881,24883,24885],'literal six graph counts')
 require(receivers==[28119]*3+[27639]*3,'six union receiver leak')
 # Hostile old receiver leak is not silently removed.
 old=[fm(pc,profile(t)) for t in (T43,T45)]
 require(old[1]['H']-old[0]['H']==2 and old[1]['V']-old[0]['V']==2,'old first leak')
 again=[fm(pc,p) for p in old];require(again[1]['H']-again[0]['H']==6,'old next count failure')
 general=0
 for h in range(1,12):
  for r in range(1,8):
   for z in range(2,22,2):
    V=h+2*r+h*z;P=2*z
    shift=lambda p:h*(g(p)-1)+(V-h)*p
    require(shift(P+1)-shift(P)==2*r,'general width doubling')
    cp=ctx_profile(P);cq=ctx_profile(P+1);n=P.bit_length()
    require(cp['Q']==cq['Q'] and cq['V']-cp['V']==-(1<<(n-1)),'general receiving cost')
    general+=1
 result=dict(status='PASS',base_profiles=profiles,sandwiches_literal=len(sandwich),unit_context_controls=unit,literal_six_counts=literal,literal_six_receivers=receivers,general_algebra_controls=general,even_levels=levels,scope='Physical rearmed three-bank and infinite six-wide intervals; no all-odd coverage or renewable six-way bank')
 Path(__file__).with_name('RESULTS.json').write_text(json.dumps(result,indent=2)+'\n')
 print(json.dumps(result,indent=2))
if __name__=='__main__':main()
