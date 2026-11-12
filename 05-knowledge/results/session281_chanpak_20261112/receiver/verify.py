"""Session281 receiver: self-contained literal seed and exact formula controls.
Recovery edition: definitions formerly imported from Session280 are copied here,
so reproduction does not depend on a missing checkout. Mathematical scope unchanged.
"""
from pathlib import Path
import json
from functools import lru_cache
P=Path(__file__).resolve().parent
def check(p,msg):
 if not p:raise RuntimeError(msg)
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
 U=a['h']*b['U']+(a['U']-a['h'])*b['Q'],d=a['h']*b['d']+(a['d']-a['h'])*b['k'],
 V=a['Q']*b['V']+(a['V']-a['Q'])*b['h'],e=a['k']*b['e']+(a['e']-a['k'])*b['h'])
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
  check(not(T[i]>>i&1),'loop')
  for j in range(i):check(((T[i]>>j)&1)+((T[j]>>i)&1)==1,'pair')
A=(6,12,24,17,3);B=(2,4,1);C=(2,4,9,3);D=(36,1,2,7,15,30);CoreE=(44,1,2,6,15,30)
def ordinal(parts):
 n=sum(map(len,parts));out=[];offset=0
 for part in parts:
  tail=((1<<n)-1)^((1<<(offset+len(part)))-1)
  out.extend((a<<offset)|tail for a in part);offset+=len(part)
 return tuple(out)
def tt(n):return [(0,)]*n
def root(D,I):return tuple([(((1<<len(D))-1)^I)<<1]+[(v<<1)|((I>>i)&1) for i,v in enumerate(D)])
def endpoint_weights(D):
 H,dp=endpoints(D);n=len(D);full=(1<<n)-1;E=[0]*n;F=[0]*n
 for s in range(1,full+1):
  sign=-1 if s.bit_count()%2 else 1
  for i in range(n):
   x=dp[s][i]*H[full^s];E[i]+=x;F[i]+=sign*x
 Q=sum(H[s]*H[full^s] for s in range(full+1))
 K=sum((-1 if s.bit_count()%2 else 1)*H[s]*H[full^s] for s in range(full+1))
 return H[full],Q,K,E,F
specs=[(tt(3)+[CoreE]+tt(1),236),([B]+tt(1)+[C]+tt(2),408),(tt(4)+[CoreE],76),([CoreE]+tt(4),326),(tt(3)+[C,B],204),(tt(1)+[C,B]+tt(2),368)]
raw=[root(ordinal(parts),I) for parts,I in specs];raw_profiles=[profile(T) for T in raw]
check([x['H'] for x in raw_profiles]==list(range(979,991,2)),'raw six counts')
check(all((x['h'],x['V'])==(15,1203) for x in raw_profiles),'raw receiver')
check([x['e'] for x in raw_profiles]==[15,15,15,13,15,15],'raw e mismatch')
bank=[fuse(context('10'),T) for T in raw];profiles=[profile(T) for T in bank]
check([x['H'] for x in profiles]==list(range(3385,3397,2)),'six counts')
check(all((x['h'],x['V'],x['e'])==(15,4797,15) for x in profiles),'six receiver')
for T in raw+bank:tournament(T)
controls=0
for parts in [(B,C),(C,B),((0,),A),(A,(0,)),(B,(0,),C),((0,),B,(0,),C)]:
 data=[endpoint_weights(D) for D in parts];out=endpoint_weights(ordinal(parts));EE=[];FF=[]
 for j,d in enumerate(data):
  leftQ=leftK=rightH=1
  for a in data[:j]:leftQ*=a[1];leftK*=a[2]
  for a in data[j+1:]:rightH*=a[0]
  EE.extend(leftQ*x*rightH for x in d[3]);FF.extend(leftK*x*rightH for x in d[4])
 check(EE==out[3] and FF==out[4],'ordinal weights');controls+=1
for (parts,I),p in zip(specs,raw_profiles):
 h,Q,K,E,F=endpoint_weights(ordinal(parts));v=h+sum(x for i,x in enumerate(E) if I>>i&1);e=h+sum(x for i,x in enumerate(F) if I>>i&1)
 check((h,v,e)==(p['h'],p['V'],p['e']),'deletion receiver')
contexts=0
for n in range(5):
 for p in range(1<<n):
  cp=ctx_profile(p,n);result=[fm(cp,x) for x in profiles]
  check(profile(context(f'{p:0{n}b}' if n else ''))==cp,'context profile')
  check([x['H']-result[0]['H'] for x in result]==list(range(0,12,2)),'renewable spacing')
  check(len({(x['h'],x['V'],x['e']) for x in result})==1,'renewable receiver')
  if n==1:
   for T,pr in zip(bank,result):check(profile(fuse(context(f'{p:01b}'),T))==pr,'literal common context')
  contexts+=1
levels=[]
for k in range(0,41,2):
 V=4792*4**k+5;H=3385+9584*(4**k-1)//3+10*k
 check(V%30==27,'carry congruence');z=(V-27)//15;P0=2*z
 shift=lambda p:15*(g(p)-1)+(V-15)*p
 vals=[H+shift(P0)+2*i for i in range(6)]+[H+shift(P0+1)+2*i for i in range(6)]
 check(vals==list(range(vals[0],vals[0]+24,2)),'twelve consecutive')
 n=P0.bit_length();check((P0+1).bit_length()==n,'equal context order')
 cp=ctx_profile(P0);cq=ctx_profile(P0+1)
 check(cq['V']-cp['V']==-(1<<(n-1)),'paid leak')
 levels.append(dict(k=k,V=str(V),P=str(P0),N=str(vals[0]),order=13+2*k+n,receiver_defect=str(-15*(1<<(n-1)))))
check(levels[0]['N']=='3051967','first interval')
preparation_controls=0
for h in range(1,32,2):
 for V in range(h,h+21):
  for r in range(1,21):
   orbit={};x=V%h;j=0
   while x not in orbit:orbit[x]=j;x=2*x%h;j+=1
   if r%h not in orbit:continue
   n=orbit[r%h]+1;period=j-orbit[x]
   while n<2 or (1<<n)*V-h<=h+2*r:n+=period
   VV=(1<<n)*V-h
   check((VV-h-2*r)%(2*h)==0,'prepared phase')
   cp=ctx_profile(1<<(n-1),n)
   check((cp['H'],cp['U'],cp['V'],cp['d'],cp['e'])==(1+(1<<(n-1)),1+(1<<(n-1)),(1<<n)-1,1,1),'preparation context')
   z=(VV-h-2*r)//h;P0=2*z
   shift=lambda p:h*(g(p)-1)+(VV-h)*p
   check(z>0 and z%2==0 and shift(P0+1)-shift(P0)==2*r,'prepared carry');preparation_controls+=1
check({pow(2,j,7) for j in range(3)}=={1,2,4},'gcd hostile')
result=dict(status='PASS',raw_rows=raw,raw_profiles=raw_profiles,seed_rows=bank,seed_profiles=profiles,core_weights={str(D):endpoint_weights(D) for D in [(0,),A,B,C,D,CoreE]},ordinal_controls=controls,unit_context_controls=contexts,preparation_algebra_controls=preparation_controls,levels=levels,scope='Renewable six-bank and infinite twelve-count families; widened union not rearmed; Chan-Pak open')
P.joinpath('RESULTS.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
