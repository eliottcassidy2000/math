import json,random,sys,pathlib,hashlib
from selector import F,count,bracket,scale
checks=0
def need(ok):
 global checks
 checks+=1
 if not ok:raise AssertionError('check '+str(checks))
def rev(x,t):return int(format(x,'0%db'%t)[::-1],2)
for x in range(256):
 need(1<=F(x)<=x+1)
 need(F(2*x+1)==F(x));need(F(2*x)==F(x)+x+(x&1))
 for t in range(1,6):
  L=1<<t
  for r in range(L):
   need(F(L*x+r)==F(x)+(L-1-rev(r,t))*x+(not(r>>(t-1)))*(x&1)+F(r)-1)
for t in range(1,8):
 L=1<<t
 for q in range(1,12):
  m=L*q-1
  for a in range(m):
   ds=abs(F(a+1)+F(m-a-1)-F(a)-F(m-a))
   need(L*ds<=L*abs(m-2*a)+2*L*L+m+5*L)
# Independent direct prefix-sum read, rather than the streaming F recurrence.
def oldF(x):
 bits=bin(x)[2:]
 return 1+sum((int(bits[:i],2) if i else 0)+(int(bits[i-1]) if i else 0)
              for i,b in enumerate(bits) if b=='0')
for a in range(1,30):
 for b in range(1,30):
  if (a+b)%2!=1:continue
  u,w=(a,b) if a%2 else (b,a)
  lo=91+5*(oldF(u)+oldF(w)-2)+142*(u+w)+256*u*w+2*(w%2)
  hi=91+5*(oldF(w)+oldF(u)-2)+142*(u+w)+256*u*w+2*(u%2)
  need(count(a,b)==lo and hi==lo+2)
rng=random.Random(27520261117);cases=[]
for bits in [19,20,23,24,30,40,60,100,200,500,1000,2000]:
 for j in range(15):
  N=rng.randrange(1<<(bits-1),1<<bits)
  z=bracket(N);lo,hi=[p['lower'] for p in z['pairs']]
  need(lo<=N<=hi);need((hi-lo)**3<=1000**3*N)
  need(z['evaluations']<=N.bit_length()+5)
  need(all(p['upper']==p['lower']+2 and p['u']%2==1 and p['w']%2==0 for p in z['pairs']))
  if j==0:cases.append({'bits':bits,'M':z['M'],'L':z['L'],'error':z['error'],'bracket_width':z['bracket_width'],'evaluations':z['evaluations']})
sectors=[]
for t in range(4,25):
 M=(1<<(2*t))-1;L=1<<t
 for N in [64*(M-1)**2+1,64*M*M-1,64*M*M]:
  z=bracket(N,L);need(z['delta']==0);need(z['bracket_width']**4<=10000**4*N)
 sectors.append({'t':t,'error':z['error'],'width':z['bracket_width']})
print(json.dumps({'status':'PASS','checks':checks,'universal_target_controls':cases,'quarter_sectors':sectors,
 'source_sha256':{p:hashlib.sha256((pathlib.Path(__file__).parent/p).read_bytes()).hexdigest() for p in ['SELECTOR.md','selector.py','check.py']}},indent=2))
