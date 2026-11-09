import json,sys,hashlib
from pathlib import Path
OUT=Path(__file__).resolve().parent
C=[4,1,18,7,75,31,175,63]; MASK=124

def literal(rows):
 n=len(rows);d=[[0]*n for _ in range(1<<n)];h=[1]+[0]*((1<<n)-1)
 for s in range(1,1<<n):
  for v in range(n):
   if s>>v&1:
    t=s^(1<<v);d[s][v]=sum(d[t][u] for u in range(n) if rows[u]>>v&1) if t else 1
  h[s]=sum(d[s])
 return h,d

def corestate():
 h,e=literal(C); rev=[sum((C[j]>>i&1)<<j for j in range(8)) for i in range(8)];_,s=literal(rev);full=255
 starts=lambda t:sum(s[t][u] for u in range(8) if MASK>>u&1)
 ends=lambda t:sum(e[t][u] for u in range(8) if MASK>>u&1)
 return (h[-1],sum(h[t]*h[full^t] for t in range(256)),sum(h[t]*starts(full^t) for t in range(256)),sum(ends(t)*h[full^t] for t in range(256)),sum(ends(t)*starts(full^t) for t in range(256)),starts(full),ends(full))
def comp(x,y):
 h,q,A,G,T,S,E=x;g,r,B,H,U,V,W=y
 return h*g,q*r,A*r+h*B,H*q+g*G,A*H+T*g+h*U+E*V,S*g,h*W

def word(p,b=None):
 if b is None:b=p.bit_length()
 z=None
 for i in range(b-1,-1,-1):
  c=p>>i&1;v=(1,2,c,c,0,c,c);z=v if z is None else comp(z,v)
 return z

def data(k,j):
 r=4096**k;J=(r-1)//4095;A=878*r+171630*J;G=577+2383200*J;q=1440*r
 Z=(13617152*r+1602510*k+717259)//24843
 S=CS[5] if j==0 else 45;E=14 if j==k else 0;T=A+45-E-Z+2*j
 return (45,q,A,G,T,S,E),Z

def value(k,j,P,R,b):
 x,Z=data(k,j);h,q,A,G,T,S,E=x;p=word(P);w=word(R,b)
 fp=0 if p is None else p[4];fr=0 if w is None else w[4];rev=0 if w is None else w[3]
 return (q*P+A)*((1<<b)-rev)+45*(R+1-fr-(R&1 if b else 0))-P*G-45*fp-T-(P&1)*S-(E*((R>>(b-1))&1) if b else E)
CS=corestate()
def check(ok,msg):
 if not ok:raise RuntimeError(msg)
def controls():
 check(CS[:4]==(45,1440,878,577),'core state')
 for k in range(3):
  for j in range(k+1):
   x,Z=data(k,j);check(value(k,j,0,0,0)==Z-2*j,'bare interval')
   for P in range(8):
    for b in range(4):
     for R in range(1<<b):
      y=x
      if P:y=comp(word(P),y)
      if b:y=comp(y,word(R,b))
      check(value(k,j,P,R,b)==y[2]-y[4]+y[0]-y[6],'context formula')
 for P in range(4):
  for b in range(3):
   for R in range(1<<b):
    a=P.bit_length();bits=[P>>i&1 for i in range(a-1,-1,-1)]+[MASK>>i&1 for i in range(8)]+[R>>i&1 for i in range(b-1,-1,-1)]
    n=len(bits);rows=[sum(1<<v for v in range(u+1,n)) for u in range(n)]
    for u in range(8):rows[a+u]=(rows[a+u]&~(255<<a))|(C[u]<<a)
    rows=[x|((not bits[u])<<n) for u,x in enumerate(rows)]+[sum(c<<u for u,c in enumerate(bits))]
    check(literal(rows)[0][-1]==value(0,0,P,R,b),'literal control')
 # Absorb X/Y into contexts: identity of the physical cut-labelled graph.
 for k in range(3):
  for j in range(k+1):
   for P,R,b in [(0,0,0),(1,2,2),(2,0,1)]:
    x=int('101000111000'*j,2) if j else 0;y=int('111011100110'*(k-j),2) if k>j else 0
    pp=(P<<(12*j))|x;rr=(y<<b)|R;bb=b+12*(k-j)
    check(value(k,j,P,R,b)==value(0,0,pp,rr,bb),'context absorption')
 # Neutral source/sink padding must not create new occupied targets.
 for k in range(2):
  for j in range(k+1):
   for P,R,b in [(0,0,0),(1,2,2),(2,0,1)]:
    n=value(k,j,P,R,b)
    check(value(k,j,P,2*R+1,b+1)==n,'neutral selected suffix sink')
    x,_=data(k,j);p=word(P);w=word(R,b)
    if p:x=comp(p,x)
    if w:x=comp(x,w)
    x=comp((1,2,0,0,0,0,0),x)
    check(x[2]-x[4]+x[0]-x[6]==n,'neutral unselected prefix source')
 hostile={str((P,R,b)):[value(1,j,P,R,b) for j in range(2)] for P,R,b in [(0,0,0),(1,0,0),(0,2,2),(2,0,1)]}
 return {'core_state':CS,'hostile_context_values':hostile,'status':'PASS'}

def census(cap):
 # target->first witness. One byte per odd target stores k=0/k>=1 membership.
 flags=bytearray((cap+1)//2);witnesses={};reps={};k=0
 while data(k,0)[1]-2*k<=cap:
  reps[k]=0
  for j in range(k+1):
   (h,q,A,G,T,S,E),Z=data(k,j)
   for b in range(0,(cap//A).bit_length()+1):
    if b==0:bound=(cap-Z+2*j+S)//(q-G-45)
    else:bound=((cap-45)//(1<<(b-1))-A)//q
    if bound<0:continue
    for P in range(bound+1):
     p=word(P);fp=p[4] if p else 0;base=q*P+A;fixed=-P*G-45*fp-T-(P&1)*S
     for R in (range(0,1<<b,2) if b else [0]):
      w=word(R,b);fr=w[4] if w else 0;rev=w[3] if w else 0
      N=base*((1<<b)-rev)+45*(R+1-fr)+fixed-(E*(R>>(b-1)&1) if b else E)
      check(N%2==1 and N>0,'odd positive')
      if N>cap:continue
      idx=N//2;flags[idx]|=1 if k==0 else 2;reps[k]+=1
      if k or len(witnesses)<20:witnesses.setdefault(N,[k,j,P,R,b,12*k+9+P.bit_length()+b])
  k+=1
 # Deterministic set-level receipt: odd n maps to bit floor(n/2), LSB-first.
 packed=bytearray((len(flags)+7)//8)
 for i,f in enumerate(flags):
  if f:packed[i//8]|=1<<(i%8)
 receipt={'cap':cap,'encoding':'LSB-first bit floor(n/2) for attained positive odd n <= cap; zero otherwise','bytes':len(packed),'sha256':hashlib.sha256(packed).hexdigest()}
 (OUT/'UNION.bin').write_bytes(packed)
 (OUT/'UNION.sha256.json').write_text(json.dumps(receipt,indent=2)+'\n')
 def stats(limit):
  count=only=high=both=0;prev=None;gaps=[];first=[];last=[]
  for i,f in enumerate(flags[:(limit+1)//2]):
   if not f:continue
   n=2*i+1;count+=1;only+=f==2;high+=bool(f&2);both+=f==3
   if len(first)<20:first.append(n)
   if prev is not None and n-prev>2:gaps.append((n-prev-2,prev,n))
   prev=n;last.append(n)
   if len(last)>10:last.pop(0)
  return dict(cap=limit,attained=count,missing=(limit+1)//2-count,missing_fraction=((limit+1)//2-count)/((limit+1)//2),k_ge_1=high,new_k_ge_1=only,overlap=both,largest_internal_gaps=sorted(gaps,reverse=True)[:10],first=first,last=last)
 result={'controls':controls(),'cap':cap,'representations_by_k':reps,'statistics':[stats(x) for x in [10000000,cap] if x<=cap],'witness_format':['k','j','prefix_integer','suffix_integer','suffix_length','order'],'witnesses':witnesses}
 (OUT/'RESULT.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:v for k,v in result.items() if k!='witnesses'},indent=2))
if __name__=='__main__':
 if len(sys.argv)>1 and sys.argv[1]=='controls':print(json.dumps(controls(),indent=2))
 else:census(int(sys.argv[1]) if len(sys.argv)>1 else 100000000)
