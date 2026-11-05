#!/usr/bin/env python3
from verify import *
from fractions import Fraction

def apply(p,w):
 for c in w:p=D(p) if c=='D' else mul(mon(1 if c=='X' else 2),p)
 return p
def moments(w,parity,n):
 p=mon(parity);out=[1]
 for _ in range(n):p=apply(p,w);out.append(value(p,1))
 return out
def matrix(w,parity,size):
 # Extra rows retain every intermediate index for the requested prefix.
 N=size+len(w)+1;T=[[int(i==j) for j in range(N)] for i in range(N)]
 current=parity
 for c in w:
  A=[[0]*N for _ in range(N)]
  if c=='X' and current==1:
   for k in range(N):A[k][k]=1
  else:
   for k in range(N):
    A[k][k]=k*k+6*k+(3 if current==1 else 5) if c=='D' else (k+1)**2
    if k+1<N:A[k+1][k]=1
    if c=='D' and k:A[k-1][k]=2*k*(k if current==1 else k+1)*(2*k+1)
  T=[[sum(A[i][k]*T[k][j] for k in range(N)) for j in range(N)] for i in range(N)]
  if c=='X':current=3-current
 return [row[:size] for row in T[:size]]
def minor_checks(A,maxk):
 total=0
 for k in range(1,maxk+1):
  for rs in combinations(range(len(A)),k):
   for cs in combinations(range(len(A[0])),k):
    require(det([[A[r][c] for c in cs] for r in rs])>=0,'TN finite minor');total+=1
 return total
def sfraction(mu):
 f=list(map(Fraction,mu));out=[]
 while len(f)>1:
  require(f[0]==1,'normalized continued fraction series')
  a=f[1];require(a>=0,'nonnegative S coefficient');out.append(str(a))
  if a==0:
   require(not any(f[1:]),'finite-window termination');break
  inv=[Fraction(1)]
  for n in range(1,len(f)):inv.append(-sum(f[k]*inv[n-k] for k in range(1,n+1)))
  f=[-v/a for v in inv[1:]]
 return out

def physical(w,parity,reps):
 a=transitive(parity)
 for _ in range(reps):
  for c in w:a=wrapper(a) if c=='D' else join(a,transitive(1 if c=='X' else 2))
 return a

def mirror(T,t,j):
 O=transitive(t);blocks=[T,O,tuple(((1<<t)-1)^(1<<v)^O[v] for v in range(t))];starts=[0,len(T),len(T)+t];n=len(T)+2*t;a=[0]*n
 for b,g in enumerate(blocks):
  for i,row in enumerate(g):a[starts[b]+i]|=row<<starts[b]
 for r in range(len(T)):
  a[r]|=((1<<t)-1)<<len(T)
  for i in range(t):a[len(T)+t+i]|=1<<r
 for u in range(t):
  for v in range(t):
   if u==j and v==j:a[len(T)+u]|=1<<(len(T)+t+v)
   else:a[len(T)+t+v]|=1<<(len(T)+u)
 validate(a);return tuple(a)
def run():
 words=[''.join(w) for k in range(5) for w in product('PD',repeat=k)];hm=pm=direct=0;sfr=[]
 for parity in (1,2):
  for w in words:
   mu=moments(w,parity,10);H=[[mu[i+j] for j in range(6)] for i in range(6)]
   hm+=minor_checks(H,3);pm+=minor_checks(matrix(w,parity,5),3)
   sf=sfraction(mu)
   if w in ('D','PD','DP'):sfr.append({'parity':parity,'word':w,'moments':mu[:5],'S_coefficients':sf})
   if len(w)<=3:
    for m in (1,2):
     a=physical(w,parity,m);require(h_subsets(a)[-1]==mu[m] and len(a)==parity+2*len(w)*m,'literal repeated physical word');direct+=1
 xwords=[''.join(w) for k in range(5) for w in product('DX',repeat=k) if w.count('X')%2==0]
 xdirect=0;xm=0;xpm=0
 for parity in (1,2):
  for w in xwords:
   mu=moments(w,parity,10);H=[[mu[i+j] for j in range(6)] for i in range(6)]
   xm+=minor_checks(H,3);xpm+=minor_checks(matrix(w,parity,5),3);sfraction(mu)
   if len(w)<=3:
    for m in (1,2):
     a=physical(w,parity,m)
     require(h_subsets(a)[-1]==mu[m] and len(a)==parity+m*(2*w.count('D')+w.count('X')),'literal parity-returning word')
     xdirect+=1
 mirror_checks=0
 for t in range(1,4):
  for j in range(t):
   for T in ((0,),transitive(2),(2,4,1)):
    a=mirror(T,t,j);p=bin_polynomial(T);want=mul(mon(2*(t-1-j)),D(mul(mon(2*j),p)))
    require(bin_polynomial(a)==want,'full physical transitive mirror identity');mirror_checks+=1
 require(moments('PD',2,2)==[1,17,913],'periodically repeated padding control')
 p=mon(4);bad=[value(p,1)]
 for _ in range(2):p=D(p);bad.append(value(p,1))
 require(bad==[1,17,241] and bad[0]*bad[2]-bad[1]**2==-48,'one-time-prefix hostile')
 print(json.dumps({'status':'PASS','words_per_parity':len(words),'Hankel_minors':hm,'production_minors':pm,'physical_word_graphs':direct,'physical_mirror_profiles':mirror_checks,'even_X_words_per_parity':len(xwords),'even_X_Hankel_minors':xm,'even_X_production_minors':xpm,'even_X_physical_graphs':xdirect,'examples':sfr,'one_time_prefix_hostile':bad},sort_keys=True,indent=2))
if __name__=='__main__':run()
