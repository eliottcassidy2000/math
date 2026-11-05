#!/usr/bin/env python3
import json, math, random
from functools import reduce

def need(x,s):
 if not x:raise RuntimeError(s)
def trim(a):
 a=list(a)
 while len(a)>1 and a[-1]==0:a.pop()
 return a
def add(*args):
 r=[0]*max(map(len,args))
 for a in args:
  for j,x in enumerate(a):r[j]+=x
 return trim(r)
def scale(a,c):return trim([c*x for x in a])
def mul(a,b):
 c=[0]*(len(a)+len(b)-1)
 for i,x in enumerate(a):
  for j,y in enumerate(b):c[i+j]+=x*y
 return trim(c)
def val(p,x):
 a=0
 for c in reversed(p):a=a*x+c
 return a
def shift(p,z):
 q=[0]*len(p)
 for i,c in enumerate(p):
  for j in range(i+1):q[j]+=c*math.comb(i,j)*z**(i-j)
 return trim(q)
def D(p):return add([0,0]+p,[0]+add(shift(p,1),scale(shift(p,-1),-1)))
def basis(k,odd):
 p=[0,1]if odd else[0,0,1]
 for j in range(1,k+1):p=mul(p,[-j*j,0,1])
 return p
def coords(p):
 odd=(len(p)-1)%2;u=[0]*((len(p))//2);q=p[:]
 for k in range(len(u)-1,-1,-1):
  b=basis(k,odd);u[k]=q[len(b)-1]if len(q)>=len(b)else 0;q=add(q,scale(b,-u[k]))
 need(q==[0],'central recovery');return u
def gamma(k,odd):return math.factorial(2*k+1)*(1 if odd else k+1)
def beta(k,odd):return 2*k*(k if odd else k+1)*(2*k+1)
def newton(p):
 v=[val(p,x)for x in range(len(p))];f=[]
 while v:f.append(v[0]);v=[b-a for a,b in zip(v,v[1:])]
 return f
def contents(p):
 f=newton(p);g=reduce(math.gcd,f);G=reduce(math.gcd,(math.factorial(j-1)*f[j]for j in range(1,len(f))));C=0;n=len(f)-1;f+=[0,0]
 for r in range(1,n+2):
  C=math.gcd(C,math.factorial(r)**2*(f[r-1]+3*f[r]+f[r+1]));C=math.gcd(C,math.factorial(r)*math.factorial(r+1)*(f[r]+f[r+1]))
 return g,G,C
def valuation(x,p=3):
 need(x!=0,'nonzero valuation');v=0
 while x%p==0:x//=p;v+=1
 return v

def tournament(n,c):
 a=[0]*n;b=0
 for i in range(n):
  for j in range(i+1,n):
   x,y=(i,j)if c>>b&1 else(j,i);a[x]|=1<<y;b+=1
 return tuple(a)
def trans(n):return tuple(sum(1<<j for j in range(i+1,n))for i in range(n))
def validate(a):
 n=len(a)
 need(all(x>=0 and x<1<<n and not x>>i&1 for i,x in enumerate(a)),'valid rows')
 need(all(((a[i]>>j)^(a[j]>>i))&1 for i in range(n)for j in range(i+1,n)),'tournament')
def join(a,b):
 n=len(a);return tuple(x|(((1<<len(b))-1)<<n)for x in a)+tuple(x<<n for x in b)
def cyclic(a,b,c):
 n,m,k=map(len,(a,b,c));return tuple(x|(((1<<m)-1)<<n)for x in a)+tuple((x<<n)|(((1<<k)-1)<<(n+m))for x in b)+tuple((x<<(n+m))|((1<<n)-1)for x in c)
def W(a):return cyclic(a,(0,),(0,))
def E(a):return cyclic(a,trans(2),(0,))
def strong(a):
 n=len(a)
 for i in range(n):
  reach=1<<i
  while True:
   z=reach
   for j in range(n):
    if reach>>j&1:z|=a[j]
   if z==reach:break
   reach=z
  if reach!=(1<<n)-1:return False
 return True
def pathcounts(a):
 n=len(a);N=1<<n;dp=[[0]*n for _ in range(N)];h=[0]*N;h[0]=1
 for i in range(n):dp[1<<i][i]=1
 for S in range(1,N):
  h[S]=sum(dp[S])
  for v,x in enumerate(dp[S]):
   if x:
    q=a[v]&~S
    while q:
     b=q&-q;q-=b;dp[S|b][b.bit_length()-1]+=x
 return h
def profile(a):
 n=len(a);h=pathcounts(a);c=[[0]*(n+1)for _ in h];c[0][0]=1
 for S in range(1,len(h)):
  least=S&-S;s=S
  while s:
   if s&least:
    for r in range(1,n+1):c[S][r]+=h[s]*c[S^s][r-1]
   s=(s-1)&S
 p=[0];fall=[1]
 for r,z in enumerate(c[-1]):p=add(p,scale(fall,z));fall=mul(fall,[-r,1])
 return p

def run():
 columns=formal=physical=strong_outputs=0;rng=random.Random(1971513)
 matrices={}
 for odd in (False,True):
  for d in range(1,9):
   L=[[0]*d for _ in range(d)]
   for k in range(d):
    p=basis(k,odd);factor=math.factorial(k)*gamma(k,odd);prod=math.prod(beta(i,odd)for i in range(1,k+1));need(prod==factor,'lowering product')
    expected=add(basis(k+1,odd),scale(p,k*k+6*k+(3 if odd else 5)),scale(basis(k-1,odd),beta(k,odd))if k else[0])
    need(D(p)==expected,'literal tridiagonal polynomial identity')
    for j in range(2*d+2):
     x=val(p,1);need(x%factor==0 and(j>=k or x==0)and(j!=k or x==factor),'all paths divisible and unit boundary')
     if j<d:L[j][k]=x//factor
     p=D(p)
    columns+=1
   matrices[odd,d]=L
   for trial in range(30):
    u=[rng.randrange(-20,21)for _ in range(d)];u[0]=2*rng.randrange(-20,21)+1;u[-1]=u[-1]or 1
    p=add(*(scale(basis(k,odd),x)for k,x in enumerate(u)));expected=[math.factorial(k)*gamma(k,odd)*x for k,x in enumerate(u)];G=reduce(math.gcd,expected,0);out=[];q=p
    for j in range(max(d,2)):out.append(val(q,1));q=D(q)
    need(reduce(math.gcd,out[:d],0)==G,'triangle bank gcd')
    recovered=[]
    for j in range(d):recovered.append(out[j]-sum(L[j][k]*recovered[k]for k in range(j)))
    need(recovered==expected,'integer inverse response reconstruction')
    strong_bank=out[1:max(1,d-1)+1]+[3*val(p,2)-val(p,1)]
    need(reduce(math.gcd,strong_bank)==G,'strong-only finite bank');formal+=1
 for n in range(1,6):
  for code in range(1<<(n*(n-1)//2)):
   a=tournament(n,code);p=profile(a);d=(n+1)//2;g,G,C=contents(p);x=a;values=[pathcounts(a)[-1]]
   for j in range(1,max(1,d-1)+1):
    x=W(x);need(strong(x),'actual repeated wrapper strong');values.append(pathcounts(x)[-1]);strong_outputs+=1
   e=E(a);need(strong(e),'actual E strong');hE=pathcounts(e)[-1];strong_outputs+=1
   need(3*values[1]-hE==4*values[0],'actual D E elimination')
   need(reduce(math.gcd,values[:d])==G and reduce(math.gcd,values[1:]+[hE])==G,'actual both finite banks')
   q=p
   for h in values:need(val(q,1)==h,'actual repeated graph count');q=D(q)
   physical+=1
 need(pathcounts(W((0,)))[-1]==3 and pathcounts(W(W((0,))))[-1]==15 and pathcounts(E((0,)))[-1]==5,'identity-exclusion boundary')
 V=(382,420,88,706,810,972,18,597,204,327);validate(V);need(strong(V),'literal V strong')
 P=profile(V);need(P==[0,0,2772,0,4020,0,1000,0,64,0,1]and contents(P)==(9,27,81),'literal physical input profile')
 ev=[val(P,x)for x in range(1,6)];need(ev==[7857,156816,1558521,10412352,52972425],'small evaluations')
 need([valuation(x)for x in ev]==[4,4,4,2,2]and ev[3]//9%3==ev[4]//9%3==2,'valuation residues')
 powers=[];p=[1]
 for m in range(1,19):
  p=mul(p,P);raw=contents(p);wrapped=contents(D(p));need(raw==(3**(2*m),3**(2*m+1),3**(2*m+2)),'all ordered-power contents')
  need(tuple(valuation(x)for x in wrapped)==(2*m,2*m+1,2*m+1),'strong wrapper valuations')
  Z=12*ev[3]**m+46*ev[2]**m-3*ev[1]**m-21*ev[0]**m
  q=D(p);need(Z==4*val(q,3)+5*val(q,2)-11*val(q,1)and valuation(Z)==2*m+1,'second cyclic escape formula')
  need(val(q,1)==ev[0]**m+ev[1]**m,'strong companion scalar ray')
  powers.append({'m':m,'join_contents':raw,'wrapped_v3':tuple(valuation(x)for x in wrapped),'escape_v3':valuation(Z)})
 for m in range(1,5):
  a=()
  for _ in range(m):a=join(a,V)
  s=W(a);z=cyclic(trans(2),trans(2),s);validate(a);validate(s);validate(z);need(strong(a)==(m==1)and strong(s)and strong(z),'physical join and wrapper scope')
 s=W(V);z=cyclic(trans(2),trans(2),s);hs=pathcounts(s)[-1];hz=pathcounts(z)[-1]
 need(hs==164673 and hz==196004745 and hz%81==54,'literal twelve and sixteen vertex counts')
 print(json.dumps({'status':'PASS','normalized_response_columns':columns,'signed_formal_reconstructions':formal,'physical_input_profiles':physical,'actual_strong_outputs':strong_outputs,'identity_boundary':[1,3,15,5],'power_controls':powers,'literal_witness':{'V':V,'P':P,'strong_wrapper_rows':s,'strong_wrapper_H':hs,'second_cyclic_rows':z,'second_cyclic_H':hz,'second_cyclic_mod81':hz%81}},indent=2,sort_keys=True))
if __name__=='__main__':run()
