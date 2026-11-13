#!/usr/bin/env python3
import json,itertools,pathlib
checks=0
def require(x):
 global checks; checks+=1
 if not x:raise RuntimeError('check failed')
def h(s):return 1+2*sum((1 if j-i<=2 else 2**(j-i-2)) for i in range(len(s)) for j in range(i+1,len(s)) if s[i]=='1' and s[j]=='0')
def dp(s):
 n=len(s)+1; z=n-1; rows=[0]*n
 for i in range(z):
  for j in range(i+1,z):rows[i]|=1<<j
  if s[i]=='1':rows[z]|=1<<i
  else:rows[i]|=1<<z
 D=[[0]*n for _ in range(1<<n)]
 for i in range(n):D[1<<i][i]=1
 for S in range(1,1<<n):
  for v in range(n):
   if D[S][v]:
    q=rows[v]&~S
    while q:
     bit=q&-q;q-=bit;D[S|bit][bit.bit_length()-1]+=D[S][v]
 return sum(D[-1])
for n in range(1,8):
 for bits in itertools.product('01',repeat=n):
  s=''.join(bits);require(h(s)==dp(s))
rows=[]
for m in range(2,15):
 plus=[];zero=0
 for bits in itertools.product('01',repeat=m):
  w=''.join(bits);x='10'+w+'01';y='01'+w+'10'
  B=int(w[:-1],2);C=int(''.join('1' if c=='0' else '0' for c in w[:0:-1]),2)
  d=h(y)-h(x);require(d==2*(B-C))
  anti=all(w[i]!=w[m-1-i] for i in range(m))
  require((d==0)==anti)
  if d==0:zero+=1
  if d==2:plus.append(w)
 require(zero==(2**(m//2) if m%2==0 else 0))
 expect=['11'] if m==2 else ['011','101','110'] if m==3 else ['1'+'0'*(m-2)+'1']
 require(plus==expect)
 if m>=3:
  w='1'+'0'*(m-2)+'1';require(h('10'+w+'01')==15*2**(m-2)+1)
 rows.append({'m':m,'zero':zero,'plus_words':plus})
# Very long exact construction checks, not floating fits.
for m in [31,64,127,200]:
 w='1'+'0'*(m-2)+'1';require(h('10'+w+'01')==15*2**(m-2)+1);require(h('01'+w+'10')==15*2**(m-2)+3)
print(json.dumps({'checks':checks,'classification':rows},indent=2))
