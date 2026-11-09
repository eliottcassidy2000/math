import json,math,random,itertools
from functools import lru_cache
from pathlib import Path
A=Path(__file__).resolve().parent

def ck(p,m):
 if not p:raise RuntimeError(m)
def val(w):return int(w,2) if w else 0
def rev(w):return val(w[::-1])
def count(w):
 # Direct positional Hamilton-count formula, independent of g recurrence.
 if not w:return 1
 return int(w[0]=='1')+int(w[-1]=='0')+sum(w[i:i+2]=='01' for i in range(len(w)-1))+sum(1<<(i-j-1) for i in range(len(w)) for j in range(i) if w[j]=='1' and w[i]=='0')
def F(w,a,b,c):return a*count(w)-b*(1-val(w)%2)+c*(val(w)+1)
def f(n):
 if not n:return 0
 u=n//2
 return f(u)+(u+u%2 if n%2 else 0)
def hp(w):
 n=len(w)+1
 def edge(i,j):
  if i==0:return w[j-1]=='1'
  if j==0:return w[i-1]=='0'
  return i<j
 @lru_cache(None)
 def paths(last,remaining):
  if not remaining:return 1
  return sum(paths(j,remaining^(1<<j)) for j in range(n) if remaining>>j&1 and edge(last,j))
 return sum(paths(i,((1<<n)-1)^(1<<i)) for i in range(n))
words=['']+[''.join(x) for k in range(1,4) for x in itertools.product('01',repeat=k)]
rows=[(1,0,1),(1,1,1),(3,1,4),(5,99,1),(189,45,432),(2,0,2)]
rng=random.Random(23102)
configs=[(u,w,s,b,g) for u in words[:7] for w in words for s in words for b in [0,1] for g in [0,1]]
configs += [(''.join(rng.choices('01',k=rng.randrange(18))),''.join(rng.choices('01',k=rng.randrange(18))),''.join(rng.choices('01',k=rng.randrange(18))),rng.randrange(2),rng.randrange(2)) for _ in range(1000)]
coeffs=rigidity=0
for U,W,S,beta,gamma in configs:
 P=val(U);M=1<<len(W);r=val(W);R=rev(W);C=M-1-R
 vals={s:U+str(beta)*s+W+str(gamma)*s+S for s in [1,2,3,5]}
 for a,b,c in rows:
  K=c*(1<<len(S))+a*((1<<len(S))-1-rev(S))
  if beta==gamma==1:A2=K*M*(P+1);A1=a*C*(P+1)-K*(M-1-r)
  elif beta==gamma==0:A2=(a+K)*M*P;A1=(a+K)*r-a*P*R
  elif beta==1:A2=(a+K)*M*(P+1);A1=-a*(R+1)*(P+1)-(a+K)*(M-r)
  else:A2=K*M*P;A1=a*P*(1+C)+K*(r+1)
  constant=F(vals[1],a,b,c)-4*A2-2*A1
  for s in [2,3,5]:ck(F(vals[s],a,b,c)==A2*4**s+A1*2**s+constant,'coefficient table')
  coeffs+=1
  if a%2 and A2>0 and A2&(A2-1)==0 and A1==0:
   ck(beta==gamma,'mixed rigidity')
   if a>1:ck(W==str(beta)*len(W),'separated rigidity')
   rigidity+=1
# Independent finite reversal-lemma control, including powers of two fractions.
reversal=0
for k in range(1,11):
 for n in range(1,1<<k):
  R=rev(format(n,'0%db'%k))
  for a in range(1,32,2):
   x=a*R;y=n
   while x%2==0:x//=2
   while y%2==0:y//=2
   if x==y:ck(a==1,'odd reversal lemma')
   reversal+=1
# Explicitly show odd-a assumption matters: a=2 and W=00 is separated success.
for s in range(1,15):ck(F('1'*s+'00'+'1'*s,2,0,2)==2**(2*s+3)-2,'even a countercontrol')
half=graphs=0
for l in range(0,9):
 for q in range(1<<l):
  x=(1<<l)-1-q;left=format(x,'0%db'%l) if l else ''
  for beta in [None,0,1]:
   w=left+('' if beta is None else str(beta))+left[::-1]
   d=1-q*q-2*f(q)-q%2 if beta is None else 1-2*(q*q+f(q)+(1-beta)*(q+q%2))
   for s in [1,3,7]:
    bits='1'*s+w+'1'*s+'0';N=2**(2*s+len(w))+d
    ck(count(bits)==N,'half constant');ck(N>0 and N%2,'physical odd');half+=1
   if l<=2:
    bits='1'+w+'10';ck(hp(bits)==2**(2+len(w))+d,'literal Hamilton count');graphs+=1
# Decoder retains exponent parity/length and rejects positive offset or bad parity.
def decode(m,d):
 if d>1:return []
 if m%2==0:
  T=1-d;q=math.isqrt(T)
  return [(q,None)] if T==q*q+2*f(q)+q%2 and m>=2*q.bit_length()+2 else []
 if d%2==0:return []
 T=(1-d)//2;q=math.isqrt(T)
 return [(q,b) for b in [0,1] if T==q*q+f(q)+(1-b)*(q+q%2) and m>=2*q.bit_length()+3]
queries=0
for m in range(2,19):
 support={}
 for q in range(1<<max(0,(m-2)//2)):
  if m%2==0: ds=[(1-q*q-2*f(q)-q%2,None)]
  else:ds=[(1-2*(q*q+f(q)+(1-b)*(q+q%2)),b) for b in [0,1]]
  for d,b in ds:support.setdefault(d,set()).add((q,b))
 for d in range(-2000,4):
  ck(set(decode(m,d))==support.get(d,set()),'clock-retaining decoder');queries+=1
 ck(not decode(m,-9),'deficit nine')
print(json.dumps(dict(status='PASS',coefficient_row_cases=coeffs,pure_dyadic_controls=rigidity,reversal_lemma_controls=reversal,halfword_checks=half,literal_graph_checks=graphs,clock_retaining_decoder_queries=queries,scope='Two-clock rigidity under odd a only; inherited centered fans with refined decoder; no coverage conclusion.'),indent=2))
