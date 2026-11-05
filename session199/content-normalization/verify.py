from pathlib import Path
from functools import reduce
from math import gcd
import hashlib,json

def need(ok,msg):
 if not ok:raise RuntimeError(msg)

def move(s,c):
 D,E,U,d=s
 if c=='O':return D,0,2*U,d
 if c=='Z':return D-E+U,d,2*U-d,d
 return D-E+U,D-E,D+U,D

def bank(cap):
 base=(1,0,2,1);out={base:''};queue=[base]
 for s in queue:
  for c in 'OZJ':
   q=move(s,c);need(q[2]>s[2],'strict increase')
   if q[2]<=cap and q not in out:out[q]=out[s]+c;queue.append(q)
 return out

def graph(word):
 T=[2,0];r=0
 for c in word:
  n=len(T);T.append(0)
  for v in range(n):
   if c!='O' and v==r:T[n]|=1<<v
   else:T[v]|=1<<n
  if c=='J':r=n
 return T,r

def count(T,r):
 n=len(T);F=(1<<n)-1;end=[[0]*n for _ in range(F+1)];root=[[0]*n for _ in range(F+1)]
 for i in range(n):end[1<<i][i]=1
 root[1<<r][r]=1
 for S in range(1,F+1):
  for j in range(n):
   if not S>>j&1:continue
   P=S^(1<<j)
   if P:
    end[S][j]=sum(end[P][i]for i in range(n)if P>>i&1 and T[i]>>j&1)
    root[S][j]=sum(root[P][i]for i in range(n)if P>>i&1 and T[i]>>j&1)
 H=[1]+[sum(row)for row in end[1:]]
 U=sum(H[S]*sum(root[F^S])for S in range(F+1))
 return (H[F],end[F][r],U,H[F^(1<<r)]),sum(root[F])

B=bank(39);need(len(B)==97,'complete cap39 count')
rows=[];fails=[]
for s,w in B.items():
 c=reduce(gcd,s)
 if c>1:
  q=tuple(x//c for x in s);physical=q in B
  rows.append(dict(profile=s,word=w,content=c,quotient=q,physical=physical))
  if not physical:fails.append((s,q))
need(len(rows)==17 and fails==[((39,9,39,9),(13,3,13,3))],'first failure')
T,r=graph('OOZOJ');profile,starts=count(T,r)
need(profile==(39,9,39,9)and starts==1 and len(T)==7,'literal witness')
s=profile
for k in range(1,101):
 s=move(s,'O');need(s==(39,0,39*2**k,9)and reduce(gcd,s)==3,'infinite-family formula')
D,E,U,d=(13,3,13,3)
need(min(D-E-d,U-D+E-d,E,d-E)>=0 and D%2==d%2==1 and (U-E)%2==0,'formal cone and parity')
need({s[2]for s in B if s[0]==3}=={3,6,12,24},'D3 controls')
here=Path(__file__).resolve().parent
out=dict(status='PASS',complete_cover_cap=39,profile_count=len(B),nonprimitive_count=len(rows),
         quotient_tests=rows,literal_witness=dict(rows=T,root=r,profile=profile,starts=starts),
         family_formula_checks=100,dependency_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest()for p in sorted((here/'dependencies').glob('*.md'))},
         scope='Primitive response division is not physical; no scalar noncoverage claim.')
print(json.dumps(out,indent=2,sort_keys=True))
