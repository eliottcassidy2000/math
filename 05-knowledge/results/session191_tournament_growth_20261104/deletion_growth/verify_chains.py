#!/usr/bin/env python3
"""Independent literal subset counts for the interval-chain deletion formulas."""
import itertools,json,pathlib,sys
P=pathlib.Path(__file__).resolve().parent
MAX=int(sys.argv[1]) if len(sys.argv)>1 else 12

def require(b,s):
 if not b: raise RuntimeError(s)
def comps(n):
 if n==0:yield ()
 for k in range(2,n+1):
  for rest in comps(n-k):yield(k,)+rest
def chain(parts):
 n=sum(parts)+1;rows=[sum(1<<j for j in range(i+1,n))for i in range(n)]
 left=0
 for p in parts:
  rows[left]^=1<<(left+p);rows[left+p]^=1<<left;left+=p
 return rows
def counts(rows):
 n=len(rows);D=[[0]*n for _ in range(1<<n)];H=[0]*(1<<n);H[0]=1
 for S in range(1,1<<n):
  for end in range(n):
   if not(S>>end&1):continue
   old=S^(1<<end)
   D[S][end]=1 if not old else sum(D[old][j]for j in range(n)if old>>j&1 and rows[j]>>end&1)
  H[S]=sum(D[S])
 return H
def strong(rows,S):
 if S.bit_count()<3:return False
 for start in range(len(rows)):
  if not(S>>start&1):continue
  seen=1<<start
  while True:
   nxt=seen
   for i in range(len(rows)):
    if seen>>i&1:nxt|=rows[i]&S
   if nxt==seen:break
   seen=nxt
  if seen!=S:return False
 return True
def prod(xs):
 a=1
 for x in xs:a*=x
 return a
def f(p):return 2**(p-1)+1
records=[];mins={}
for n in range(4,MAX+1):
 best=None;witness=None
 for cc in comps(n-1):
  rows=chain(cc);H=counts(rows);full=(1<<n)-1
  require(H[-1]==prod(map(f,cc)),'Moon product')
  start=0;interiors={}
  for k,p in enumerate(cc):
   for v in range(start+1,start+p):interiors[v]=(k,p)
   start+=p
  for v in range(n):
   isstrong=strong(rows,full^(1<<v));require(isstrong==(v in interiors),'strong deletion classification')
   if not isstrong:continue
   k,p=interiors[v]
   if p>=3: predicted=2**(p-2)*prod(f(q)for j,q in enumerate(cc)if j!=k)
   else:
    left,right=cc[:k],cc[k+1:];a=prod(map(f,left));b=prod(map(f,right))
    ea=prod(map(f,left[:-1])) if left else 1
    sb=prod(map(f,right[1:])) if right else 1
    predicted=a*b+ea*sb
   delta=H[-1]-H[full^(1<<v)];require(delta==predicted,'marked block formula')
   records.append({'parts':cc,'root':v,'H':H[-1],'Hminusroot':H[full^(1<<v)],'Delta':delta})
   if best is None or delta<best:best=delta;witness=records[-1]
 mins[n]={'minimum':best,'witness':witness}
# Exact finite minimization is additionally checked over many orders without graph DP.
def M(t):
 if t<0 or t==1:return None
 if t==0:return 1
 if t%3==0:return 5**(t//3)
 if t%3==1:return 9*5**((t-4)//3)
 return 3*5**((t-2)//3)
coeffs=[(2,3),(4,4),(6,5),(10,6),(16,7),(26,8),(46,9),(82,10)]
formulas=[]
for n in range(4,401):
 ans=min(c*M(n-1-cost)for c,cost in coeffs if M(n-1-cost)is not None)
 if n in mins:require(ans==mins[n]['minimum'],'finite minimum formula')
 if n>=9:
  seed={0:(9,26),1:(10,46),2:(11,78)}[n%3]
  require(ans==seed[1]*5**((n-seed[0])//3),'period three exact formula')
 formulas.append({'n':n,'minimum':ans})
receipt={'status':'PASS','maximum_literal_order':MAX,'literal_marked_pairs':len(records),'minima':mins,'records':records,'finite_minimization_through400':formulas}
(P/'chain_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps({'status':'PASS','maximum_literal_order':MAX,'literal_marked_pairs':len(records),'minima':mins},sort_keys=True))
