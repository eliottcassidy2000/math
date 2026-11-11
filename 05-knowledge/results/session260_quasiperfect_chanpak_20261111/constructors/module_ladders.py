#!/usr/bin/env python3
"""Exact new module-substitution ladders; no inherited computation imports."""
import json
from math import factorial
from pathlib import Path
HERE=Path(__file__).resolve().parent

def endpoints(T):
 n=len(T);dp=[[0]*n for _ in range(1<<n)]
 for v in range(n):dp[1<<v][v]=1
 for S in range(1,1<<n):
  for v,a in enumerate(dp[S]):
   if not a:continue
   left=T[v]&~S
   while left:
    b=left&-left;left-=b;dp[S|b][b.bit_length()-1]+=a
 return dp

def data(T,O):
 n=len(T);full=(1<<n)-1;e=endpoints(T)
 rev=tuple(sum(1<<u for u in range(n) if T[u]>>v&1) for v in range(n));s=endpoints(rev)
 h=sum(e[full]);w=u=0
 for A in range(1<<n):
  B=full^A
  if not B:continue
  bs=sum(s[B][v] for v in range(n) if O>>v&1)
  ae=sum(e[A][v] for v in range(n) if O>>v&1)
  w+=ae*bs;u+=(sum(e[A]) if A else 1)*bs
 return h,w,u

def inflate(T,O,v,m,mark='none'):
 """Replace v by m-chain; m=0 deletes v. Preserve all other marks."""
 old=[i for i in range(len(T)) if i!=v];k=len(old);out=[0]*(k+m)
 for i,x in enumerate(old):
  for j,y in enumerate(old):
   if T[x]>>y&1:out[i]|=1<<j
  if T[x]>>v&1:out[i]|=((1<<m)-1)<<k
 for a in range(m):
  for j,y in enumerate(old):
   if T[v]>>y&1:out[k+a]|=1<<j
  out[k+a]|=sum(1<<(k+b) for b in range(a+1,m))
 sel=sum(1<<i for i,x in enumerate(old) if O>>x&1)
 if m and mark=='first':sel|=1<<k
 if m and mark=='last':sel|=1<<(k+m-1)
 return tuple(out),sel

def stirling(D):
 S=[[0]*(D+1) for _ in range(D+1)];S[0][0]=1
 for m in range(1,D+1):
  for r in range(1,m+1):S[m][r]=S[m-1][r-1]+r*S[m-1][r]
 fact=1
 for r in range(D+1):
  if r:fact*=r
  for m in range(D+1):S[m][r]*=fact
 return S

def coefficients(vals):
 M=stirling(len(vals)-1);ans=[]
 for m,x in enumerate(vals):
  numerator=x-sum(M[m][r]*ans[r] for r in range(m))
  if numerator%M[m][m]:raise RuntimeError('nonintegral skeleton coefficient')
  ans.append(numerator//M[m][m])
 return ans

def make_family(name,T,O,v,mark,expected):
 # Unselected replacement: k=n-1 outside vertices, degree k+2.
 # Selected-first replacement: fix the first clone as outside; replace remaining
 # unselected clones. Parameter t=m-1, k=n, degree k+2.
 D=len(T)+1 if mark=='none' else len(T)+2
 samples=[]
 for t in range(D+1):
  m=t if mark=='none' else t+1
  B,P=inflate(T,O,v,m,mark);h,w,u=data(B,P)
  samples.append((h,h-w,u))
 coef={key:coefficients([s[i] for s in samples]) for i,key in enumerate(['H','D','U'])}
 # Separate direct skeleton count: r mutually nonadjacent twins represent r runs.
 # Divide by r! to remove their freely permutable labels.
 direct={key:[] for key in ['H','D','U']}
 for r in range(D+1):
  B,P=inflate(T,O,v,r if mark=='none' else r+1,mark)
  rows=list(B);start=len(T)-1+(mark!='none');mask=((1<<r)-1)<<start
  for a in range(start,len(rows)):rows[a]&=~mask
  h,w,u=data(tuple(rows),P)
  for key,value in zip(['H','D','U'],[h,h-w,u]):
   if value%factorial(r):raise RuntimeError('skeleton labelling divisibility')
   direct[key].append(value//factorial(r))
 if coef!=direct:raise RuntimeError('independent skeleton coefficient disagreement')
 # Positive-t exponential identity certificate in the surjection basis.
 # Evaluate candidate also at t=0; deletion boundary allowed separately.
 expected_samples=[expected(t) for t in range(D+1)]
 expected_coef={key:coefficients([s[i] for s in expected_samples]) for i,key in enumerate(['H','D','U'])}
 if coef!=expected_coef:raise RuntimeError((name,coef,expected_coef))
 return dict(name=name,rows=T,selected_mask=O,vertex=v,mark=mark,parameter=('chain length' if mark=='none' else 'chain length minus one'),degree_bound=D,skeleton_coefficients=coef,samples=samples)

def main():
 A=(20,1,2,7,14);B=(24,1,3,6,14)
 families=[
 make_family('binary_q14_plus6',A,6,3,'none',lambda t:(3*2**t+3,1,7*2**t+6) if t else (5,1,11)),
 make_family('mixed_q4_3_plus9_2_plus2',A,10,2,'none',lambda t:(3*2**t+3,1,4*3**t+9*2**t+2) if t else (5,-1,15)),
 make_family('mixed_q8_3_plus12_2_minus4',A,24,2,'none',lambda t:(3*2**t+3,1,8*3**t+12*2**t-4) if t else (5,1,16)),
 make_family('ternary_q30_minus1',B,6,1,'first',lambda t:(12*2**t-1,1,30*3**t-1)),
 make_family('ternary_q50_minus3',B,24,3,'first',lambda t:(12*2**t-1,1,50*3**t-3)),
 make_family('relocated_binary_q24_plus16',(40,5,1,6,15,30),10,4,'none',lambda t:(7*2**t+5,1,24*2**t+16) if t else (11,1,37)),
 ]
 out={'status':'PASS','scope':'Exact finite skeleton-coefficient identities imply all nonnegative parameter values; strong seed statements exclude deletion boundaries','families':families}
 (HERE/'module_ladders_receipt.json').write_text(json.dumps(out,indent=2)+'\n')
 print(json.dumps(out,indent=2))
if __name__=='__main__':main()
