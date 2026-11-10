#!/usr/bin/env python3
"""Synchronized physical outer blocks; standalone exact compiler and verifier."""
import argparse,hashlib,itertools,json
from pathlib import Path
T=(1,0,2,1)
def need(ok,msg):
 if not ok:raise RuntimeError(msg)
def step(v,c):
 D,E,U,d=v
 if c=='O':return D,0,2*U,d
 if c=='Z':return D-E+U,d,2*U-d,d
 if c=='J':return D-E+U,D-E,D+U,D
 raise ValueError(c)
def state(w,v=T):
 for c in w:v=step(v,c)
 return v
def coefficients(w,s=1):
 need(s>=1,'positive final O padding')
 b=state(w,(1,0,0,1));a=state(w,(0,0,2,0));z=2**s
 return z*a[2],a[0]+z*b[2],b[0]
def recipe(w,k,s=1):
 need(k>=0 and s>=1,'padding guards')
 return 'O'*k+w+'O'*(k+s)+'J'
def compile_graph(w):
 rows=[2,0];root=0
 for c in w:
  n=len(rows);rows.append(0)
  for v in range(n):
   if c!='O' and v==root:rows[n]|=1<<v
   else:rows[v]|=1<<n
  if c=='J':root=n
 return rows,root
def literal(rows):
 n=len(rows);full=(1<<n)-1
 dp=[[0]*n for _ in range(full+1)];fixed=[[0]*n for _ in range(full+1)]
 for i in range(n):dp[1<<i][i]=1
 fixed[1][0]=1
 for mask in range(1,full+1):
  for j in range(n):
   if not (mask>>j&1):continue
   rest=mask^(1<<j)
   if rest:
    for i in range(n):
     if rest>>i&1 and rows[i]>>j&1:
      dp[mask][j]+=dp[rest][i];fixed[mask][j]+=fixed[rest][i]
 return sum(dp[full]),sum(fixed[full][j]for j in range(1,n)if rows[j]&1)
def bits(w):
 eps=[int(c=='Z')for c in w];l=len(w)
 A=sum(e*2**i for i,e in enumerate(eps));R=sum(e*2**(l-1-i) for i,e in enumerate(eps))
 C=sum(2**(i-j-1)+(i==j+1)for i in range(l)for j in range(i)if eps[i]and eps[j])
 return A,R,C

def palindrome_inverse(h):
 if h==5:return ('',0)
 if h<9 or h%2==0:return None
 m=(h-2).bit_length() # least m with 2**m >= h-1
 if 3*(2**m)>4*(h-3):return None
 defect=2**m+1-h
 for l in range(m-1):
  if (m-l-2)%2:continue
  for t in itertools.product('OZ',repeat=(l+1)//2):
   v=''.join(t);w=v+v[:l//2][::-1]
   if bits(w)[2]==defect:return w,(m-l-2)//2
 return None

def verify():
 count=0;pal=0;offsets={};literals=[]
 for l in range(13):
  for t in itertools.product('OZ',repeat=l):
   w=''.join(t);A,R,C=bits(w)
   need(coefficients(w)==(2**(l+2),2*(A-R),1-C),'exact bit specialization')
   need((A==R)==(w==w[::-1]),'palindrome iff no linear term')
   count+=1
   if A==R:
    pal+=1;offsets.setdefault(C-1,w)
 for l in range(7):
  for t in itertools.product('OZJ',repeat=l):
   w=''.join(t)
   for s in (1,2,3):
    q,b,c=coefficients(w,s)
    for k in (0,1,2,7):
     need(state(recipe(w,k,s))[0]==q*4**k+b*2**k+c,'general physical quadratic identity')
 # Independent forward image vs the target-first full palindrome inverse.
 forward=set()
 for l in range(11):
  for t in itertools.product('OZ',repeat=(l+1)//2):
   v=''.join(t);w=v+v[:l//2][::-1]
   for k in range(6):
    h=state(recipe(w,k))[0]
    if h<=5001:forward.add(h)
 for h in range(1,5002,2):
  answer=palindrome_inverse(h)
  need((answer is not None)==(h in forward),'complete target-first palindrome inverse')
  if answer is not None:need(state(recipe(*answer))[0]==h,'inverse physical witness')
 fixed=['ZZ','ZOZ','ZZZ','ZOZOZ','ZZOZZ','ZZOJOOZ','ZZJOOOZ','ZZJJOOZ']
 for w in fixed:
  q,b,c=coefficients(w);need(b==0,'selected cancellation')
  for k in range(2):
   word=recipe(w,k);rows,root=compile_graph(word);H,cycles=literal(rows)
   need(H==q*4**k+c and cycles==1,'literal count and unique Hamilton cycle')
   literals.append(dict(middle=w,k=k,count=H,cycles=cycles,rows=rows,root=root))
 for l in range(9):
  for t in itertools.product('OZ',repeat=l):
   v=''.join(t);A,R,C=bits(v);e=int(bool(v) and v[-1]=='Z')
   need(bits(v+v[::-1])[2]==2*C+R*R+e,'even mirror norm')
   need(bits(v+'O'+v[::-1])[2]==2*C+2*R*R,'zero center norm')
   need(bits(v+'Z'+v[::-1])[2]==2*C+2*R*R+2*R+2*e,'one center norm')
 for k in range(1001):
  need(state(recipe('ZZ',k))[0]==2**(2*k+4)-1,'inherited even Mersenne')
  need(state(recipe('ZOZ',k))[0]==2**(2*k+5)-1,'inherited odd Mersenne')
 # Hostile controls: a palindrome need not have every even defect; nonpalindrome has nonzero cross term.
 need(coefficients('OZ')[1]!=0,'asymmetry control')
 need(9 not in offsets,'offset9 is absent in tested palindrome length<=12 only')
 need(state(recipe('ZZZ',0))[0]==27,'non-Mersenne defect5')
 return dict(status='PASS',target_inverse_cap=5001,target_inverse_accepted=len(forward),mirror_norm_words=511,binary_words=count,palindromes=pal,general_words=sum(3**l for l in range(7)),general_substitutions=sum(3**l for l in range(7))*12,literal_graphs=literals,distinct_offsets_length_at_most_12=len(offsets),selected_coefficients={w:coefficients(w)for w in fixed},mersenne_symbolic_evaluations=2002,source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),scope='Universal identities are proved in PROOF.md; finite tests do not prove all-target coverage.')
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--target',type=int);ap.add_argument('--middle');ap.add_argument('--k',type=int,default=0);ap.add_argument('--s',type=int,default=1);ap.add_argument('--fail',action='store_true');a=ap.parse_args()
 if a.fail:need(False,'intentional failing sentinel')
 if a.target is not None:
  answer=palindrome_inverse(a.target)
  if answer is None:print(json.dumps(dict(status='REJECT',target=a.target,scope='synchronized palindrome family only')))
  else:
   w,k=answer;word=recipe(w,k);rows,root=compile_graph(word)
   print(json.dumps(dict(status='ACCEPT',target=a.target,middle=w,k=k,word=word,rows=rows,root=root),indent=2))
 elif a.middle is None:print(json.dumps(verify(),indent=2))
 else:
  w=recipe(a.middle,a.k,a.s);q,b,c=coefficients(a.middle,a.s);rows,root=compile_graph(w)
  print(json.dumps(dict(word=w,coefficients=[q,b,c],count=q*4**a.k+b*2**a.k+c,rows=rows,root=root),indent=2))
