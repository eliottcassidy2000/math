#!/usr/bin/env python3
"""Independent direct path-DP checks; no imports from construction lanes."""
from itertools import combinations, permutations
from functools import lru_cache
import json
checks=0

def require(x,label):
 global checks
 checks+=1
 if not x: raise RuntimeError(label)

def tt(n):return tuple(sum(1<<j for j in range(i+1,n)) for i in range(n))
def all_t(n):
 pairs=list(combinations(range(n),2))
 for k in range(1<<len(pairs)):
  r=[0]*n
  for p,(a,b) in enumerate(pairs):
   if k>>p&1:r[a]|=1<<b
   else:r[b]|=1<<a
  yield tuple(r)
def join(*blocks):
 n=sum(map(len,blocks));out=[];off=0
 for r in blocks:
  for row in r:out.append((row<<off)|sum(1<<j for j in range(off+len(r),n)))
  off+=len(r)
 return tuple(out)
@lru_cache(None)
def endpoints(r,mask):
 n=len(r);out=[[0]*n for _ in range(n)]
 for start in range(n):
  if not mask>>start&1:continue
  d={(1<<start,start):1}
  for size in range(1,mask.bit_count()):
   nxt={}
   for (used,last),num in d.items():
    choices=r[last]&mask&~used
    while choices:
     bit=choices&-choices;choices-=bit;j=bit.bit_length()-1
     key=(used|bit,j);nxt[key]=nxt.get(key,0)+num
   d=nxt
  for (used,last),num in d.items():out[start][last]+=num
 return tuple(map(tuple,out))
def H(r,mask=None):
 if mask is None:mask=(1<<len(r))-1
 return 1 if mask==0 else sum(map(sum,endpoints(r,mask)))
def extend(r,cut):
 n=len(r);return tuple(row|((0 if cut>>i&1 else 1)<<n) for i,row in enumerate(r))+(cut,)
def response(r,a,start):
 n=len(r);full=(1<<n)-1;ans=0
 for mask in range(1,full+1):
  if mask>>a&1:
   e=endpoints(r,mask)
   val=sum(e[a]) if start else sum(row[a] for row in e)
   ans+=val*H(r,full^mask)
 return ans

def apex(word):
 n=len(word);r=list(tt(n));r.append(0)
 for i,b in enumerate(word):
  if b=='1':r[n]|=1<<i
  else:r[i]|=1<<n
 return tuple(r)
def formula(word):
 return 1+2*sum((1 if j-i<=2 else 2**(j-i-2)) for i in range(len(word)) for j in range(i+1,len(word)) if word[i]=='1' and word[j]=='0')

words=direct=0;zero_counts={}
for m in range(2,15):
 zero=[];plus=[];minus=[]
 for val in range(1<<m):
  w=f'{val:0{m}b}';X='10'+w+'01';Y='01'+w+'10'
  B=int(w[:-1],2);C=int(''.join('1' if z=='0' else '0' for z in w[:0:-1]),2)
  require(formula(Y)-formula(X)==2*(B-C),'two-switch full invoice')
  require((B==C)==all(w[i]!=w[m-1-i] for i in range(m)),'zero iff anti-palindrome')
  if B==C:
   zero.append(w)
   require(''.join('1' if z=='0' else '0' for z in X[:-1][::-1])==Y[1:],'zero converse cores')
  if B-C==1:plus.append(w)
  if B-C==-1:minus.append(w)
  if m<=5:
   require(H(apex(X))==formula(X) and H(apex(Y))==formula(Y),'independent direct full H')
   direct+=1
  words+=1
 expected=['11'] if m==2 else ['011','101','110'] if m==3 else ['1'+'0'*(m-2)+'1']
 require(plus==expected,'complete positive carry language')
 require(set(minus)=={''.join('1' if b=='0' else '0' for b in w) for w in expected},'negative carry language')
 require(len(zero)==(2**(m//2) if m%2==0 else 0),'zero language size')
 zero_counts[m]=len(zero)
 if m>=3:
  w='1'+'0'*(m-2)+'1'
  require(formula('10'+w+'01')==15*2**(m-2)+1,'dyadic X')
  require(formula('01'+w+'10')==15*2**(m-2)+3,'dyadic Y')

# Both sides arbitrary through order3; feedback through m=0,1,2 seams.
activation=0
blocks=[t for n in range(1,4) for t in all_t(n)]
for A in blocks:
 for B in blocks:
  for a in range(len(A)):
   for b in range(len(B)):
    L=response(A,a,True);R=response(B,b,False)
    for m in range(3):
     T=list(join(A,tt(m),B));bb=len(A)+m+b
     T[a]^=1<<bb;T[bb]^=1<<a;T=tuple(T)
     correction=0
     if m==0:
      EA=endpoints(A,(1<<len(A))-1);EB=endpoints(B,(1<<len(B))-1)
      correction=sum(row[a] for row in EA)*sum(EB[b])
     require(H(T)==H(A)*H(B)+(2**m)*L*R-correction,'seam activation including empty boundary')
     activation+=1
# Preserve ALL old endpoint entries, not just starts/ends separately.
A=(2,4,1);B=tt(2);base=join(A,B);E=endpoints(base,31)
for m in range(1,6):
 T=join(A,tt(m),B);N=len(T);F=endpoints(T,(1<<N)-1)
 old=[0,1,2,3+m,4+m]
 for i in range(5):
  for j in range(5):require(F[old[i]][old[j]]==E[i][j],'old endpoint matrix invariance')
 for s in range(3,3+m):require(sum(F[s])==0 and sum(row[s] for row in F)==0,'seam global endpoint invisibility')
 require(H(T)==3,'neutral chain full H')
# Sector-renewal hostile: find an endpoint column with multiple starts in B4.
Bnext=None
for B4 in all_t(4):
 E4=endpoints(B4,15)
 for t in range(4):
  e=sum(row[t] for row in E4)
  for b in range(4):
   if b!=t and 0<E4[b][t]<e:
    Bnext=(B4,b,t,E4[b][t],e);break
  if Bnext:break
 if Bnext:break
require(Bnext is not None,'column guard hostile exists')
# Direct permutations count selected predecessor/successor sectors.
B4,b,t,hbt,et=Bnext
base=join(A,B4)
def sector(r,start,end,z,a,b):
 return sum(1 for p in permutations(range(len(r))) if p[0]==start and p[-1]==end and all(r[x]>>y&1 for x,y in zip(p,p[1:])) and 0<p.index(z)<len(p)-1 and p[p.index(z)-1]==a and p[p.index(z)+1]==b)
T1=join(A,tt(1),B4);T2=join(A,tt(2),B4)
# a=2->v=0 in C3. Newly added newest seam is index3; former seam index4.
s1=sector(T1,0,4+t,3,2,4+b);s2=sector(T2,0,5+t,3,2,4)
require((s1,s2)==(hbt,et),'initial sector differs after unguarded port renewal')
print(json.dumps({'status':'PASS','checks':checks,'words_m2_to14':words,'direct_H_word_pairs':direct,'activation_cases':activation,'zero_counts':zero_counts,'sector_guard_hostile':{'B_rows':B4,'b':b,'t':t,'first':s1,'renewed':s2}},sort_keys=True,indent=2))
