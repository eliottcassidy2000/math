#!/usr/bin/env python3
import itertools as it,json,math
from common import *
def skeleton(r):
 counts={}
 def walk(last,left,i,j):
  if not left and i and j:counts[i,j]=counts.get((i,j),0)+1
  for x in ['A','L',*sorted(left)]:
   if x=='A' and i==r+1 or x=='L' and j==r+1:continue
   ok=last is None or (last=='A' and x=='L') or (last=='L' and isinstance(x,int)) or (isinstance(last,int) and x=='A') or (isinstance(last,int) and isinstance(x,int) and last<x)
   if ok:walk(x,left-{x} if isinstance(x,int) else left,i+(x=='A'),j+(x=='L'))
 walk(None,set(range(r)),0,0)
 return counts
def kernel(C,r):
 K=[[0]*(r+1) for _ in range(r+1)]
 for (i,j),c in C.items():
  for s in range(1,i+1):
   for t in range(1,j+1):K[s-1][t-1]+=c*(-1)**(i-s+j-t)*math.comb(i,s)*math.comb(j,t)
 return K
require(all(strong(b) for b in B),'strong seeds')
require(all(any(((B[0][i]>>j)&1)!=((B[1][perm[i]]>>perm[j])&1) for i in range(5) for j in range(5)) for perm in it.permutations(range(5))),'nonisomorphic seeds')
require(cover(B[0])==cover(B[1]),'equal seed profiles')
require(all(evalF(cover(B[0]),j)==f(j) for j in range(7)),'seed polynomial')
records=[];kernels={}
for r in range(1,5):
 K=kernel(skeleton(r),r);kernels[r]=K
 require(all(K[i][j]==K[j][i] for i in range(r+1) for j in range(r+1)),'symmetric kernel')
 for a in [(0,),(2,0),(2,4,1),*B]:
  g=wrap(a,r);n=len(a);pi=list(range(n,2*n))+list(range(n))+list(reversed(range(2*n,2*n+r)))
  require(all(((g[i]>>j)&1)==((g[pi[j]]>>pi[i])&1) for i in range(len(g)) for j in range(len(g))),'explicit anti-automorphism')
  p=cover(a);F=[evalF(p,j) for j in range(1,r+2)]
  expected=sum(K[i][j]*F[i]*F[j] for i in range(r+1) for j in range(r+1));actual=paths(g)[-1]
  require(expected==actual,'literal Hamilton count')
  records.append({'r':r,'module':list(a),'order':len(g),'H':actual})
require(kernels[1]==[[3,-1],[-1,1]],'r1 kernel')
require(kernels[2]==[[15,-11,2],[-11,13,-4],[2,-4,2]],'r2 kernel')
profiles={}
for r in (1,2):
 profiles[r]=[cover(wrap(b,r)) for b in B]
 require(profiles[r][0]==profiles[r][1],'profile-blind mirrored seeds')
for r in range(1,5):
 a=wrap(B[0],r)
 for k in range(r):require(all(survive(a,set(S)) for S in it.combinations(range(len(a)),k)),'resilience')
 cuts=[S for S in it.combinations(range(len(a)),r) if not survive(a,set(S))]
 require(cuts==[tuple(range(10,10+r))],'unique cut')
clocks=[]
for M in range(1,51):
 h1=3*81**M-2*756**M+7056**M
 h2=15*81**M-22*756**M+4*3699**M+13*7056**M-8*34524**M+2*168921**M
 for h in (h1,h2):require(h%3**(2*M)==0 and h%3**(2*M+1)!=0,'exact valuation')
 require(h1%10==7 and h1%8==3,'r1 residue clocks')
 require(h2%10==(7 if M%2 else 9),'r2 last-digit clock')
 clocks.append({'M':M,'H1_mod10':h1%10,'H2_mod10':h2%10,'v3':2*M})
pal_tests=0
for m in range(2,41):
 for center in range(10):
  for digits in it.product((1,3,9),repeat=3):
   p=center;a=100%m;s=str(center)
   for d in digits:
    p=(10*p+d*(a+1))%m;a=100*a%m;s=str(d)+s+str(d)
    require(p==int(s)%m,'palindrome finite-state encoder');pal_tests+=1
for k in range(1,101):
 p=10**(2*k+1)-10**k-1
 require(str(p)=='9'*k+'8'+'9'*k,'large palindrome identity')
 require(p%3==2 and p%11==(-2-(-1)**k)%11,'central-eight clocks')
require(1001==7*11*13 and 2222==2*11*101,'even palindrome examples')
print(json.dumps({'kernels':kernels,'literal_counts':records,'full_profiles':profiles,'clocks':clocks,'palindrome_transition_checks':pal_tests,'large_expression_digits':2*944264+1},sort_keys=True,indent=2))
