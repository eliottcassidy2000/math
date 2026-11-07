#!/usr/bin/env python3
"""Deterministic exact fixed-cap checks; no external dependencies."""
from itertools import permutations
from math import comb,factorial
import json
from fractions import Fraction
from probe import tournament,profile,q

def check(x,msg):
 if not x:raise RuntimeError(msg)
def transitive(n):return [[int(i<j) for j in range(n)] for i in range(n)]
def join(blocks,cyclic=False):
 sizes=list(map(len,blocks));n=sum(sizes);a=[[0]*n for _ in range(n)];off=[0]
 for z in sizes:off.append(off[-1]+z)
 for i,B in enumerate(blocks):
  for u in range(sizes[i]):
   for v in range(sizes[i]):a[off[i]+u][off[i]+v]=B[u][v]
 for i in range(len(blocks)):
  for j in range(i+1,len(blocks)):
   edge=not(cyclic and i==0 and j==2)
   for u in range(off[i],off[i+1]):
    for v in range(off[j],off[j+1]):a[u][v]=int(edge);a[v][u]=int(not edge)
 return a

def hp(a):
 n=len(a);N=1<<n;dp=[[0]*n for _ in range(N)]
 for v in range(n):dp[1<<v][v]=1
 out=[sum(1<<w for w in range(n) if a[v][w]) for v in range(n)]
 for s in range(1,N):
  for v,x in enumerate(dp[s]):
   if not x:continue
   todo=out[v]&~s
   while todo:
    bit=todo&-todo;todo-=bit;dp[s|bit][bit.bit_length()-1]+=x
 return sum(dp[-1])
def perm_profile(a):
 n=len(a);f=[0]*(n+1)
 for p in permutations(range(n)):
  b=sum(not a[p[i]][p[i+1]] for i in range(n-1));g=n-1-b
  for k in range(b+1,n+1):f[k]+=comb(g,k-b-1)
 return f

def coeff(A,C):
 fa,fc=profile(A),profile(C);m=min(len(A),len(C));c=[0]*(m+2)
 for x in range(1,len(fa)):
  for y in range(1,m+2):
   for z in range(1,len(fc)):
    if max(x,y,z)-min(x,y,z)<=1:c[y]+=(3 if x==y==z else 1)*fa[x]*fc[z]
 return c
def beta(c):return [0]+[sum(c[k]*(-1)**(k-j)*comb(k,j) for k in range(j,len(c))) for j in range(1,len(c))]
def predicted(A,B,C):
 c=coeff(A,C);f=profile(B)
 return sum(c[k]*f[k] for k in range(1,min(len(c),len(f))))
def recover(n,obs):
 f=[0,obs[0]]
 for m in range(1,n):
  c=coeff(transitive(m),transitive(m));num=obs[m]-sum(c[k]*f[k] for k in range(1,m+1))
  check(num%factorial(m)**2==0,'integrality');f.append(num//factorial(m)**2)
 return f

def short_recover(n,obs):
 d=(n-1)//2
 if d==0:return profile(transitive(n))
 f=recover(d,obs)
 # Unknown powers have order parity; leading coefficient is 1.
 powers=list(range(1 if n%2 else 2,n,2))
 mat=[[Fraction(j**e) for e in powers]+[Fraction(q(f,j)-j**n)] for j in range(1,d+1)]
 for i in range(d):
  pivot=next(r for r in range(i,d) if mat[r][i]);mat[i],mat[pivot]=mat[pivot],mat[i]
  scale=mat[i][i];mat[i]=[v/scale for v in mat[i]]
  for r in range(d):
   if r!=i:
    scale=mat[r][i];mat[r]=[v-scale*w for v,w in zip(mat[r],mat[i])]
 coef=[row[-1] for row in mat]
 vals=[Fraction(j**n)+sum(c*j**e for c,e in zip(coef,powers)) for j in range(n+1)]
 out=[sum((-1)**(k-j)*comb(k,j)*vals[j] for j in range(k+1)) for k in range(n+1)]
 check(all(x.denominator==1 for x in out),'short recovery integral')
 return list(map(int,out))

def main():
 receipt={};pc_checks=0;physical=0;tomography=0;short_checks=0
 cases=[]
 for n in range(1,5):
  for b in range(1<<(n*(n-1)//2)):cases.append(tournament(n,b))
 for B in cases:
  f=profile(B);check(f==perm_profile(B),'permutation cover');pc_checks+=1
  for m in (1,2):
   A=C=transitive(m);v=hp(join([A,B,C],True));check(v==predicted(A,B,C),'physical cap');physical+=1
  obs=[f[1]]+[predicted(transitive(m),B,transitive(m)) for m in range(1,len(B))]
  check(recover(len(B),obs)==f,'tomography');tomography+=1
 # Deep tomography uses independently computed physical path counts.
 D=tournament(5,12);f=profile(D);check(f==perm_profile(D),'seed profile')
 obs=[hp(D)]
 for m in range(1,5):obs.append(hp(join([transitive(m),D,transitive(m)],True)));physical+=1
 check(recover(5,obs)==f,'physical tomography');tomography+=1
 for B in cases+[D,tournament(6,1729),tournament(7,19231),tournament(8,177777)]:
  n=len(B);d=(n-1)//2;fb=profile(B);ob=[fb[1]]+[predicted(transitive(m),B,transitive(m)) for m in range(1,d)]
  check(short_recover(n,ob)==fb,'parity short recovery');short_checks+=1
 receipt['short_recovery_checks']=short_checks
 bridge=beta(coeff(tournament(3,2),transitive(3)))
 check(bridge==[0,87,-189,72,36],'descent compiler bridge')
 receipt['multi_exit_bridge_beta']=bridge
 receipt['seed']={'order':5,'upper_triangle_mask':12,'adjacency_rows':[''.join(map(str,r)) for r in D],'F':f,'Q2':q(f,2),'physical_cap_observations':obs}
 # Asymmetric caps and mixed word; beta can be negative.
 words=[[D],[D,D],[transitive(1),tournament(3,2),D]]
 for A,C in [(transitive(1),transitive(1)),(transitive(2),tournament(3,2))]:
  c=coeff(A,C);b=beta(c)
  for word in words:
   B=join(word);value=sum(b[j]*__import__('math').prod(q(profile(E),j) for E in word) for j in range(1,len(b)))
   actual=hp(join([A,B,C],True));check(value==actual,'mixed-word eigenmodes');physical+=1
 receipt['pump']={'formula':'100^r + 11^r','values_r1_to8':[100**r+11**r for r in range(1,9)],'order_formula':'5r+2','recurrence':'u_(r+2)=111*u_(r+1)-1100*u_r','physical_r':[1,2]}
 # The sign error 4h-Q2 instead of 3h-Q2 is a live hostile.
 c=coeff(transitive(2),transitive(1));b=beta(c);check(b==[0,-1,3],'asymmetric signed coefficients')
 same=[tournament(5,b) for b in (69,76)]
 sf=[profile(B) for B in same];sv=[hp(join([transitive(1),B,transitive(1)],True)) for B in same]
 check([f[1] for f in sf]==[15,15] and sv==[123,135],'inherited equal-order hostile')
 receipt['inherited_same_order_hostile']={'masks':[69,76],'F':sf,'capped_H':sv,'credit':'Session176 information PROOF.md section 5, relabelled witnesses'}
 receipt['hostiles']={'TT1_and_TT2_same_H':1,'singleton_cap_values':[hp(join([transitive(1),transitive(n),transitive(1)],True)) for n in (1,2)],'TT2_TT1_cap_beta':b,'no_empty_block':'All blocks nonempty; r >= 1. Recurrence may algebraically extend to r=0 but no empty C3 block is asserted.'}
 receipt['checks']={'permutation_vs_subset_profiles':pc_checks+1,'physical_H_compositions':physical,'profile_recovery':tomography,'all_labelled_orders':[1,2,3,4],'large_seed_cases':'mask12 order5, two copies order10, cap order up to4','status':'PASS'}
 print(json.dumps(receipt,indent=2,sort_keys=True))
if __name__=='__main__':main()
