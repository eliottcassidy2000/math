#!/usr/bin/env python3
import math,json
from pathlib import Path

def require(x,msg):
 if not x:raise RuntimeError(msg)
N=20
T=[1]+[2**((n-1)*(n-2)//2) for n in range(1,N+1)]
s=[0]*(N+1)
for n in range(1,N+1):s[n]=T[n]-sum(s[j]*T[n-j] for j in range(1,n))
# Ordered run-weight convolutions Q[k][r].
Q=[[0]*(N+1) for _ in range(N+1)];Q[0][0]=1
for k in range(2,N+1):
 for r in range(1,N+1):Q[k][r]=sum(s[l+1]*Q[k-l][r-1] for l in range(2,min(k,N-1)+1))
def formula(n,k):
 if k==0:return 1
 return sum(Q[k][r]*math.comb(n-k,r) for r in range(1,k//2+1) if n-k>=r)
# Independent interval-union enumeration via union DP, then graph reachability.
results={}
for n in range(1,9):
 tiles=[(i,j) for i in range(n) for j in range(i+2,n)]
 intervals=[((1<<(j-i))-1)<<i for i,j in tiles]
 size=1<<len(tiles);unions=[0]*size;bucket=[0]*n
 for mask in range(size):
  if mask:
   b=mask&-mask;unions[mask]=unions[mask^b]|intervals[b.bit_length()-1]
  k=unions[mask].bit_count();bucket[k]+=1
  if n<=7:
   reach=[sum(1<<j for j in range(i+1,n))|(1<<i) for i in range(n)]
   for t,(i,j) in enumerate(tiles):
    if mask>>t&1:reach[i]^=1<<j;reach[j]|=1<<i
   for z in range(n):
    for i in range(n):
     if reach[i]>>z&1:reach[i]|=reach[z]
   reps=[]
   for i in range(n):
    if not any((reach[i]>>j&1) and (reach[j]>>i&1) for j in reps):reps.append(i)
   require(k==n-len(reps),'good count = SCC defect')
 require(bucket==[formula(n,k) for k in range(n)],'brute/formula agreement')
 require(bucket[-1]==s[n] if n>1 else True,'strong connected count')
 results[n]=bucket
for n in range(1,N+1):require(sum(formula(n,k) for k in range(n))==T[n],'all counts sum to total')
# Correct vs erroneous vertex-count positive-good GF, integer coefficients.
# Denominator1-x-xA(xu): recurrence G[n,k]=G[n-1,k]+sum_l s[l+1]G[n-l-1,k-l].
G=[[0]*(N+1) for _ in range(N+1)];G[0][0]=1
for n in range(1,N+1):
 for k in range(N+1):
  G[n][k]=G[n-1][k]+sum(s[l+1]*G[n-l-1][k-l] for l in range(2,min(k,n-1)+1))
wrong=[[0]*(N+1) for _ in range(N+1)]
for n in range(1,N+1):
 for k in range(2,n):wrong[n][k]=sum(s[l+1]*G[n-l-1][k-l] for l in range(2,min(k,n-1)+1))
for n in range(1,N+1):
 for k in range(2,n):require(sum(wrong[j][k] for j in range(n+1))==formula(n,k),'missing geometric factor repair')
require(wrong[4][2]==1 and formula(4,2)==2,'minimal GF hostile')
require(s[6]==903,'reject874')
out={'status':'FINITE-EXACT','strong_fixed_path_counts_n1_to10':s[1:11],'brute_buckets_n1_to8':results,'SCC_checks_through_n':7,'formula_total_checks_through_n':20,'wrong_GF_first_hostile':{'n':4,'good_cuts':2,'printed':wrong[4][2],'correct':formula(4,2)},'Q_k6':Q[6][1:4],'Q_k7':Q[7][1:4],'F_n10_k6':formula(10,6)}
Path(__file__).with_name('goodcuts_receipt.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
