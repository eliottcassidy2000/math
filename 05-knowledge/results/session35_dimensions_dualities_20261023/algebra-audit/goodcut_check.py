"""Direct interval-subset census versus the printed/corrected THM341 GFs."""
from collections import Counter
from itertools import combinations
from math import comb, factorial
import json

def require(c,label):
 if not c:raise RuntimeError(label)
def direct(N):
 intervals=[sum(1<<j for j in range(a,b)) for a in range(N) for b in range(a+2,N+1)]
 unions=[0]
 for interval in intervals:unions += [u|interval for u in unions]
 return Counter(u.bit_count() for u in unions),Counter(unions)
cover={}
for L in range(2,7):
 _,unions=direct(L); cover[L]=unions[(1<<L)-1]
require([cover[k] for k in range(2,6)]==[1,5,50,903],'connected covers')
maxn=7
printed=[[0]*(maxn+1) for _ in range(maxn+1)]
for n in range(maxn+1):
 for k in range(maxn+1):
  printed[n][k]=(printed[n-1][k] if n else 0)
  for L,c in cover.items():
   if n>=L+1 and k>=L:
    printed[n][k]+=c*printed[n-L-1][k-L]
    if n==L+1 and k==L:printed[n][k]+=c
corrected=[[sum(printed[j][k] for j in range(n+1)) for k in range(maxn+1)] for n in range(maxn+1)]
census=[]
for n in range(1,maxn+1):
 bins,_=direct(n-1)
 for k in range(2,n):require(corrected[n][k]==bins[k],('corrected GF',n,k))
 require(sum(bins.values())==2**comb(n-1,2),'all configurations')
 census.append({'n':n,'direct_buckets':dict(sorted(bins.items()))})
require(printed[4][2]==1 and corrected[4][2]==2,'canon hostile')
# Independent maximal-run coefficient and leading-order checks by ordered compositions.
def compositions(k,r):
 if r==0:
  if k==0:yield ()
  return
 for a in range(2,k-2*(r-1)+1):
  for tail in compositions(k-a,r-1):yield (a,)+tail
lead=[]
for k in range(2,13):
 r=k//2;profiles=list(compositions(k,r))
 require(all(all(a==2 for a in p) for p in profiles) if k%2==0 else all(p.count(3)==1 and p.count(2)==r-1 for p in profiles),'maximal morphology')
 weight=sum(5**p.count(3) for p in profiles)
 require(weight==(1 if k%2==0 else 5*r),'leading weight')
 lead.append({'k':k,'degree':r,'leading_numerator':weight,'leading_denominator':factorial(r),'maximal_profiles':len(profiles)})
print(json.dumps({'connected_cover_counts':cover,'canon_hostile':{'n':4,'good_cuts':2,'printed':printed[4][2],'corrected_and_actual':corrected[4][2]},'direct_census':census,'fixed_k_leading':lead},indent=2))
