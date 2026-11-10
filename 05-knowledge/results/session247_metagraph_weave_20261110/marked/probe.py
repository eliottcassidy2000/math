#!/usr/bin/env python3
"""Small probes of the NEW marked formula, not an inherited-result audit."""
from itertools import combinations,product
from math import factorial,gcd
from functools import reduce,lru_cache
import json

def require(v,msg):
 if not v: raise RuntimeError(msg)

def adj(n,code):
 a=[[0]*n for _ in range(n)]
 for k,(i,j) in enumerate(combinations(range(n),2)):
  a[i][j]=(code>>k)&1; a[j][i]=1-a[i][j]
 return a

def tables(a):
 n=len(a);N=1<<n;dp=[[0]*n for _ in range(N)];sp=[[0]*n for _ in range(N)];H=[0]*N;H[0]=1
 for i in range(n):dp[1<<i][i]=sp[1<<i][i]=1
 for m in range(1,N):
  for j in range(n):
   if m>>j&1:
    rest=m^(1<<j)
    if rest:
     dp[m][j]=sum(dp[rest][i] for i in range(n) if a[i][j])
     sp[m][j]=sum(sp[rest][i] for i in range(n) if a[j][i])
  H[m]=sum(dp[m])
 return H,dp,sp

@lru_cache(None)
def state(code,cut):
 H,e,s=tables(adj(6,code));N=63
 E=[sum(e[m][v] for v in range(6) if not(cut>>v&1)) for m in range(64)]
 S=[sum(s[m][v] for v in range(6) if cut>>v&1) for m in range(64)]
 return dict(h=H[N],Q=sum(H[m]*H[N^m] for m in range(64)),E=E[N],S=S[N],
   D=sum(E[m]*S[N^m] for m in range(1,63)),
   A=sum(E[m]*H[N^m] for m in range(1,64)),B=sum(S[m]*H[N^m] for m in range(1,64)))

def formula(ss):
 k=len(ss)
 def hp(ids):
  v=1
  for i in ids:v*=ss[i]['h']
  return v
 out=ss[0]['S']*hp(range(1,k))+ss[-1]['E']*hp(range(k-1))
 for i in range(k):out+=ss[i]['D']*hp(j for j in range(k) if j!=i)
 for i in range(k-1):out+=ss[i]['E']*ss[i+1]['S']*hp(j for j in range(k) if j not in (i,i+1))
 for r in range(k):
  for l in range(r+1,k):
   v=ss[r]['B']*ss[l]['A']*hp(j for j in range(k) if j<r or j>l)
   for j in range(r+1,l):v*=ss[j]['Q']
   out+=v
 return out

def physical(codes,cuts):
 k=len(codes);n=6*k+1;x=n-1;a=[[0]*n for _ in range(n)]
 for u in range(k):
  core=adj(6,codes[u])
  for i in range(6):
   for j in range(6):a[6*u+i][6*u+j]=core[i][j]
   a[x][6*u+i]=(cuts[u]>>i)&1;a[6*u+i][x]=1-a[x][6*u+i]
 for u in range(k):
  for v in range(u+1,k):
   for i in range(6):
    for j in range(6):a[6*u+i][6*v+j]=1
 return a

def count(a):
 n=len(a);N=1<<n;dp=[[0]*n for _ in range(N)]
 pred=[sum(1<<i for i in range(n) if a[i][j]) for j in range(n)]
 for i in range(n):dp[1<<i][i]=1
 for m in range(1,N):
  bits=m
  while bits:
   jbit=bits&-bits;bits-=jbit;j=jbit.bit_length()-1;prev=m^jbit
   if prev:
    b=prev&pred[j];s=0
    while b:
     ibit=b&-b;b-=ibit;s+=dp[prev][ibit.bit_length()-1]
    dp[m][j]=s
 return sum(dp[-1])

def main():
 receipts=[]
 for k in (1,2):
  vals={}
  for bits in product((0,1),repeat=k):
   codes=[18+4*b for b in bits];ss=[state(c,5) for c in codes]
   theory=formula(ss);literal=count(physical(codes,[5]*k));require(theory==literal,'one-root formula')
   vals[''.join(map(str,bits))]=literal
  receipts.append({'k':k,'literal_counts':vals})
 require(receipts[0]['literal_counts']=={'0':41,'1':43},'seed')
 require(receipts[1]['literal_counts']=={'00':5605,'01':5635,'10':6005,'11':6035},'cube')
 # A nondivisible-by-3 two-module example shows the k-2 exponent cannot
 # be raised uniformly to k-1 for the one-root construction.
 val=count(physical([18,18],[1,2]));require(val==formula([state(18,1),state(18,2)])==3395,'boundary example')
 symbolic=[]
 for k in range(1,9):
  base=formula([state(18,5)]*k);gs=[]
  for i in range(k):
   ss=[state(18,5)]*k;ss[i]=state(22,5);g=formula(ss)-base;gs.append(g)
   exact=2*15**(k-1) if i==k-1 else 43*15**(k-2)+357*sum(15**(k-d-1)*216**(d-1) for d in range(1,k-i))
   require(g==exact,'closed coefficients')
  for bits in product((0,1),repeat=k):
   f=formula([state(18+4*b,5) for b in bits]);require(f==base+sum(g*b for g,b in zip(gs,bits)),'affine response')
   require(f%(3**max(0,k-2))==0,'interface divisibility')
  symbolic.append({'k':k,'baseline':base,'increments':gs})
 result={'scope':'Literal subset-DP only for k=1,2 and one boundary example. Larger rows evaluate the proved formula, not independent tournaments.',
   'inherited_seed_profiles_not_reaudited':True,'states':{str(c):state(c,5) for c in (18,22)},
   'literal':receipts,'sharp_zero_exponent_example':{'root_masks':[1,2],'H':val},'formula_rows':symbolic,
   'ordered_profile_gcd':reduce(gcd,[factorial(c)*v for c,v in enumerate([15,93,138,73,15,1],1)])}
 print(json.dumps(result,indent=2))
if __name__=='__main__':main()
