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
def state_n(n,code,cut):
 H,e,s=tables(adj(n,code));N=(1<<n)-1
 E=[sum(e[m][v] for v in range(n) if not(cut>>v&1)) for m in range(N+1)]
 S=[sum(s[m][v] for v in range(n) if cut>>v&1) for m in range(N+1)]
 return dict(h=H[N],Q=sum(H[m]*H[N^m] for m in range(N+1)),E=E[N],S=S[N],
   D=sum(E[m]*S[N^m] for m in range(1,N)),
   A=sum(E[m]*H[N^m] for m in range(1,N+1)),B=sum(S[m]*H[N^m] for m in range(1,N+1)))

def state(code,cut): return state_n(6,code,cut)

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

def physical_n(m,codes,cuts):
 k=len(codes);n=m*k+1;x=n-1;a=[[0]*n for _ in range(n)]
 for u in range(k):
  core=adj(m,codes[u])
  for i in range(m):
   for j in range(m):a[m*u+i][m*u+j]=core[i][j]
   a[x][m*u+i]=(cuts[u]>>i)&1;a[m*u+i][x]=1-a[x][m*u+i]
 for u in range(k):
  for v in range(u+1,k):
   for i in range(m):
    for j in range(m):a[m*u+i][m*v+j]=1
 return a

def physical(codes,cuts): return physical_n(6,codes,cuts)

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

def profile(n,code):
 H,_,_=tables(adj(n,code))
 @lru_cache(None)
 def rec(mask):
  if not mask:return (1,)+(0,)*n
  first=mask&-mask;out=[0]*(n+1);sub=mask
  while sub:
   if sub&first:
    prev=rec(mask^sub)
    for c in range(1,n+1):out[c]+=H[sub]*prev[c-1]
   sub=(sub-1)&mask
  return tuple(out)
 return rec((1<<n)-1)[1:]

def chargefree():
 profiles=[profile(7,c) for c in (568,696)]
 require(profiles[0]==profiles[1]==(101,585,915,570,162,21,1),'new charge-free profile')
 vals={}
 for bits in product((0,1),repeat=2):
  codes=[568+128*b for b in bits];ss=[state_n(7,c,53) for c in codes]
  got=count(physical_n(7,codes,[53,53]));require(got==formula(ss),'charge-free literal formula')
  vals[''.join(map(str,bits))]=got
 require(vals=={'00':410445,'01':411051,'10':412729,'11':413335},'charge-free cube')
 rows=[]
 for k in range(2,9):
  base=formula([state_n(7,568,53)]*k);gs=[]
  for i in range(k):
   ss=[state_n(7,568,53)]*k;ss[i]=state_n(7,696,53);g=formula(ss)-base
   exact=6*101**(k-1) if i==k-1 else 661*101**(k-2)+1623*sum(101**(k-d-1)*1372**(d-1) for d in range(1,k-i))
   require(g==exact,'charge-free increments');gs.append(g)
  require(reduce(gcd,[g//2 for g in gs])==1,'primitive increments')
  for bits in product((0,1),repeat=k):
   require(formula([state_n(7,568+128*b,53) for b in bits])==base+sum(g*b for g,b in zip(gs,bits)),'charge-free affine')
  rows.append({'k':k,'baseline':base,'increments':gs,'normalized_gcd':1})
 return {'seed_profiles':profiles,'states':{str(c):state_n(7,c,53) for c in (568,696)},'literal_k2':vals,'formula_rows':rows}

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
 # Sharp degree-two physical response of the invisible pair.
 sharp={}
 for bits in product((0,1),repeat=2):
  codes=[18+4*b for b in bits]
  got=count(physical(codes,[1,1]));require(got==formula([state(c,1) for c in codes]),'degree-two formula')
  sharp[''.join(map(str,bits))]=got
 require(sharp=={'00':4161,'01':4297,'10':4933,'11':5085},'degree-two literal values')
 require(sharp['11']-sharp['10']-sharp['01']+sharp['00']==16,'degree-two sharpness')
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
 result={'degree_two_sharpness':sharp,'charge_free_repair':chargefree(),'scope':'Literal subset-DP for the stated k=1,2 six-vertex cases, charge-free k=2, and one boundary example. Larger rows evaluate the proved formula, not independent tournaments.',
   'inherited_seed_profiles_not_reaudited':True,'states':{str(c):state(c,5) for c in (18,22)},
   'literal':receipts,'sharp_zero_exponent_example':{'root_masks':[1,2],'H':val},'formula_rows':symbolic,
   'ordered_profile_gcd':reduce(gcd,[factorial(c)*v for c,v in enumerate([15,93,138,73,15,1],1)])}
 print(json.dumps(result,indent=2))
if __name__=='__main__':main()
