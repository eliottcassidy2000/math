import json,subprocess,itertools,math,hashlib
from pathlib import Path
ROOT=Path(__file__).resolve().parent

def need(c,s):
 if not c:raise RuntimeError(s)
def tt(n):return [sum(1<<j for j in range(i+1,n)) for i in range(n)]
def blow(Q,ns,a,c):
 T=[];off=[]
 for n in ns:off.append(len(T));T.extend([0]*n)
 n=len(T)
 for i,size in enumerate(ns):
  for u in range(size):
   T[off[i]+u]=sum(1<<(off[i]+v) for v in range(u+1,size))|sum(((1<<ns[j])-1)<<off[j] for j in range(7) if Q[i]>>j&1)
   if a>>i&1:T[off[i]+u]|=1<<n
   if c>>i&1:T[off[i]+u]|=1<<(n+1)
 T+=[(1<<(n+1))|sum(((1<<ns[j])-1)<<off[j] for j in range(7) if not a>>j&1),sum(((1<<ns[j])-1)<<off[j] for j in range(7) if not c>>j&1)]
 need(len(T)<=13,'physical maximum')
 return T

def profile(n):
 S=[[0]*6 for _ in range(6)];S[0][0]=1
 for m in range(1,6):
  for k in range(1,m+1):S[m][k]=S[m-1][k-1]+k*S[m-1][k]
 return [math.factorial(k)*S[n][k] for k in range(6)]
models=[x for x in itertools.product(range(1,5),repeat=7) if sum(x)-7<=3]
models += [tuple(5 if i==j else 1 for i in range(7)) for j in range(7)]
models.sort(key=lambda x:(sum(x),x))
need(len(models)==127,'all 127 downclosed models')
R=[{'rows': r} for r in json.load(open(ROOT/'rows.json'))]
qs=json.load(open(ROOT/'even_queries.json'));raw=[list(map(int,x.split())) for x in open(ROOT/'even_bank.txt')]
selections=[(0,1,0),(1,0,3),(1,4,68),(2,0,3),(2,4,68),(3,3,1)]
cases=[];meta=[]
for ri,a,c in selections:
 for side in range(2):
  Q=R[ri]['rows'];Q=Q if side==0 else [127^(1<<i)^r for i,r in enumerate(Q)]
  for ns in models:cases.append(blow(Q,ns,a,c));meta.append((ri,a,c,side,ns))
subprocess.run(['g++','-std=c++17','-O2','-Wall','-Wextra','-Werror',str(ROOT/'physical.cpp'),'-o',str(ROOT/'physical')],check=True)
inp=str(len(cases))+'\n'+''.join(str(len(T))+'\n'+' '.join(map(str,T))+'\n' for T in cases)
out=subprocess.check_output([str(ROOT/'physical')],input=inp.encode()).decode();ls=[list(map(int,x.split())) for x in out.splitlines()]
need(len(ls)==len(meta),'all physical counts')
records=[];coefs=[]
for b,(ri,a,c) in enumerate(selections):
 for side in range(2):
  start=(b*2+side)*len(models);values=ls[start:start+len(models)];K={}
  for ns,vs in zip(models,values):
   v=vs.copy()
   for ks,kv in K.items():
    if all(k<=n for k,n in zip(ks,ns)):
     w=math.prod(profile(n)[k] for k,n in zip(ks,ns))
     for h in range(5):v[h]-=w*kv[h]
   den=math.prod(math.factorial(n) for n in ns);need(all(x%den==0 for x in v),'integral Stirling inversion')
   K[ns]=[x//den for x in v];need(all(x>=0 for x in K[ns]),'nonnegative cut-aware kernel')
  baseline=[0]*5
  for ks,kv in K.items():
   w=math.prod([0,1,14,4,8,0][k] for k in ks)
   for h in range(5):baseline[h]=(baseline[h]+w*kv[h])%16
  correction=[]
  for k in (2,3,5):
   for i in range(7):correction.extend(x%2 for x in K[tuple(k if j==i else 1 for j in range(7))])
  observed=baseline+correction;qi=qs.index([a,c]);expected=raw[2*(ri*len(qs)+qi)+side];need(observed==expected,'independent full-cover physical interpolation')
  records.append(dict(quotient_index=ri,query=[a,c],side=side,values=values,interpolated=[list(k)+v for k,v in K.items()],baseline_mod16=baseline,correction_parities=correction))
result=dict(status='PASS',scope='Independent physical subset DP and exact triangular Stirling inversion reproduce every baseline and every q2/q3/q5 profile derivative for the six selected queries, both orientations, first five ordered-cover coordinates.',models=[list(x) for x in models],records=records,physical_cases=len(cases),physical_order_max=max(map(len,cases)))
print(json.dumps(result,sort_keys=True,indent=2))
