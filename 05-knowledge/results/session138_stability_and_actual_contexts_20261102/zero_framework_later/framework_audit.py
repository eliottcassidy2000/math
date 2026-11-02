"""Actual fixed-module and independent literal context controls."""
import json,math,itertools,hashlib
from pathlib import Path
from fractions import Fraction
from certify import certify_actual,cyclic_kernel,quadratic,negative_certificate,validate_tournament
from profile_tools import covers
P=Path(__file__).resolve().parent

def literal_profile(a):
 n=len(a);f=[0]*n
 for p in itertools.permutations(range(n)):
  d=sum(not(a[p[i]]>>p[i+1]&1) for i in range(n-1))
  for r in range(d+1,n+1):f[r-1]+=math.comb(n-1-d,r-1-d)
 return f

def h_forward(a):
 n=len(a);dp=[[0]*n for _ in range(1<<n)]
 for i in range(n):dp[1<<i][i]=1
 for S in range(1,1<<n):
  for i,v in enumerate(dp[S]):
   if v:
    z=a[i]&~S
    while z:
     b=z&-z;z-=b;j=b.bit_length()-1;dp[S|b][j]+=v
 return sum(dp[-1])

def cyclic(blocks):
 sizes=list(map(len,blocks));off=[0,sizes[0],sizes[0]+sizes[1]];out=[]
 for i,A in enumerate(blocks):
  j=(i+1)%3;cross=((1<<sizes[j])-1)<<off[j]
  out.extend((v<<off[i])|cross for v in A)
 return tuple(out)

def strong(a):
 n=len(a)
 for r in range(n):
  seen=1<<r
  while True:
   nxt=seen
   for i in range(n):
    if seen>>i&1:nxt|=a[i]
   if nxt==seen:break
   seen=nxt
  if seen!=(1<<n)-1:return False
 return True
seeds=json.loads((P/'seed_controls.json').read_text());base=(58,60,57,16,32,8);accepted=[]
for row in seeds['rows']:
 r=certify_actual(seeds['original_masks'],row['child_masks'],base)
 if r['actual_profile_difference']!=row['difference']:raise RuntimeError('seed profile')
 if row['cut']==996:
  if r['positive']['status']!='SUFFICIENT_CERTIFICATE' or r['positive']['lower_bound_coefficient']!=135618:raise RuntimeError('positive certificate')
 else:
  if r['negative']['status']!='SUFFICIENT_CERTIFICATE' or r['negative']['pivot']!={867:4,995:5}[row['cut']]:raise RuntimeError('negative certificate')
 accepted.append({'cut':row['cut'],'certificate':r['positive'] if row['cut']==996 else r['negative']})
D=json.loads((P/'actual_row_hostile.json').read_text());old=tuple(D['old_module_masks']);new=tuple(D['new_module_masks']);B=D['parent_masks']
for A in [old,new]:
 validate_tournament(A)
 if not strong(A) or any((A[i]&63)!=B[i] for i in range(6)):raise RuntimeError('strong same-card hostile')
f=literal_profile(old);g=literal_profile(new);d=[b-a for a,b in zip(f,g)]
if f!=D['old_profile'] or g!=D['new_profile'] or f!=covers(old) or g!=covers(new):raise RuntimeError('literal module profiles')
K=cyclic_kernel(d)
if negative_certificate(d,[1])['status']!='NOT_CERTIFIED' or quadratic(K,[1])!=-36:raise RuntimeError('hostile rejected by row criterion')
contexts=[]
for side in [(0,),(2,0),(2,4,1)]:
 A=cyclic((old,side,side));C=cyclic((new,side,side));h0=h_forward(A);h1=h_forward(C);v=covers(side);pred=quadratic(K,v)
 if h1-h0!=pred:raise RuntimeError('independent literal cyclic context')
 contexts.append({'side_order':len(side),'old_H':h0,'new_H':h1,'gap':pred})
if [r['gap'] for r in contexts]!=[-36,828,67932]:raise RuntimeError('actual missing-row hostile')
v=[2,1,2];L=[r*((v[r-1] if r<=len(v) else 0)+(v[r-2] if r>=2 else 0)) for r in range(1,len(v)+2)];ratios=[str(Fraction(L[i],v[i])) for i in range(3)]
if ratios!=['1','6','9/2']:raise RuntimeError('formal missing-log-concavity control')
report={'certified_actual_controls':accepted,'hostile':{'module_order':7,'difference':d,'base_profile':[1],'TT2_multiplier_power_coefficients':[0,1],'C3_multiplier_power_coefficients':[2,0,1],'contexts':contexts},'formal_log_concavity_boundary':{'base':v,'multiply_by_m':L,'ratios':ratios}}
print(json.dumps(report,sort_keys=True));print('PASS sha256='+hashlib.sha256(json.dumps(report,separators=(',',':')).encode()).hexdigest())
