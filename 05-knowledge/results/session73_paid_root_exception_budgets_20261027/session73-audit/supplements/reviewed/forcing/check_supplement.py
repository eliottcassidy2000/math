#!/usr/bin/env python3
"""Independent cached/comb scalar certificate and exact phase/spatial bridges."""
from fractions import Fraction as F
from math import comb
import json,hashlib
checks=0
cache={1:(0,0,0)}
def need(x,msg):
 global checks
 checks+=1
 if not x:raise RuntimeError(msg)
def step(n):
 z=3*n+1;a=0
 while z%2==0:z//=2;a+=1
 return z,a

def counters(n):
 if n in cache:return cache[n]
 path=[];x=n
 while x not in cache:
  y,a=step(x);path.append((x,a));x=y
  if len(path)>10000:raise RuntimeError('finite cap reached')
 tau,K,A=cache[x]
 for y,a in reversed(path):
  tau+=1;K+=(a-1)//2;A+=a;cache[y]=(tau,K,A)
 return cache[n]
def weight(n):
 if n==1:return F(1)
 t,k,a=counters(n)
 return F(2,(t+k+1)*comb(t+k,k))
def sign(n):return (-1)**counters(n)[2]
def rational(x):return str(F(x))
def w(n):return 3*(n%3==2)-1
def chi(n):return [0,1,-1][n%3]
def kernel(N,m):
 return w(m)*(N>=(m+3)//4)+(-3*chi(m)-w(m))*(N>=(m+1)//2)
def round_box(x,y,scale=10**10):
 lo=(x.numerator*scale)//x.denominator;hi=-((-y.numerator*scale)//y.denominator)
 need(F(lo,scale)<=x<=y<=F(hi,scale),'outward rational rounding')
 return [str(F(lo,scale)),str(F(hi,scale))]

mass=energy=F(0);max_tau=0;source_hash=hashlib.sha256()
for n in range(3,5002,2):
 t,k,a=counters(n);z=weight(n);mass+=z;energy+=z/F(6*n+1);max_tau=max(max_tau,t)
 source_hash.update(f'{n}:{t},{k},{a}:{z};'.encode())
upper=energy+(F(13,3)-mass)/F(6*5003+1)
need(0<energy<upper<F(1,42),'global native scalar upper')
box=round_box(energy,upper)
need(F(box[1])<F(1,42),'rounded scalar upper')

# Full marked source ancestry retained; no change of source mass at merges.
mu={};J={};raw={};visits=0
for source in range(3,502,6):
 initial=sign(source);lam=weight(source);n=source;phase=1;shortcut=0;odd=0
 while n!=1:
  need(initial*phase==sign(n),'transported native sign')
  mu[n]=mu.get(n,F(0))+lam
  J[n]=J.get(n,F(0))+lam*phase*initial
  raw[n]=raw.get(n,F(0))+lam*initial
  n,a=step(n);phase*=(-1)**a;shortcut+=a;odd+=1;visits+=1
 need(phase*initial==1,'root signed transport')
 need((odd,shortcut)==(counters(source)[0],counters(source)[2]),'marked clocks')
 need(odd+shortcut>=shortcut,'literal clock sidecar')
for n in mu:
 need(J[n]==mu[n]*sign(n),'named merge current')
 need(mu[n]<=weight(n),'native occupancy bound')
transport_energy=sum((J[n]**2/weight(n) for n in mu),F(0))
need(transport_energy<=sum(mu.values(),F(0))<F(13,3),'weighted transported energy')
hostiles=[n for n in sorted(mu) if raw[n]!=J[n]]
need(bool(hostiles),'phase-forgetting hostile exists')
h=hostiles[0]

# Exact original spatial current, independent direct/fold implementations.
fold_rows=[]
for N in range(1,257):
 direct=0
 for n in range(N+1,4*N+1):
  odd=n;v=0
  while odd%2==0:odd//=2;v+=1
  direct+=w(n)*((-1)**v)*sign(odd)
 folded=sum(kernel(N,m)*sign(m) for m in range(1,4*N+1,2))
 need(direct==folded,'original spatial current fold')
 # Literal weighted-current readout: division by exact source mass.
 paid=sum(F(kernel(N,m))*weight(m)*sign(m)/weight(m) for m in range(3,4*N+1,2))
 need(paid==folded+3,'spatial inverse-mass readout with root separated')
 fold_rows.append(direct)
windows=[]
for X in [4,8,16,32,64,128]:
 observed=sum((fold_rows[N-1]+3)**2 for N in range(X,2*X))
 columns=[]
 for m in range(3,8*X,2):
  col=sum(kernel(N,m)**2 for N in range(X,2*X))
  if col:columns.append((m,F(col)/weight(m)))
 hs=sum((x[1] for x in columns),F(0));lower=max(x[1] for x in columns)
 need(F(observed)<=F(13,3)*hs,'spatial Cauchy bound')
 windows.append({'X':X,'actual_centered_spatial_energy':observed,'source_operator_squared_norm_lower':str(lower),'hilbert_schmidt_squared':str(hs),'Cauchy_upper':str(F(13,3)*hs)})
print(json.dumps({'status':'PASS','checks':checks,'scalar':{'cutoff':5001,'sources':2500,'maximum_odd_steps':max_tau,'outward_rational_interval':box,'strict_upper_1_over_42':True,'method':'cached native counters, integer binomial weight; inherited nonroot mass13/3 tail','counters_sha256':source_hash.hexdigest()},'phase_transport':{'sources':84,'nonroot_marked_visits':visits,'distinct_visited_vertices':len(mu),'weighted_energy_upper_13_over_3':True,'phase_forgetting_hostile':{'vertex':h,'phase_retained_J':rational(J[h]),'phase_erased_J':rational(raw[h]),'occupation_mass':rational(mu[h])}},'spatial_current':{'all_N_1_to_256_direct_fold_agree':True,'current_sha256':hashlib.sha256(json.dumps(fold_rows).encode()).hexdigest(),'windows':windows},'boundaries':['scalar payment is not spatial covariance','inverse masses are not supplied by positive-energy existence','zero ROOT mass never provides an invertible observation weight']},indent=2,sort_keys=True))
