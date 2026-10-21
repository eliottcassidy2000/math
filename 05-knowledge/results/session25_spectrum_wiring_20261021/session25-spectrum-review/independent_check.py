from itertools import product, permutations
from math import gcd
import json

def require(ok,msg):
 if not ok: raise RuntimeError(msg)

def adjacency(w):
 n=len(w); a=[[False]*(n+1) for _ in range(n+1)]
 for i in range(n):
  for j in range(i+1,n): a[i][j]=True
  a[n][i]=bool(w[i]);a[i][n]=not w[i]
 return a

def dp(w):
 a=adjacency(w);n=len(a);D=[{} for _ in range(1<<n)]
 for v in range(n):D[1<<v][v]=1
 for s in range(1,1<<n):
  for v,c in D[s].items():
   for u in range(n):
    if not(s>>u&1) and a[v][u]:
     z=D[s|1<<u];z[u]=z.get(u,0)+c
 return sum(D[-1].values())

def strong(w):
 a=adjacency(w);n=len(a)
 for v in range(n):
  seen={v};todo=[v]
  while todo:
   i=todo.pop()
   for j in range(n):
    if a[i][j] and j not in seen:seen.add(j);todo.append(j)
  if len(seen)!=n:return False
 return True

def split_count(w):
 # Independently enumerate all assignments before/after apex, not recurrence.
 n=len(w);total=0
 for mask in range(1<<n):
  left=[i for i in range(n) if mask>>i&1]
  right=[i for i in range(n) if not(mask>>i&1)]
  total+=bool((not left or not w[left[-1]]) and (not right or w[right[0]]))
 return total

def pair_formula(w):
 if not w:return 1
 return w[0]+1-w[-1]+sum(w[i]==0 and w[i+1]==1 for i in range(len(w)-1))+sum(1<<(j-i-1) for i in range(len(w)) for j in range(i+1,len(w)) if w[i] and not w[j])

def recurrence(w,mod=None):
 h,a,p=1,0,0
 for bit in w:
  if not bit:h+=a+p
  a,p=2*a+bit,bit
  if mod:h%=mod;a%=mod
 return h

rows=0;literals=0
for n in range(9):
 for w in product((0,1),repeat=n):
  x=dp(w)
  require(x==split_count(w)==pair_formula(w)==recurrence(w),('count',w,x))
  require(x%2==1,('parity',w))
  if n>=1:require(strong(w)==(n>=2 and w[0]==1 and w[-1]==0),('strong',w))
  if n<=5:
   a=adjacency(w)
   z=sum(all(a[p[i]][p[i+1]] for i in range(n)) for p in permutations(range(n+1)))
   require(z==x,('literal',w));literals+=1
  rows+=1

# Independent binomial/exponential closed form, evaluated as integer before reducing.
closed=0
for L in range(1,8):
 q=1<<L
 for t in range(0,16):
  numerator=(q+1)**2+q*q*(q**(2*t)-1)+2*t*(q+1)
  require(numerator%(q+1)**2==0,('integral',L,t))
  exact=numerator//(q+1)**2
  require(exact==recurrence((1,)*L+(0,)*L) if t==1 else exact==recurrence(((1,)*L+(0,)*L)*t),('closed',L,t))
  closed+=1

# CRT multiplier checks across all moduli <=1024. Coverage is tested independently
# using block state dynamics, never asserting recurrence from observations.
modrows=0;wordchecks=0;maxL=0
for M in range(1,1025):
 odd=M;a=0
 while odd%2==0:odd//=2;a+=1
 L=max(1,a)
 while pow(2,L,odd)!=1%odd:L+=1
 maxL=max(maxL,L);q=pow(2,L,M);d=((q-1)**2+1)%M;g=gcd(2,M)
 require(gcd(d,M)==g,('gcd',M,L,d))
 h,A=1%M,0;found=set();period=M//g
 for t in range(1,period+1):
  h=(h+q*(q-1)*A+(q-1)**2+1)%M
  A=(q*q*A+q*(q-1))%M
  require(A==0,('reset',M,t))
  require(h==(1+t*d)%M,('linear',M,t))
  found.add(h);modrows+=1
  if M<=40 and t in {1,period,period//2}:
   w=((1,)*L+(0,)*L)*t
   require(recurrence(w,M)==h,('wordmod',M,t));wordchecks+=1
 require(found=={r for r in range(M) if (r-1)%g==0},('coverage',M))

# New sparse-family and exact untuned-modulus corollaries.
for L in range(1,11):
 for t in range(1,31):
  h=recurrence(((1,)*L+(0,)*L)*t)
  require(3*h>=4**(L*t),('sparse lower bound',L,t))
  if L==1: require(h%3==(1+2*t*t)%3,('quadratic mod3',t))

hostiles={
 'L0_not_fvs1': {'word':[], 'H':dp(()), 'feedback_vertex_number':0},
 't0_not_fvs1': {'word':[], 'H':dp(()), 'feedback_vertex_number':0},
 'insufficient_dyadic_L': {'M':8,'L':1,'t':2,'actual':recurrence((1,0)*2,8),'claimed_linear':(1+2*2)%8},
 'missing_odd_order_condition': {'M':3,'L':1,'t':2,'actual':recurrence((1,0)*2,3),'claimed_linear':(1+2*2)%3},
 'residue7_not_exact7': {'M':8,'L':3,'t':3,'H':recurrence(((1,)*3+(0,)*3)*3)}
}
require(hostiles['insufficient_dyadic_L']['actual']!=hostiles['insufficient_dyadic_L']['claimed_linear'],'hostile8')
require(hostiles['missing_odd_order_condition']['actual']!=hostiles['missing_odd_order_condition']['claimed_linear'],'hostile3')
require(hostiles['residue7_not_exact7']['H']%8==7 and hostiles['residue7_not_exact7']['H']!=7,'residue hostile')
print(json.dumps({'status':'PASS','all_word_DP_split_pair_recurrence_controls':rows,'literal_controls':literals,'closed_form_controls':closed,'sparse_lower_bound_controls':300,'untuned_mod3_controls':30,'moduli':1024,'block_state_controls':modrows,'direct_modular_word_controls':wordchecks,'max_L':maxL,'hostiles':hostiles},sort_keys=True,indent=2))
