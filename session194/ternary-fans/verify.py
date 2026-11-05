"""Independent certificate for the Session194 ternary join fans.
Run in this directory: python3 verify.py and python3 -O verify.py.
No source-bank counting code is imported. All output is deterministic.
"""
import itertools,json,math
from pathlib import Path

def require(c,m):
 if not c: raise RuntimeError(m)

def validate(out):
 n=len(out)
 require(all(0<=a<(1<<n) and not(a>>i&1) for i,a in enumerate(out)),'bad row')
 require(all(((out[i]>>j&1)+(out[j]>>i&1))==1 for i in range(n) for j in range(i)),'bad edge')

def dp(out):
 n=len(out);full=(1<<n)-1
 ends=[[0]*n for _ in range(full+1)]
 hs=[0]*(full+1);hs[0]=1
 incoming=[sum(1<<j for j in range(n) if out[j]>>i&1) for i in range(n)]
 for S in range(1,full+1):
  vs=S
  while vs:
   z=vs&-vs;v=z.bit_length()-1;vs-=z;T=S^z
   if not T: ends[S][v]=1
   else:
    us=T&incoming[v]
    while us:
     u=us&-us;us-=u;ends[S][v]+=ends[T][u.bit_length()-1]
  hs[S]=sum(ends[S])
 return hs,ends

def triple(out,r):
 h,e=dp(out);full=(1<<len(out))-1
 rev=[sum(1<<j for j,x in enumerate(out) if x>>i&1) for i in range(len(out))]
 _,s=dp(rev)
 R=sum(h[S]*s[full^S][r] for S in range(full+1) if not(S>>r&1))
 return [h[full],e[full][r],R]

def literal_profile(out,r):
 h=e=R=0;n=len(out)
 for p in itertools.permutations(range(n)):
  ok=[bool(out[p[i]]>>p[i+1]&1) for i in range(n-1)]
  if all(ok): h+=1;e+=int(p[-1]==r)
  cut=p.index(r)
  if all(ok[i] for i in range(n-1) if i!=cut-1):R+=1
 return [h,e,R]

def strong(out):
 n=len(out);full=(1<<n)-1
 for start in range(n):
  old=0;seen=1<<start
  while old!=seen:
   old=seen
   for i in range(n):
    if seen>>i&1:seen|=out[i]
  if seen!=full:return False
 return True

def g(v):
 h=1
 while v:
  u=v//2
  if not(v&1):h+=u+(u&1)
  v=u
 return h

def apex(v):
 bits=list(map(int,bin(v)[2:])) if v else []
 n=len(bits);out=[0]*(n+1)
 for i,b in enumerate(bits):
  out[i]=sum(1<<j for j in range(i+1,n))
  if b:out[n]|=1<<i
  else:out[i]|=1<<n
 return out,n

def join(A,a,B,b):
 n=len(A);m=len(B);out=list(A)+[0]*m
 for i in range(n):out[i]|=((1<<m)-1)<<n
 for i in range(m):
  for j in range(m):
   if B[j]>>i&1:out[n+i]|=1<<(n+j)
 out[a]^=1<<(n+b);out[n+b]|=1<<a
 return out

def birth(B,mask):
 n=len(B)
 return [x+(0 if mask>>i&1 else 1<<n) for i,x in enumerate(B)]+[mask]

fans=json.loads(Path('witnesses.json').read_text())
profile_checks=0;direct_checks=0;literal_permutations=0
for fan in fans:
 p=fan['parent'];h,e,R=p['triple']
 for rec in [p,*fan['children'].values()]:
  out=rec['out'];r=rec['root'];validate(out)
  require(triple(out,r)==rec['triple'],'DP profile')
  require(literal_profile(out,r)==rec['triple'],'literal profile')
  literal_permutations+=math.factorial(len(out));profile_checks+=1
 for d in [2,4]:
  c=fan['children'][str(d)]
  require(c['triple']==[3*h,3*e-d,3*R],'coefficient match')
 for v in [0,1,2,3,4,5,6]:
  A,a=apex(v);J=join(A,a,p['out'],p['root']);validate(J)
  N=dp(J)[0][-1]
  require(N==h*g(v)-e*int(v%2==0)+R*(v+1),'parent join')
  direct_checks+=1
  for d in [2,4]:
   c=fan['children'][str(d)];K=join(A,a,c['out'],c['root']);validate(K)
   require(dp(K)[0][-1]==3*N+d*int(v%2==0),'child join')
   direct_checks+=1
special=next(f for f in fans if f['parent']['triple']==[57,15,57])
p=special['parent'];require(p['root']==1,'root convention')
require(strong(p['out']),'parent strong')
for d in [2,4]:
 c=special['children'][str(d)]
 require(c['out']==birth(p['out'],c['birth_outmask']),'physical birth')
 require(c['root']==p['root'],'mark preserved')
 require(strong(c['out']),'child strong')
# Exact all-digit arithmetic and induced graph extension controls.
for v in range(0,65536,2):
 N=57*g(v)-15+57*(v+1)
 require(N==57*g(2*v)+42,'even simplification')
 for d in [2,4]:
  require(171*g(v)-(45-d)+171*(v+1)==3*N+d,'all-digit arithmetic')
special_cases=[]
for v in range(16):
 A,a=apex(v);J=join(A,a,p['out'],p['root']);h,ends=dp(J);N=h[-1]
 require(strong(J),'joined parent strong')
 for d in [2,4]:
  c=special['children'][str(d)];K=join(A,a,c['out'],c['root']);hk,_=dp(K)
  require([x&((1<<len(J))-1) for x in K[:-1]]==J,'old arcs preserved')
  require(hk[-1]==3*N+d*int(v%2==0),'induced count')
  require(strong(K),'joined child strong');direct_checks+=1
 special_cases.append({'v':v,'N':N,'lower':3*N+2*int(v%2==0),'upper':3*N+4*int(v%2==0),'parent_endpoint_counts':ends[-1]})
# Non-apex inputs: the identity has digit d*e_a for every rooted A.
input_graphs=[[0],[2,4,1],[2,8,3,5],[16,9,11,49,38,71,31]]
general_checks=0
for A in input_graphs:
 ha,ea=dp(A)
 for a in range(len(A)):
  J=join(A,a,p['out'],p['root']);N=dp(J)[0][-1]
  for d in [2,4]:
   c=special['children'][str(d)];K=join(A,a,c['out'],c['root'])
   require(dp(K)[0][-1]==3*N+d*ea[-1][a],'general input identity')
   general_checks+=1
report={'status':'PROVED elementary infinite families with FINITE-EXACT certificate checks; global spectrum OPEN','physical_fans':len(fans),'independently_checked_profiles':profile_checks,'literal_permutation_iterations':literal_permutations,'direct_one_apex_join_checks':direct_checks,'general_rooted_input_checks':general_checks,'even_arithmetic_parameters':32768,'special_parent':p,'special_children':special['children'],'special_cases':special_cases,'coverage_limit':'No all-order union/prime-lane coverage or return-interface claim. Source bank completeness is not needed. The 2^32 claim is not used.'}
Path('receipt.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2))
