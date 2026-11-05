"""Supplement to verify.py: C2 explanation and a paired large-atom certificate.
Run: python3 verify_structure_and_atoms.py (also -O).
The large graph counts use proved join/apex identities, not exponential subset DP.
"""
import itertools,json,math
from pathlib import Path

def require(c,m):
 if not c: raise RuntimeError(m)

def profile(o,r):
 n=len(o);h=e=R=0
 for p in itertools.permutations(range(n)):
  bad=[i for i in range(n-1) if not(o[p[i]]>>p[i+1]&1)]
  if not bad:h+=1;e+=int(p[-1]==r)
  cut=p.index(r)-1
  if not bad or bad==[cut]:R+=1
 return [h,e,R]

def make_c2(U,r):
 # Output order y,x, followed by U; x's sole old out-neighbor is r.
 n=len(U);o=[2,1<<(r+2)]+[0]*n
 for i in range(n):o[i+2]=(U[i]<<2)|1|(0 if i==r else 2)
 return o

U=[24,9,3,20,6]
structure=[]
for d,mask in [(0,None),(2,14),(4,6)]:
 V=U if mask is None else [x+(0 if mask>>i&1 else 32) for i,x in enumerate(U)]+[mask]
 P=profile(V,0);want={0:[15,3,21],2:[43,6,64],4:[41,5,65]}[d]
 require(P==want,'base profile')
 C=make_c2(V,0);pc=profile(C,1)
 require(pc==[P[0]+2*P[2],P[0],P[0]+2*P[2]],'C2 profile')
 structure.append({'digit':d,'out':V,'root':0,'profile':P,'C2_out':C,'C2_profile':pc})
fans=json.loads(Path('witnesses.json').read_text());f=next(x for x in fans if x['parent']['triple']==[57,15,57])
require(structure[0]['C2_out']==f['parent']['out'],'parent names')
for x in structure[1:]:require(x['C2_out']==f['children'][str(x['digit'])]['out'],'child names')

def recurrence(v):
 h=1;u=0
 for bit in bin(v)[2:]:
  if bit=='0':h+=u+(u&1)
  u=2*u+int(bit)
 return h

def pairs(v):
 w=list(map(int,bin(v)[2:]));n=len(w)
 return 1+sum(w[i]*(1-w[i+1]) for i in range(n-1))+sum(w[i]*(1-w[j])*(1<<(j-i-1)) for i in range(n) for j in range(i+1,n))

v=33554594;gv=recurrence(v);g2v=recurrence(2*v)
require(gv==pairs(v) and g2v==pairs(2*v),'independent apex formula')
require(g2v==gv+v,'even recurrence')
N=57*gv-15+57*(v+1);lo=3*N+2;hi=3*N+4
require((N,lo,hi)==(3309715245,9929145737,9929145739),'large counts')
require(lo==7*1418449391,'seven-prime factor')

def trial_prime(n):
 require(n>=3 and n%2==1,'prime input')
 bound=math.isqrt(n);cnt=0
 for d in range(3,bound+1,2):
  cnt+=1
  require(n%d!=0,('composite',n,d))
 return {'n':n,'method':'exhaustive odd trial division from 3 through floor(sqrt(n))','square_root_floor':bound,'trial_divisors_checked':cnt}

primality=[trial_prime(1418449391),trial_prime(hi)]
w=list(map(int,bin(v)[2:]));n=len(w);A=[[False]*(n+1) for _ in range(n+1)]
for i in range(n):
 for j in range(i+1,n):A[i][j]=True
 A[n][i]=bool(w[i]);A[i][n]=not bool(w[i])

def concrete_join(B):
 a=n;b=1;na=n+1;nb=len(B);m=na+nb;M=[[False]*m for _ in range(m)]
 for i in range(na):
  for j in range(na):M[i][j]=A[i][j]
  for j in range(nb):M[i][na+j]=True
 for i in range(nb):
  for j in range(nb):M[na+i][na+j]=bool(B[j]>>i&1)
 M[a][na+b]=False;M[na+b][a]=True
 require(all(not M[i][i] for i in range(m)),'loop')
 require(all(M[i][j]!=M[j][i] for i in range(m) for j in range(i)),'edge')
 for start in range(m):
  reached={start};todo=[start]
  while todo:
   k=todo.pop()
   for j in range(m):
    if M[k][j] and j not in reached:reached.add(j);todo.append(j)
  require(len(reached)==m,'not strong')
 return [sum(1<<j for j in range(m) if M[i][j]) for i in range(m)]
J=concrete_join(f['parent']['out']);children={}
for d in [2,4]:
 K=concrete_join(f['children'][str(d)]['out'])
 require([r&((1<<len(J))-1) for r in K[:-1]]==J,'induced large parent')
 children[str(d)]={'order':len(K),'H':3*N+d,'out':K}
endpoint_certificates=[]
for digit,rec in [(0,f['parent']),*( (d,f['children'][str(d)]) for d in [2,4])]:
 o=rec['out'];b=rec['root'];m=len(o);ss=[0]*m;cc=[0]*m;DD=[0]*m
 for p in itertools.permutations(range(m)):
  bad=[i for i in range(m-1) if not(o[p[i]]>>p[i+1]&1)]
  if not bad:
   ss[p[0]]+=1
   if p[-1]==b:cc[p[0]]+=1
  if p[0]==b:
   for cut in range(1,m):
    if not bad or bad==[cut-1]:DD[p[cut]]+=1
 q=3 if digit==0 else 7
 require(ss[b]==q,'left endpoint factor')
 require((ss[2],cc[2],DD[2])==(q,q,q),'exceptional endpoint')
 require(all(ss[x]-cc[x]>=q for x in range(m) if x!=2),'right endpoint lower bound')
 endpoint_certificates.append({'digit':digit,'q':q,'s':ss,'c':cc,'D':DD})
report={'uniform_no_endpoint_one_certificates':endpoint_certificates,'status':'FINITE-EXACT physical cores and primality; PROVED algebraic graph counts','C2_structure':structure,'large_pair':{'v':v,'binary_word':bin(v)[2:],'g_v':gv,'g_2v':g2v,'parent_order':len(J),'parent_H':N,'parent_out':J,'children':children,'prime_proofs':primality,'count_method':'The exact rooted join law, independently certified core profiles, and the all-word one-apex formula; not a 35-vertex subset-DP computation.','scope':'Two explicit atom witnesses beyond 2^32, not an interval-coverage claim.'}}
Path('structure_and_atoms_receipt.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
