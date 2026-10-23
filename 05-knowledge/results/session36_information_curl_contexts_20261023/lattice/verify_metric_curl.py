"""Exact certificates. Python standard library + SymPy. No assert dependence.
Run from workspace: python3 session36-lattice/verify_metric_curl.py
"""
from itertools import combinations,permutations
from collections import Counter,defaultdict
from math import factorial
from random import Random
from pathlib import Path
import json,hashlib
from sympy import Matrix,ones,eye
from explore import matrix,Kkey,hamilton,scc,rootpair

def check(cond,msg):
 if not cond: raise RuntimeError(msg)
def switch(M,mask):
 n=len(M);d=[1]+[-1 if mask>>(i-1)&1 else 1 for i in range(1,n)]
 return [[d[i]*M[i][j]*d[j] for j in range(n)] for i in range(n)]
def conference(q):
 n=q+1;Q={i*i%q for i in range(1,q)}
 M=[[0]*n for _ in range(n)]
 for i in range(n):
  for j in range(i+1,n):
   M[i][j]=1 if i==0 or (j-i)%q in Q else -1;M[j][i]=-M[i][j]
 check(Matrix(M)**2==-(n-1)*eye(n),('conference square',q))
 return M

def brute_H(M):
 n=len(M);return sum(all(M[p[i]][p[i+1]]==1 for i in range(n-1)) for p in permutations(range(n)))

def root_data(M):
 n=len(M);roots=list(permutations(range(n),2))
 G=[[(i==k)-(i==l)-(j==k)+(j==l) for k,l in roots] for i,j in roots]
 W=[[M[i][k]-M[i][l]-M[j][k]+M[j][l] for k,l in roots] for i,j in roots]
 return roots,G,W

def triple_coefficient(M):
 roots,G,W=root_data(M);N=len(roots);total=0
 for a in range(N):
  for b in range(N):
   if G[a][b]!=-1 or W[a][b]!=-3:continue
   for c in range(N):
    total+=G[a][c]==0 and G[b][c]==0 and W[a][c]==-4 and W[b][c]==2
 return total

def direct_five_vertex(M):
 # Independent path: roots x=e_a-e_b,y=e_b-e_c,z=e_d-e_e, and
 # roots x=e_a-e_b,y=e_c-e_a,z=e_d-e_e. Enumerate vertices, not roots.
 n=len(M);out=0
 for a,b,c,d,e in permutations(range(n),5):
  wab=M[a][b]-M[a][c]+M[b][c]
  wac=M[a][d]-M[a][e]-M[b][d]+M[b][e]
  wbc=M[b][d]-M[b][e]-M[c][d]+M[c][e]
  out+=(wab,wac,wbc)==(-3,-4,2)
  # second geometry: y=e_c-e_a
  wab2=M[a][c]-M[b][c]+M[b][a]
  wbc2=M[c][d]-M[c][e]-M[a][d]+M[a][e]
  out+=(wab2,wac,wbc2)==(-3,-4,2)
 return out

def pair_formula(n,c3,c4):
 A=n*(n-1)*(n-2);z={(-2,0):n*(n-1),(2,0):n*(n-1)}
 for g in [-1,1]:
  for w in [-3,3]:z[g,w]=6*c3
  for w in [-1,1]:z[g,w]=A-6*c3
 for w in [-4,4]:z[0,w]=4*c4
 for w in [-2,2]:z[0,w]=A*(n-3)//6+8*c3*(n-3)-16*c4
 z[0,0]=2*A*(n-3)//3-16*c3*(n-3)+24*c4
 return {k:v for k,v in z.items() if v}
report={};census={};total_t=0;all_pair_checks=0
for n in range(1,7):
 fibres=defaultdict(list)
 for m in range(1<<(n*(n-1)//2)):
  M=matrix(n,m);C=scc(M);h=hamilton(M);key=Kkey(M)
  if n<=5:
   check(h==brute_H(M),('independent H',n,m))
   scores=[sum(x==1 for x in row) for row in M]
   c3=n*(n-1)*(n-2)//6-sum(x*(x-1)//2 for x in scores)
   c4=sum(all(M[p[i]][p[(i+1)%4]]==1 for i in range(4)) for p in permutations(range(n),4))//4
   check(dict(rootpair(M))==pair_formula(n,c3,c4),('general pair formula',n,m,c3,c4))
   all_pair_checks+=1
  fibres[key].append((m,len(C),h))
  # Every SCC rotation has the exact same labelled curl.
  for k in range(len(C)):
   X=set(v for c in C[:k] for v in c)
   U=[[(-M[i][j] if (i in X)!=(j in X) else M[i][j]) for j in range(n)] for i in range(n)]
   check(Kkey(U)==key,('rotation',n,m,k))
 for f in fibres.values():
  check(len({x[2] for x in f})==1,('H fibre',n,f))
  check(all(x[1]==len(f) for x in f),('size fibre',n,f))
  # Every distinct collision is exactly a two-level cut difference.
  M=matrix(n,f[0][0])
  for m,_,_ in f[1:]:
   U=matrix(n,m);a=[M[i][-1]-U[i][-1] for i in range(n)]
   check(len(set(a))==2 and max(a)-min(a)==2,('two levels',n,f[0][0],m))
   check(all(M[i][j]-U[i][j]==a[i]-a[j] for i in range(n) for j in range(n)),('gradient',n,m))
 census[n]={'tournaments':sum(map(len,fibres.values())),'forms':len(fibres),'fibre_sizes':dict(sorted(Counter(map(len,fibres.values())).items()))}
 total_t+=census[n]['tournaments']
 print('census',n,census[n],flush=True)
report['labelled_census']=census;report['total_tournaments']=total_t;report['all_tournament_pair_controls']=all_pair_checks

# Actual ordered root-star frame enumeration, without presuming classification.
frame_counts={}
for n in range(2,7):
 roots,G,W=root_data(matrix(n,0));N=len(roots);count=[0]
 def visit(chosen):
  if len(chosen)==n-1:count[0]+=1;return
  for a in range(N):
   if all(G[a][b]==1 for b in chosen):visit(chosen+[a])
 visit([])
 expected=2 if n==2 else 2*factorial(n)
 check(count[0]==expected,('frame count',n,count[0],expected))
 frame_counts[n]=count[0]
report['root_basis_counts']=frame_counts

# All128 specified switches, not inferred from representative types.
M0=conference(7);base_pairs=rootpair(M0);cs=[]
for mask in range(128):
 M=switch(M0,mask);pairs=rootpair(M);c=triple_coefficient(M);direct=direct_five_vertex(M);h=hamilton(M)
 check(pairs==base_pairs,('genus2',mask))
 check(c==direct,('two genus3 paths',mask,c,direct))
 universal=any(abs(sum(row))==7 for row in M)
 check(c==(84 if universal else 36),('coefficient',mask,c,universal))
 check(h==441-3*c,('H linear within family',mask,h,c))
 cs.append({'mask':mask,'H':h,'coefficient':c,'universal':universal})
report['paley8_switches']=cs
report['paley8_genus2_distribution']=[{'inner_product':k[0],'alternating_pairing':k[1],'count':v} for k,v in base_pairs]
print('all128 switches: genus2 identical; genus3 counts',dict(Counter(x['coefficient'] for x in cs)),flush=True)

# General conference distribution controls at multiple orders and switching gauges.
conference_checks=[]
for q in [3,7,11,19]:
 n=q+1;M=conference(q);A=n*(n-1)*(n-2)
 expected={(-2,0):n*(n-1),(2,0):n*(n-1)}
 for g in [-1,1]:
  for w in [-3,3]:expected[g,w]=A//4
  for w in [-1,1]:expected[g,w]=3*A//4
 for w in [-4,4]:expected[0,w]=A*(n-4)//16
 for w in [-2,2]:expected[0,w]=A*(n-2)//4
 expected[0,0]=3*A*(n-4)//8
 expected={k:v for k,v in expected.items() if v}
 for mask in [0,1,3,(1<<(n-1))-1]:
  got=dict(rootpair(switch(M,mask)))
  check(got==expected,('universal rootpair',n,mask,got,expected))
  conference_checks.append({'n':n,'mask':mask})
report['general_conference_controls']=conference_checks
print('general conference controls',len(conference_checks),flush=True)

# A hostile to extrapolating the selected single genus3 coefficient.
M12=conference(11);single12=[]
for mask in [0,1]:
 M=switch(M12,mask);c=triple_coefficient(M);h=hamilton(M)
 check(c==660,('single coefficient order12',mask,c))
 single12.append({'mask':mask,'coefficient':c,'H':h,'fourth_row_sum_moment':sum(sum(row)**4 for row in M)})
check([x['H'] for x in single12]==[95095,196815],'order12 H hostile')
report['single_coefficient_order12_hostile']=single12

# Exact polynomial identity by enough evaluation points; all determinant arithmetic integer.
rng=Random(36035);pencil_checks=0
for n in range(2,10):
 for rep in range(3):
  M=Matrix(matrix(n,rng.randrange(1<<(n*(n-1)//2))))
  B=Matrix(n,n-1,lambda i,j:1 if i==j else (-1 if i==n-1 else 0))
  G=B.T*B;K=B.T*M*B
  for t in range(n):
   X=t*eye(n)-M
   # rank-one determinant identity gives q(t) with no symbolic adjugate.
   q=(X+ones(n)).det()-X.det()
   check((t*G-K).det()==q,('pencil',n,rep,t))
   pencil_checks+=1
report['pencil_integer_evaluations']=pencil_checks

# Metric-free Smith collision inherited from Session35, recomputed.
from sympy.matrices.normalforms import smith_normal_form
from sympy import ZZ
smith=[]
for out in [[18,4,9,19,6],[30,4,8,18,6]]:
 M=[[0 if i==j else (1 if out[i]>>j&1 else -1) for j in range(5)] for i in range(5)]
 B=Matrix(5,4,lambda i,j:1 if i==j else (-1 if i==4 else 0));K=B.T*Matrix(M)*B
 snf=smith_normal_form(K,domain=ZZ)
 ds=[abs(int(snf[i,i])) for i in range(4)]
 check(ds==[1,1,1,1],('Smith witness',out,ds))
 smith.append({'outmasks':out,'H':hamilton(M),'Smith':ds,'K':list(K)})
report['Smith_hostile']=smith
check([x['H'] for x in smith]==[11,5],'Smith H')
report['status']='PASS'
p=Path(__file__).with_suffix('.json');p.write_text(json.dumps(report,indent=2,default=int)+'\n')
print('PASS',total_t,'tournaments;',pencil_checks,'pencil evaluations; all128 two-path genus3 checks',flush=True)
