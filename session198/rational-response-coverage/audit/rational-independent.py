"""Independent integer g evaluation, clock coefficients, and literal controls."""
from pathlib import Path
import importlib.util,json,math
HERE=Path(__file__).resolve().parent
P=HERE/'rational-replay'
spec=importlib.util.spec_from_file_location('frozen_rational_constructor',P/'verify.py');m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)

def g(n):
 ans=1
 while n>1:
  if n%2==0:n//=2
  else:
   n=(n+1)//2;ans+=n-n%2
 return ans

def clock(a,b):
 k=1
 while pow(2,k,b)!=1:k+=1
 B=1<<k;A=a*(B-1)//b;R=int(format(A,f'0{k}b')[::-1],2)
 C=A-g(A+1)+1;K=a*R-b*C
 assert A%2==0 and B//2<=A<B-1 and K>0 and K%2==0
 for n in range(2,20):
  nn=B*n-R
  assert g(nn)==g(n)+A*n-C
  assert a*nn-b*g(nn)==a*n-b*g(n)-K
 return k,B,A,R,C,K
clocks=0
for b in range(3,102,2):
 for a in range(2,b,2):
  if math.gcd(a,b)==1 and 2*a>b:clock(a,b);clocks+=1

seed_checks=materialized=0
receipt=json.loads((P/'VALIDATION_NORMAL.json').read_text())
for row in receipt['compressed_exact_recipes']:
 a,b,r=row['a'],row['b'],row['target'];q=m.construct(a,b,r);L=q['seed_L'];D,E,U,d=q['seed_profile'];h=q['repetitions']
 assert E==d==1 and U==2*L-1 and D==g(U)==g(L)+L-L%2
 assert 2*a*U-b*D-q['clock']['K']*h==r
 seed_checks+=1
 k,B,A,R,C,K=clock(a,b)
 if h*k<=20000:
  power=pow(B,h);nn=power*(2*U)-R*((power-1)//(B-1))
  assert nn%2==1 and a*nn-b*g(nn)==r
  materialized+=1

def source(n):
 letters=[]
 while n>2:
  if n%2==0:letters.append('O');n//=2
  else:letters.append('Z');n=(n+1)//2
 assert n==2
 T=[2,0]
 for c in reversed(letters):
  size=len(T);A=T.copy()+[0]
  for i in range(size):
   if i==0 and c=='Z':A[size]|=1
   else:A[i]|=1<<size
  T=A
 return T

def inv(T):
 n=len(T);size=1<<n;P=[[0]*n for _ in range(size)];Q=[[0]*n for _ in range(size)]
 for i in range(n):P[1<<i][i]=1
 Q[1][0]=1
 for S in range(1,size):
  for u in range(n):
   if not P[S][u]:continue
   rest=T[u]&~S
   while rest:
    bit=rest&-rest;rest-=bit;v=bit.bit_length()-1
    P[S|bit][v]+=P[S][u];Q[S|bit][v]+=Q[S][u]
 H=[sum(row)for row in P];H[0]=1
 U=sum(H[(size-1)^S]*sum(Q[S])for S in range(1,size));C=sum(Q[-1][v]for v in range(n)if T[v]&1)
 return (H[-1],P[-1][0],U,H[size-2]),sum(Q[-1]),C
literal=0
for row in receipt['literal_clocks']:
 D,E,U,d=row['profile'];T=source(U);p,s,C=inv(T)
 assert p==(g(U),1,U,1)==(D,E,U,d)and s==C==1
 assert all(T[i]>>j&1 for i in range(1,len(T))for j in range(i+1,len(T)))
 literal+=1
for n in range(3,20000,2):
 r=3*n-4*g(n)
 assert r!=3 and (r<0 or r%4==1)
out={'status':'PASS','independent_clock_slopes_denominator_through101':clocks,'clock_identity_inputs_per_slope':18,'all_frozen_seed_recipes_checked_via_original_g':seed_checks,'full_exact_arithmetic_n_witnesses_materialized':materialized,'literal_strong_unit_sources':literal,'hostile_three_quarters_inputs':9999,'methods':'K via digit-complement g identity; direct g(B^h n-R sum B^j); graph recovered from cover integer'}
(HERE/'RATIONAL_INDEPENDENT.json').write_text(json.dumps(out,indent=2,sort_keys=True)+'\n');print(json.dumps(out,sort_keys=True))

# Independently construct the C2 offset pairs, including the nontrivial J root.
def birth(T,c,g):
 n=len(T);A=T.copy()+[0]
 for u in range(n):
  if g!='O'and u==c:A[n]|=1<<u
  else:A[u]|=1<<n
 return A,n if g=='J'else c

def graft(W,c,kind):
 T,a=birth([2,0],0,'J');entries={a,0}
 if kind=='M':T,u=birth(T,a,'J');entries={u,0}
 n=len(T);A=T.copy()+[0]*len(W)
 for i,row in enumerate(W):
  for j in range(len(W)):
   if i!=j and not row>>j&1:A[n+i]|=1<<(n+j)
 for i in range(n):
  for j in range(len(W)):
   if j==c and i in entries:A[n+j]|=1<<i
   else:A[i]|=1<<(n+j)
 return A
pairs=[]
for cover,parent,child,offset in [(11,157,313,-1),(13,187,375,1)]:
 q=source(cover);p,s,c=inv(q);assert p==(g(cover),1,cover,1)and s==c==1
 W,r=birth(q,0,'O');T=graft(W,r,'N');p,s,c=inv(T);assert p[0]==parent and c==2
 W,r=birth(q,0,'J');T=graft(W,r,'M');p,s,c=inv(T);assert p[0]==child and c==2
 assert child==2*parent+offset and parent%34==(19-2*offset)%34
 pairs.append([parent,child,offset])
out['independent_literal_C2_offset_pairs']=pairs
(HERE/'RATIONAL_INDEPENDENT.json').write_text(json.dumps(out,indent=2,sort_keys=True)+'\n');print('C2_OFFSET_CONTROLS_PASS',pairs)
