"""Literal seed reconstruction, finite closure inverse audit, exact digital controls."""
from pathlib import Path
from itertools import permutations,product,combinations
from fractions import Fraction
from collections import deque
import json,sys,hashlib
SOURCE=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(SOURCE));import digital as d
raw=json.loads((SOURCE/'actual_seed_roots.json').read_text())
roots=[(x['state'][0],x['state'][1],x['state'][2],x['state'][0]-x['state'][3]) for x in raw]
def graph(n,code):
 A=[[False]*n for _ in range(n)]
 for i,(u,v) in enumerate(combinations(range(n),2)):
  A[u][v]=bool(code>>i&1);A[v][u]=not A[u][v]
 return A

def literal(A):
 n=len(A);h=0;start=[0]*n;end=[0]*n;W=[[0]*n for _ in range(n)]
 for p in permutations(range(n)):
  bad=[i for i in range(n-1) if not A[p[i]][p[i+1]]]
  if not bad:h+=1;start[p[0]]+=1;end[p[-1]]+=1
  for split in range(1,n):
   if all(i==split-1 for i in bad):W[p[split-1]][p[split]]+=1
 return h,start,end,W

def profile(A,cut):
 h,s,t,W=literal(A);n=len(A)
 q=sum(s[v]+sum(W[u][v] for u in range(n)) for v in range(n) if cut>>v&1)
 F=sum(s[v] for v in range(n) if cut>>v&1)+sum(t[u] for u in range(n) if not cut>>u&1)+sum(W[u][v] for u in range(n) for v in range(n) if not cut>>u&1 and cut>>v&1)
 g=sum(t[u] for u in range(n) if cut>>u&1)
 return h,F,q,g

def actual(seed,blocks):
 n,code,cut=seed;A=graph(n,code);xs=[bool(cut>>i&1) for i in range(n)]
 for b in blocks:
  m=len(A);size=3 if isinstance(b,str) else 1
  for row in A:row.extend([True]*size)
  A.extend([[False]*(m+size) for _ in range(size)])
  if size==3:
   for i in range(3):A[m+i][m+(i+1)%3]=True
   xs.extend(i<int(b[1]) for i in range(3))
  else:xs.append(bool(b))
 return [row+[not xs[i]] for i,row in enumerate(A)]+[xs+[False]]

def paths(A):
 n=len(A);out=[sum(1<<v for v in range(n) if A[u][v]) for u in range(n)];D=[[0]*n for _ in range(1<<n)]
 for u in range(n):D[1<<u][u]=1
 for mask in range(1,1<<n):
  for u,x in enumerate(D[mask]):
   if x:
    todo=out[u]&~mask
    while todo:
     bit=todo&-todo;todo-=bit;D[mask|bit][bit.bit_length()-1]+=x
 return sum(D[-1])

def step(s,b):
 h,F,q,g=s
 if b==0:return h,F+q+g,2*q,0
 if b==1:return h,F,2*q+h,h
 j=int(b[1]);return 3*h,3*F+(9-3*j)*q+(3-j)*g+(2*h if j in (1,2) else 0),12*q+3*j*h,j*h

def run(s,w):
 for b in w:s=step(s,b)
 return s
allprofiles=set();labelled=0
for n in range(1,6):
 for code in range(1<<(n*(n-1)//2)):
  A=graph(n,code);h,s,t,W=literal(A)
  for cut in range(1<<n):
   q=sum(s[v]+sum(W[u][v] for u in range(n)) for v in range(n) if cut>>v&1)
   F=sum(s[v] for v in range(n) if cut>>v&1)+sum(t[u] for u in range(n) if not cut>>u&1)+sum(W[u][v] for u in range(n) for v in range(n) if not cut>>u&1 and cut>>v&1)
   g=sum(t[u] for u in range(n) if cut>>u&1);allprofiles.add((h,F,q,g));labelled+=1
if len(roots)!=273 or len(set(roots))!=273 or set(roots)!=allprofiles:raise RuntimeError('actual seed bank')
actual_checks=0
for rec,s in zip(raw,roots):
 n,code,cut=rec['witness']
 if profile(graph(n,code),cut)!=s or s[0]%2!=1:raise RuntimeError('seed witness')
 if s[2]==0 and not(s[1]==s[0] and s[3]==0 and cut==0):raise RuntimeError('zero-q realization')
 for w in ((),(0,),(1,),(0,1),(1,0),('C0',),('C1',),('C2',),('C3',)):
  if paths(actual(rec['witness'],w))!=run(s,w)[1]:raise RuntimeError(('actual adjacency',s,w))
  actual_checks+=1
formula=0
for s in roots:
 for length in range(7):
  for w in product((0,1),repeat=length):
   final=run(s,w);B=sum((1-b)*2**i for i,b in enumerate(w))
   if d.output(s,w)!=final[1]:raise RuntimeError('digital formula')
   if w:
    first=d.kappa_sources(s)[w[0]]
    if final[1]-d.phi(s[0],final[2])!=first[1]:raise RuntimeError('potential')
   formula+=1
# Exact rational identity and the actual run-count error bound.
for n in range(4096):
 bits=[(n>>i)&1 for i in range(n.bit_length())]
 f=1+sum((n>>(j+1))+((n>>(j+1))&1) for j in range(max(0,len(bits)-1)) if not bits[j])
 rho=sum((Fraction(bits[i],2**(i+1)) for i in range(len(bits))),Fraction())
 error=Fraction(f)-n*(1-rho)
 r=sum(z and (i==0 or not bits[i-1]) for i,z in enumerate(bits))
 if f!=d.f(n) or f!=d.f_pair_formula(n) or error!=d.error_runs(n):raise RuntimeError('exact digital/error')
 if n and not 0<error<=Fraction(1-(n&1))+Fraction(3*r,4):raise RuntimeError('run error bound')
# Independent unbounded-depth closure, using only global safe final cap4M.
# q is nondecreasing under every block; any target<=M has a reduced final q<4M.
M=63;inverse_checks=0;closure_states=0;returned_dp=0
for rec,s in zip(raw,roots):
 found={s[1]};seen={s};queue=deque([s])
 while queue:
  v=queue.popleft()
  if v[1]>M or v[2]>=4*M:continue
  for b in (0,1,'C0','C1','C2','C3'):
   z=step(v,b)
   if z[1]>M or z[2]>=4*M or z in seen:continue
   seen.add(z);queue.append(z);found.add(z[1])
 closure_states+=len(seen)
 for N in range(1,M+1,2):
  if N%3==0:continue
  answer=d.inverse_non3(N,s)
  if (answer is not None)!=(N in found):raise RuntimeError(('inverse vs forward closure',s,N,answer,sorted(found)))
  if answer is not None:
   if run(s,answer)[1]!=N:raise RuntimeError('inverse reconstruction')
   order=rec['witness'][0]+1+sum(3 if isinstance(b,str) else 1 for b in answer)
   if order<=9 and rec in raw[:12]:
    if paths(actual(rec['witness'],answer))!=N:raise RuntimeError('inverse actual witness')
    returned_dp+=1
  inverse_checks+=1
# Full-length scalar negative control: odd-address stripping forces finite G interval.
G=[d.f(n)+n+n%2 for n in range(6,13)]
if 13 in G or d.inverse_addresses(13,(1,1,0,0)) or d.inverse_addresses(13,(1,1,1,1)):raise RuntimeError('13 singleton control')
# Leading-zero exception has an actual indispensable witness.
s=(1,1,0,0);final=run(s,('C1',));Q=final[2]//3
if d.ancestor_word(1,0,Q)!=(1,) or ('C1',) not in d.macro_decodings(s,Q,1):raise RuntimeError('leading-zero exception')
# Frontier-miss repairs: inspect known Q only, never scan up to each large target.
repairs=[]
file=SOURCE/'frontier_repairs.jsonl'
for line in file.read_text().splitlines():
 rec=json.loads(line);s=profile(graph(rec['seed'][0],rec['seed'][1]),rec['seed'][2]);w=tuple(int(b) if b in ('0','1') else b for b in rec['blocks']);final=run(s,w);r=sum(isinstance(b,str) for b in w)
 if final[1]!=rec['target'] or final[2]>=4*final[1] or final[2]%(3**r):raise RuntimeError('repair profile')
 Q=final[2]//3**r;decoded=d.macro_decodings(s,Q,r)
 if not any(run(s,v)==final for v in decoded):raise RuntimeError('repair inverse')
 repairs.append(dict(target=rec['target'],q=final[2],normalized_Q=Q,triangles=r,actual_order=len(actual(rec['seed'],w))))
print(json.dumps(dict(status='PASS',labelled_seed_cut_profiles=labelled,actual_seed_states=273,literal_adjacency_checks=actual_checks,singleton_formula_checks=formula,rational_digit_checks=4096,forward_closure_target_cap=M,forward_closure_states=closure_states,inverse_against_closure_checks=inverse_checks,inverse_actual_dp=returned_dp,negative13_G_values=G,frontier_repairs=repairs,source_sha256=hashlib.sha256((SOURCE/'digital.py').read_bytes()).hexdigest(),seed_sha256=hashlib.sha256((SOURCE/'actual_seed_roots.json').read_bytes()).hexdigest()),sort_keys=True))
