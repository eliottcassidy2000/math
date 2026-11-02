"""Hostile independent Session136 checks. Imports own kernel only."""
from kernel import *
from itertools import combinations
from fractions import Fraction as Fr
import json,hashlib,pathlib

out={}
counts=[]
for n in range(6):
 for mask in range(1<<(n*(n-1)//2)):
  a=tournament(n,mask);validate(a)
  require(h_literal(a)==h_recursive(a)==ocf(a),('kernel',n,mask))
 counts.append(1<<(n*(n-1)//2))
out['kernel_all_labelled_n0_to_n5']=counts

a=tournament(6,83);values={s:h_literal(extend(a,s)) for s in range(1,63)}
require(h_literal(a)==17 and min(values.values())==27, 'ear hostile H')
require([s for s,h in values.items() if h==27]==[14], 'ear hostile argmin')
out['ear_balanced_minimum']={'rows':a,'H':17,'min_child':27,'min_cuts':[14],'singleton_complement_children':{s:h for s,h in values.items() if s.bit_count() in(1,5)}}

pairs=[]
for a,env in [((48,1,3,7,14,30),11),((34,8,3,5,15,30),13)]:
 validate(a);ds=[h_literal(induced(a,[u for u in range(6) if u!=v]))for v in range(6)]
 require(h_literal(a)==23 and max(ds)==env and len(scc(a))==1,'obstruction collision')
 pairs.append({'rows':a,'H':23,'deletions':ds,'envelope':env})
out['obstruction_collision']=pairs

cyclic=[]
for n,expected,expected_rho in [(9,3267,105),(11,93027,2125)]:
 a=tuple(sum(1<<((i+j)%n) for j in range(1,(n+1)//2))for i in range(n));validate(a)
 h=h_recursive(a);rho=min(h_recursive(induced(a,[v for v in range(n) if v not in pair]))for pair in combinations(range(n),2))
 require((h,rho)==(expected,expected_rho),'cyclic control')
 cyclic.append({'n':n,'H':h,'rho':rho,'certified_bottom_band_holes':max(0,(rho+1)//2-(4*n-6))})
out['cyclic_ear_hole_controls']=cyclic

# Full OCF even-path cut kernel independent of boundary/endpoint producer.
def cut_kernel(a):
 n=len(a);r=[[0]*n for _ in a]
 for k in range(2,n+1,2):
  for p in permutations(range(n),k):
   if all(a[u]>>v&1 for u,v in zip(p,p[1:])):
    r[p[0]][p[-1]]+=h_recursive(induced(a,[u for u in range(n) if u not in p]))
 return r

ear_count=0
for n in range(1,5):
 for code in range(1<<(n*(n-1)//2)):
  a=tournament(n,code);r=cut_kernel(a);h=h_literal(a)
  f=[]
  for s in range(1<<n):
   expected=h+2*sum(r[i][j]for i in range(n) for j in range(n)if s>>i&1 and not s>>j&1)
   require(h_literal(extend(a,s))==expected,('cut-kernel',n,code,s));f.append(expected);ear_count+=1
  for i,j in combinations(range(n),2):
   require(r[i][j]+r[j][i]>=h_recursive(induced(a,[v for v in range(n)if v not in(i,j)])), 'positive pair floor')
  for s in range(1<<n):
   for t in range(1<<n):
    rhs=2*sum(r[i][j]+r[j][i]for i in range(n)if s>>i&1 and not t>>i&1 for j in range(n)if t>>j&1 and not s>>j&1)
    require(f[s]+f[t]-f[s&t]-f[s|t]==rhs,'uncrossing')
out['literal_ear_cut_and_uncrossing_count']=ear_count

# One-hole quotient words, computed literally; do not import substitution law.
def context_coefficients(q,x,maxc):
 m=len(q)-1;other=tuple(i for i in range(len(q))if i!=x);cs=[]
 for c in range(1,maxc+1):
  count=0
  for positions in combinations(range(m+c),c):
   pos=set(positions)
   for outside in permutations(other):
    it=iter(outside);w=tuple(x if j in pos else next(it)for j in range(m+c))
    if all(q[u]>>v&1 for u,v in zip(w,w[1:])):count+=1
  cs.append(factorial(c)*count)
 return cs

blocks=[tournament(1,0),tournament(2,0)]+[tournament(3,c) for c in range(8)]
glue_count=0
for qn in range(1,5):
 for code in range(1<<(qn*(qn-1)//2)):
  q=tournament(qn,code)
  for x in range(qn):
   cs=context_coefficients(q,x,4)
   require(all(v==0 for c,v in enumerate(cs,1)if 2*(c-1)>qn-1),'context return depth')
   interacting=any(x in comp and len(comp)>1 for comp in scc(q))
   require((cs[1]>0)==interacting,'interaction iff a2')
   for b in blocks:
    pc=path_cover_literal(b);expected=sum(cs[c-1]*pc[c]for c in range(1,len(pc)))
    actual=h_recursive(substitute(q,[b if i==x else tournament(1,0)for i in range(qn)]))
    require(actual==expected,('context law',qn,code,x,b))
    if interacting:require(actual>=2**len(b)+1,'module barrier')
    glue_count+=1
for s in range(1,5):
 # C3[T,Tr_s,Tr_s] written as a one-hole quotient of order2s+1
 q=substitute(tournament(3,5),[tournament(1,0),tournament(s,(1<<(s*(s-1)//2))-1),tournament(s,(1<<(s*(s-1)//2))-1)])
 # direct return coefficient can be huge to enumerate at s4; exploit only triangular top via words at s<=2 here
 if s<=2:
  cs=context_coefficients(q,0,s+2)
  require(cs[s]==factorial(s+1)*factorial(s)**2 and cs[s+1]==0,'top detector coefficient')
out['literal_context_cases']=glue_count

# Exact arithmetic certificate for proposed universal witness compression.
for n,f in [(23,140625),(24,234375),(25,390625)]:require(f**22>=80267**n,('compression base',n))
require(5**22>=80267**3 and 285**2>80267 and 283**2<80265,'compression step')
out['compression_integer_checks']='PASS; residual first-order table audited separately'

# Affine odd-modulus operations and dyadic initialization, checked independently.
def A(p,m):d,t=p;return(d%m,2*t%m)
def B(p,m):d,t=p;return((d-t)%m,(4*t+1)%m)
for n in range(1,100,2):
 # Breadth-first orbit under legal positive actual blocks from B(seed).
 start=B((0,0),n);todo=[start];seen={start}
 for p in todo:
  for q in(A(p,n),B(p,n)):
   if q not in seen:seen.add(q);todo.append(q)
 require({(d+t)%n for d,t in seen}==set(range(n)),('odd modular coverage',n))
for a in range(9):
 mod=1<<a;lo=max(1,(a+1)//2)
 ds={((3*r+1-4**r)//9)%mod for r in range(lo,lo+mod)}
 require(ds==set(range(mod)),('dyadic initialization',a))
out['positive_affine_odd_orbits_moduli']='every odd modulus1..99'
out['dyadic_initialization']='all2^a,a0..8'

path=pathlib.Path(__file__).parent
out['own_sources_sha256']={p.name:hashlib.sha256(p.read_bytes()).hexdigest()for p in(path/'kernel.py',path/'audit.py')}
print(json.dumps(out,indent=2,sort_keys=True))
