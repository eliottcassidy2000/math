"""Independent reduced forms, Dirichlet composition, prime characters and edge sets."""
from pathlib import Path
from math import gcd,isqrt
from itertools import product,combinations
from collections import Counter
import json,re,hashlib,argparse
P=argparse.ArgumentParser();P.add_argument('--proof',type=Path);args=P.parse_args();base=Path(__file__).resolve().parent
proof=args.proof or base.parent/'IDONEAL_INFORMATION_CHANNELS.md'
if not proof.exists():proof=base.parent/'session18-idoneal'/'IDONEAL_INFORMATION_CHANNELS.md'
text=proof.read_text();checks=Counter()
def need(ok,key):
 if not ok:raise RuntimeError(key)
 checks[key]+=1
def forms(D):
 out=set();bound=isqrt((-D)//3)
 for b in range(-bound,bound+1):
  n=b*b-D
  if n%4:continue
  n//=4
  for a in range(max(1,abs(b)),bound+1):
   if n%a:continue
   c=n//a
   if a>c or gcd(gcd(a,abs(b)),c)>1:continue
   if b<0 and (abs(b)==a or a==c):continue
   out.add((a,b,c))
 return sorted(out)
def reduce(f):
 a,b,c=f
 for _ in range(1000):
  q=(a-b)//(2*a);a,b,c=a,b+2*a*q,a*q*q+b*q+c
  if a>c:a,b,c=c,-b,a;continue
  if a==c and b<0:b=-b
  need(-a<b<=a and a<=c and not(a==c and b<0),'reduction_output')
  return(a,b,c)
 raise RuntimeError('reduction did not terminate')
def xgcd(a,b):
 r0,r1=a,b;s0,s1=1,0;t0,t1=0,1
 while r1:q=r0//r1;r0,r1=r1,r0-q*r1;s0,s1=s1,s0-q*s1;t0,t1=t1,t0-q*t1
 if r0<0:return-r0,-s0,-t0
 return r0,s0,t0
def transformed(f,D,coprime_to=1):
 a,b,c=f;options=[]
 for u,v in product(range(-8,9),repeat=2):
  if gcd(u,v)!=1:continue
  A=a*u*u+b*u*v+c*v*v
  if gcd(A,abs(D)*coprime_to)!=1:continue
  options.append((A,u,v))
 need(bool(options),'coprime_equivalent_exists');A,u,v=min(options);_,s,t=xgcd(u,v);r=-t;B=2*a*u*r+b*(u*s+v*r)+2*c*v*s;C=a*r*r+b*r*s+c*s*s
 need(u*s-v*r==1 and B*B-4*A*C==D and reduce((A,B,C))==f,'proper_equivalence');return A,B,C
def compose(f,g,D):
 a,b,c=transformed(f,D);aa,bb,cc=transformed(g,D,a)
 need(b%2==bb%2==0 and gcd(a,aa)==1,'coprime_composition_domain')
 h=b//2+a*((bb//2-b//2)*pow(a,-1,aa)%aa);B=2*h;A=a*aa
 need((B*B-D)%(4*A)==0,'composite_integrality');return reduce((A,B,(B*B-D)//(4*A)))
def chars(n):
 need(n>0 and gcd(n,1092)==1,'genus_character_domain')
 return(1 if n%4==1 else -1,1 if n%3==1 else -1,1 if n%7 in(1,2,4) else -1)
expected63=[(1,0,63),(7,0,9),(8,-6,9),(8,6,9)];need(forms(-252)==expected63,'literal_63_forms');need(forms(-7)==[(1,1,2)],'maximal_order_minus7');need(-252==36*(-7),'conductor6')
parsed=[]
for a,b,c,n,x,y,z in re.findall(r'^\| \((\d+),(-?\d+),(\d+)\) \|(\d+)\| \(([+-]),([+-]),([+-])\) \|$',text,re.M):
 f=tuple(map(int,(a,b,c)));n=int(n);ch=tuple(1 if v=='+' else -1 for v in(x,y,z));need(chars(n)==ch,'printed_character_values');need(any(gcd(u,v)==1 and f[0]*u*u+f[1]*u*v+f[2]*v*v==n for u,v in product(range(-20,21),repeat=2)),'printed_primitive_representative');parsed.append((f,n,ch))
F=forms(-1092);need(sorted(f for f,n,ch in parsed)==F and len(parsed)==8,'printed_complete_forms');code={f:tuple(int(c==-1) for c in ch) for f,n,ch in parsed};need(len(set(code.values()))==8,'genus_code_bijection')
multiplication={}
for D in(-252,-1092):
 fs=forms(D);identity=fs[0];tab=[]
 for f,g in product(fs,repeat=2):
  h=compose(f,g,D);need(h in fs,'composition_closed');tab.append([fs.index(f),fs.index(g),fs.index(h)])
  if D==-1092:need(code[h]==tuple(a^b for a,b in zip(code[f],code[g])),'actual_class_composition_XOR')
 for f in fs:need(compose(f,reduce((f[0],-f[1],f[2])),D)==identity,'actual_inverse_composition')
 multiplication[str(D)]={'forms':fs,'table':tab}
need(compose(expected63[2],expected63[2],-252)==(7,0,9),'order4_square_nontrivial')
need(all(reduce((a,-b,c))==(a,b,c) for a,b,c in F),'all273_self_inverse')
for f,n,ch in parsed:
 for u,v in product(range(-12,13),repeat=2):
  if gcd(u,v)!=1:continue
  n=f[0]*u*u+f[1]*u*v+f[2]*v*v
  if gcd(n,1092)==1:need(chars(n)==ch,'independent_representation_characters')
matchings=[frozenset(((0,1),(2,3))),frozenset(((0,2),(1,3))),frozenset(((0,3),(1,2)))];E=frozenset(combinations(range(4),2));G=lambda bits:frozenset(e for i,b in enumerate(bits) if b for e in matchings[i]);rho=lambda g:frozenset(tuple(sorted((3-a,3-b))) for a,b in g)
regular=[]
for bits in product(range(2),repeat=6):
 g=frozenset(e for e,b in zip(sorted(E),bits) if b);degrees=[sum(i in e for e in g) for i in range(4)]
 if len(set(degrees))==1:regular.append(g)
need(set(regular)=={G(b) for b in code.values()},'whole_regular_graph_family');eta=(6,6,47);need(code[eta]==(1,1,1),'actual_eta')
for f in F:
 g=G(code[f]);need(rho(g)==g,'graph_reflection_identity');need(G(code[compose(eta,f,-1092)])==E-g,'actual_eta_complement');t=frozenset((j,i) if j>i+1 and (i,j-1) in g else(i,j) for i,j in combinations(range(5),2));need(all((i,i+1) in t for i in range(4)),'marked_path');need(frozenset((4-b,4-a) for a,b in t)==t,'actual_blue_tournament')
def U(n):
 t=3*n+1;a=(t&-t).bit_length()-1;return t>>a,a
def P(n):
 a=2 if n%8==1 else 1;t=(3*n+1)>>a
 while t%8==5:t=(t-1)//4;a+=2
 return t,a
need((63%15,63%35)==(273%15,273%35)==(3,28),'joint_register_equal');need(P(63)==(95,1) and P(273)==(51,4),'unstopped_P_clocks');need(U(273)==(205,2),'odd_U_distinct_clock');need([U(n)[0] for n in(3,7,17)]==[5,11,13],'U_hostiles')
known=[]
for d in range(1,1849):
 fs=forms(-4*d);need(bool(fs),'bounded_nonempty');flag=all(reduce((a,-b,c))==(a,b,c) for a,b,c in fs)
 if flag:known.append(d)
need(len(known)==65 and known[-1]==1848,'bounded65_only')
need(all(d in known for d in(3,5,7,13,273)) and all(d not in known for d in(11,17,63)),'idoneal_flags')
print(json.dumps({'status':'PASS','method':'Independent b-first reduction, proper SL2 transforms, coprime Dirichlet composition, modular prime characters, literal table and set-edge graphs; no producer import','proof_sha256':hashlib.sha256(proof.read_bytes()).hexdigest(),'total_checks':sum(checks.values()),'checks':dict(sorted(checks.items())),'class_composition':multiplication,'printed_table':parsed,'bounded_idoneal_candidates':known,'genus_homomorphism_scope':'CITED classical theory; direct finite composition independently agrees','global_completeness_audited':False},indent=2,sort_keys=True))
