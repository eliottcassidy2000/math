from pathlib import Path
from math import gcd,isqrt
import json,itertools,re
C={}
def check(v,k):
 if not v:raise RuntimeError(k)
 C[k]=C.get(k,0)+1
def forms(d):
 D=-4*d;out=[]
 for a in range(1,isqrt((-D)//3)+1):
  for b in range(-a,a+1):
   if (b*b-D)%(4*a):continue
   c=(b*b-D)//(4*a)
   if c<a or gcd(gcd(a,b),c)!=1:continue
   if (abs(b)==a or a==c) and b<0:continue
   out.append((a,b,c))
 return out
def ambiguous(f):a,b,c=f;return b==0 or abs(b)==a or a==c
def jacobi(a,n):
 if n<=0 or n%2==0:raise ValueError(n)
 a%=n;s=1
 while a:
  while a%2==0:
   a//=2
   if n%8 in (3,5):s=-s
  a,n=n,a
  if a%4==n%4==3:s=-s
  a%=n
 return s if n==1 else 0
expected={63:[(1,0,63),(7,0,9),(8,-6,9),(8,6,9)],273:[(1,0,273),(2,2,137),(3,0,91),(6,6,47),(7,0,39),(13,0,21),(14,14,23),(17,8,17)]}
rows={};codes={}
for d in (63,273):
 fs=forms(d);check(fs==expected[d],'literal_complete_forms');rr=[]
 for a,b,c in fs:
  check(b*b-4*a*c==-4*d and gcd(gcd(a,b),c)==1,'form_valid')
  ns=sorted({a*x*x+b*x*y+c*y*y for x in range(-16,17) for y in range(-16,17) if (x or y) and gcd(a*x*x+b*x*y+c*y*y,4*d)==1})
  ds=(-3,) if d==63 else (-4,-3,-7)
  chars=tuple(jacobi(delta,ns[0]) for delta in ds)
  for n in ns:check(tuple(jacobi(delta,n) for delta in ds)==chars,'representation_character_controls')
  rr.append({'form':[a,b,c],'representative':ns[0],'characters':chars,'self_inverse':ambiguous((a,b,c))})
  if d==273:codes[a,b,c]=tuple(int(v==-1) for v in chars)
 rows[d]=rr
check(len(set(codes.values()))==8 and all(ambiguous(f) for f in forms(273)),'lossless_genus_273')
check(sum(ambiguous(f) for f in forms(63))==2 and len(forms(63))==4,'cyclic_four_63')
check(-252==6**2*(-7),'conductor_63')
# Fundamental discriminant -7: general reduction checked directly.
f7=[]
for a in range(1,isqrt(7//3)+1):
 for b in range(-a,a+1):
  if (b*b+7)%(4*a):continue
  c=(b*b+7)//(4*a)
  if c>=a and gcd(gcd(a,b),c)==1 and not((abs(b)==a or a==c) and b<0):f7.append((a,b,c))
check(f7==[(1,1,2)],'maximal_order_hostile')
def mask(bits):return (33 if bits[0] else 0)^(18 if bits[1] else 0)^(12 if bits[2] else 0)
es=list(itertools.combinations(range(4),2));ix={e:i for i,e in enumerate(es)}
def rev(g):
 return sum(1<<ix[tuple(sorted((3-a,3-b)))] for i,(a,b) in enumerate(es) if g>>i&1)
for u,v in itertools.product(codes.values(),repeat=2):
 check(mask(tuple(a^b for a,b in zip(u,v)))==mask(u)^mask(v),'graph_group_code')
for u in codes.values():
 check(rev(mask(u))==mask(u),'inversion_reflection')
 check(mask(tuple(1-a for a in u))==(mask(u)^63),'translation_complement')
check(codes[(6,6,47)]==(1,1,1),'named_eta')
regular=[]
for g in range(64):
 ds=[sum(bool(g>>ix[tuple(sorted((a,b)))]&1) for b in range(4) if b!=a) for a in range(4)]
 if len(set(ds))==1:regular.append(g)
check(set(regular)=={mask(b) for b in codes.values()},'all_regular_four_graphs')
# Actual fixed-path tournament arc dictionary, no imported producer.
def beat(g,a,b):
 if a>b:return not beat(g,b,a)
 return b==a+1 or not(g>>ix[a,b-1]&1)
for g in map(mask,codes.values()):
 check(all(beat(g,i,i+1) for i in range(4)),'hamiltonian_path')
 for a,b in itertools.combinations(range(5),2):
  check(beat(g,a,b)==beat(g,4-b,4-a),'actual_blue_path')
  check(beat(g^63,a,b)==(beat(g,a,b) if b==a+1 else not beat(g,a,b)),'actual_path_complement')
check((63%15,63%35)==(273%15,273%35)==(3,28),'joint_register_hostile')
check(all(s%2 and s%8!=5 for s in (63,273)),'primitive_source_membership')
def U(n):
 y=3*n+1;a=0
 while y%2==0:y//=2;a+=1
 return y,a
def P(n):
 a=2 if n%8==1 else 1;y=(3*n+1)//2**a;clock=a
 while y%8==5:y=(y-1)//4;clock+=2
 return y,clock
check(P(63)==(95,1) and P(273)==(51,4),'literal_primitive_step')
check(U(273)==(205,2),'clock_type_boundary')
check([U(n)[0] for n in(3,7,17)]==[5,11,13],'actual_odd_hostiles')
for d,want in [(3,True),(5,True),(7,True),(11,False),(13,True),(17,False)]:
 check(all(ambiguous(f) for f in forms(d))==want,'idoneal_small_flag')
known=[]
for d in range(1,1849):
 fs=forms(d);check(bool(fs),'reduction_nonempty')
 if all(ambiguous(f) for f in fs):known.append(d)
check(len(known)==65 and known[-1]==1848,'bounded_65_census')
# Independently parse every printed 273 form/representative/character row.
text=Path(__file__).with_name('IDONEAL_INFORMATION_CHANNELS.md').read_text()
parsed=[]
for a,b,c,n,s1,s2,s3 in re.findall(r'^\| \((\d+),(-?\d+),(\d+)\) \|(\d+)\| \(([+-]),([+-]),([+-])\) \|$',text,re.M):
 f=tuple(map(int,(a,b,c)));n=int(n);ss=tuple(1 if z=='+' else -1 for z in(s1,s2,s3))
 check(f in codes and tuple(1-2*z for z in codes[f])==ss,'printed_character_row')
 check(gcd(n,1092)==1 and any(gcd(x,y)==1 and f[0]*x*x+f[1]*x*y+f[2]*y*y==n for x in range(-16,17) for y in range(-16,17)),'printed_representative')
 parsed.append(f)
check(len(parsed)==8,'printed_table_complete')
print(json.dumps({'status':'PASS','checks':C,'total_checks':sum(C.values()),'class_rows':rows,'graph_code_rows':[{'form':f,'bits':b,'graph_code':mask(b)} for f,b in codes.items()],'known_up_to_1848':known,'scope':'Finite reduction checks; genus homomorphism uses cited theory; no global completeness or Collatz inference.'},indent=2,sort_keys=True))
