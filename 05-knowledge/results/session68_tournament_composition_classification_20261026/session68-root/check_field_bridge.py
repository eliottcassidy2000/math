from itertools import product
import json

def need(x,s):
 if not x:raise RuntimeError(s)
def mul(a,b,f):
 n=f.bit_length()-1;r=0
 while b:
  if b&1:r^=a
  b>>=1;a<<=1
  if a>>n:a^=f
 return r
def powf(a,k,f):
 r=1
 while k:
  if k&1:r=mul(r,a,f)
  k>>=1;a=mul(a,a,f)
 return r
def rem(a,b):
 while a and a.bit_length()>=b.bit_length():a^=b<<(a.bit_length()-b.bit_length())
 return a
def gcd(a,b):
 while b:a,b=b,rem(a,b)
 return a
def irred(f):
 n=f.bit_length()-1;x=2;t=x
 for k in range(1,n+1):
  t=mul(t,t,f)
  if k<=n//2 and gcd(t^x,f)!=1:return False
 return t==x
def rank(vs):
 piv={}
 for x in vs:
  while x:
   b=x.bit_length()-1
   if b in piv:x^=piv[b]
   else:piv[b]=x;break
 return len(piv)
def normal(f):
 n=f.bit_length()-1
 for a in range(1,1<<n):
  v=[a]
  for _ in range(n-1):v.append(mul(v[-1],v[-1],f))
  if rank(v)==n:return a,v
 raise RuntimeError('normal missing')
def field(mask,bas):
 r=0
 for i,b in enumerate(bas):
  if mask>>i&1:r^=b
 return r
def degree(x,f):
 y=mul(x,x,f);d=1
 while y!=x:y=mul(y,y,f);d+=1
 return d
def act(mask,p):return sum(((mask>>i)&1)<<p[i] for i in range(len(p)))
def orbit(mask,n):return {((mask<<j)|(mask>>(n-j)))&((1<<n)-1) for j in range(n)}

f=11;need(irred(f),'F8 modulus');basis=[3,5,7];need(rank(basis)==3,'normalbasis rank')
need(all(mul(basis[i],basis[i],f)==basis[(i+1)%3] for i in range(3)),'Frobenius basis')
need(rank([2,4,6])==2,'bad tracezero basis hostile')
map8={mask:field(mask,basis) for mask in range(8)};inv={v:k for k,v in map8.items()}
J=lambda x:map8[act(inv[x],[0,2,1])]
for mask,x in map8.items():
 need(map8[act(mask,[1,2,0])]==mul(x,x,f),'rotation Frobenius')
 need(degree(x,f)==len(orbit(mask,3)),'field degree word period')
 need(J(mul(x,x,f))==powf(J(x),4,f),'reflection conjugates Frobenius')
 need(1^J(x)!=x,'complemented reflection no fixed')
need(J(mul(2,2,f))!=mul(J(2),J(2),f),'reflection is not field automorphism')

f9=next(v for v in range(513,1024,2) if irred(v));alpha,basis9=normal(f9)
W=[]
for rot in range(3):
 for shifts in product(range(3),repeat=3):W.append([3*((b+rot)%3)+(j+shifts[b])%3 for b in range(3) for j in range(3)])
need(len({tuple(p) for p in W})==81,'W2 action size')
cyc=[0,3,6,1,4,7,2,5,8];g=[0]*9
for i,u in enumerate(cyc):g[u]=cyc[(i+1)%9]
need(g in W,'odometer belongs to wreath action')
witness=[]
for chars in ['000000001','001001011']:
 mask=sum(int(c)<<i for i,c in enumerate(chars));word=sum(((mask>>u)&1)<<i for i,u in enumerate(cyc));x=field(word,basis9)
 stab=sum(act(mask,p)==mask for p in W)
 need(degree(x,f9)==9,'same field degree 9')
 witness.append({'colors_in_index_order':chars,'normal_basis_word':word,'field_element':x,'field_degree':9,'W2_stabilizer':stab})
need([r['W2_stabilizer'] for r in witness]==[9,1],'field degree loses hierarchical symmetry')

# Relative phases of two fixed rotation orbits, not independent free labels.
crt=[]
for m,d,e in [(3,3,3),(15,3,5)]:
 a=sum(1<<i for i in range(0,m,d));b=sum(1<<i for i in range(0,m,e));oa=orbit(a,m);ob=orbit(b,m);pairs={(x,y) for x in oa for y in ob};num=0
 while pairs:
  x,y=min(pairs);component={((((x<<j)|(x>>(m-j)))&((1<<m)-1)),(((y<<j)|(y>>(m-j)))&((1<<m)-1))) for j in range(m)};pairs-=component;num+=1
 from math import gcd as igcd,lcm
 need(num==igcd(d,e),'gcd relative phase');need(len(component)==lcm(d,e),'lcm joint orbit')
 crt.append({'m':m,'periods':[d,e],'relative_phase_orbits':num,'joint_orbit_size':len(component)})
print(json.dumps({'status':'PASS','F8_modulus_bits':f,'normal_basis_F8':basis,'F8_map':map8,'reflection_nonmultiplicative_witness':[2,2,J(mul(2,2,f)),mul(J(2),J(2),f)],'F512_modulus_bits':f9,'normal_element_F512':alpha,'normal_basis_F512':basis9,'hierarchy_witness':witness,'CRT':crt},sort_keys=True,indent=2))
