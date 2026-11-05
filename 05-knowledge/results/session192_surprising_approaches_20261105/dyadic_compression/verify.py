#!/usr/bin/env python3
"""Sharp insulated-switch length; exact carry ledger; finite nonlinear banks.
All gates are active under python -O. Python standard library only.
"""
from fractions import Fraction
import hashlib,itertools,json,math
from pathlib import Path
HERE=Path(__file__).resolve().parent
CHECKS={}
def require(ok,msg,key):
 if not ok: raise RuntimeError(msg)
 CHECKS[key]=CHECKS.get(key,0)+1

def bits(A,n):return tuple((A>>j)&1 for j in range(n-1,-1,-1))
def integer(w):
 A=0
 for b in w:A=2*A+b
 return A

def count(w):
 h,A,p=1,0,0
 for b in w:
  if b==0:h+=A+p
  A=2*A+b;p=b
 return h

def slot(w):
 if not w:return 1
 return w[0]+1-w[-1]+sum(a==0 and b==1 for a,b in zip(w,w[1:]))+sum(1<<(i-j-1) for j in range(len(w)) for i in range(j+1,len(w)) if w[j] and not w[i])
def dp(w):
 n=len(w)+1;rows=[0]*n
 for i,b in enumerate(w,1):
  if b:rows[0]|=1<<i
  else:rows[i]|=1
  for j in range(i+1,n):rows[i]|=1<<j
 D=[[0]*n for _ in range(1<<n)]
 for i in range(n):D[1<<i][i]=1
 for mask in range(1,1<<n):
  for v in range(n):
   if D[mask][v]:
    avail=rows[v]&~mask
    while avail:
     bit=avail&-avail;avail-=bit
     D[mask|bit][bit.bit_length()-1]+=D[mask][v]
 return sum(D[-1])

def pair(e):return (1,0) if e else (0,1)
def spacers(k):
 if k==1:return []
 D=[(0,)*k]
 for j in range(1,k):
  q=((1<<(j-1))-1-integer(D[-1])//2)%(1<<(k-1))
  D.append((1,)+tuple(1-((q>>t)&1) for t in range(k-1)))
 return D

def bank(k,e):
 if k==1:return ()
 D=spacers(k);w=D[0]
 for j in range(1,k):w+=pair((e>>(j-1))&1)+D[j]
 return w

def length(k):return 0 if k==1 else k*k+2*k-2

def mixed(u,v,z):
 return count(u+pair(1)+v+pair(1)+z)-count(u+pair(1)+v+pair(0)+z)-count(u+pair(0)+v+pair(1)+z)+count(u+pair(0)+v+pair(0)+z)

def ends(n,L,h):
 a=(n-L)//2;b=n-L-a
 def root(power):return math.isqrt(h<<power) if power>=0 else math.isqrt(h>>(-power))
 A=min((1<<a)-1,root(2*a-n));B=min((1<<b)-1,root(2*b-n))
 return bits(A,a),tuple(1-x for x in bits(B,b)[::-1]),(1<<(n-a))+(1<<(n-b))+n//2+1

def select(n,k,r,h):
 u,v,E=ends(n,length(k),h);base=count(u+bank(k,0)+v)
 e=((r-base)%(1<<k))//2
 return u+bank(k,e)+v,E

# Independent slot count and vertex-subset DP.
for n in range(13):
 for A in range(1<<n):
  w=bits(A,n)
  require(count(w)==slot(w),'slot identity','slot')
  if n<=7:require(count(w)==dp(w),'graph path count','dp')

# Exact mixed derivative: touching pairs are the special hostile boundary.
short=[w for n in range(4) for w in itertools.product((0,1),repeat=n)]
for u,v,z in itertools.product(short,repeat=3):
 expected=2 if len(v)==0 else 1<<len(v)
 require(mixed(u,v,z)==expected,'mixed derivative','mixed')
require(mixed((),(),())==2,'0101 touching-pair hostile','touching_hostile')
for g in range(1,21):
 for v in ((0,)*g,(1,)*g,tuple(i%2 for i in range(g))):
  require(mixed((1,0),v,(1,0))==1<<g,'long mixed derivative','mixed_long')

banks=[]
for k in range(2,13):
 M=1<<k;L=length(k);D=spacers(k)
 require(L==len(bank(k,0)),'bank length','length')
 for j in range(1,k):
  left=integer(D[j-1]);right=sum((1-b)<<t for t,b in enumerate(D[j]))
  eps=1-D[j][0]
  delta=left+2-D[j-1][-1]+right-eps
  require(delta%M==1<<j,'spacer recursion','spacer')
 for u,v in [((),()),((1,),(0,)),((0,1,1),(1,0,0)),((1,0,1,1),(0,1,1,0))]:
  base=count(u+bank(k,0)+v)
  deltas=[count(u+bank(k,1<<j)+v)-base for j in range(k-1)]
  linear_carries=[(deltas[j]-(1<<(j+1)))//M for j in range(k-1)]
  require(all(x>=0 for x in linear_carries),'nonnegative linear carries','carry_sign')
  residues=[]
  for e in range(1<<(k-1)):
   value=count(u+bank(k,e)+v);residues.append(value%M)
   require((value-base)%M==2*e,'bank exact residue code','bank')
   digits=[(e>>j)&1 for j in range(k-1)]
   carry=sum(c*b for c,b in zip(linear_carries,digits))+sum((1<<((j-i-1)*(k+2)))*digits[i]*digits[j] for i in range(k-1) for j in range(i+1,k-1))
   require(value==base+2*e+M*carry,'exact carry ledger','carry')
  require(sorted(residues)==list(range(1,M,2)),'onto odd residues','onto')
 banks.append({'k':k,'old_length':2*k*k-2,'sharp_independent_switch_length':L,'spacers':[''.join(map(str,d)) for d in D]})

joint=[]
for k in range(1,17):
 L=length(k)
 for extra in (2,3,12,25):
  n=L+extra
  for h in (1<<(n-2),3*(1<<(n-2)),(1<<n)-1,1<<n):
   residues=range(1,1<<k,2) if k<=7 else (1,3,(1<<k)-3,(1<<k)-1)
   for r in residues:
    w,E=select(n,k,r,h)
    require(len(w)==n and w[0]==1 and w[-1]==0,'strong selector type','joint_type')
    require(count(w)%(1<<k)==r and abs(count(w)-h)<=E,'simultaneous selector','joint')
  joint.append({'k':k,'n':n,'L':L,'error':E})

# Generic separated two-ended certificate: information budget forbids E<2^k.
for k in range(2,101):
 for L in (k-1,k,k*k+2*k-2):
  for extra in (2,3,10,100):
   n=L+extra;a=(n-L)//2;b=n-L-a
   E=(1<<(n-a))+(1<<(n-b))+n//2+1
   require(E>=(1<<(L+2))>=(1<<(k+1)),'information/error barrier','generic_barrier')

# Nonlinear finite automaton bank: H(binary A) is unchanged by leading zeros.
# We enumerate ALL binary words, not only strong words. Minimal lengths are exact.
H=[1];ks=range(2,15);witness={k:{1:0} for k in ks};found={};table=[];previous={}
for n in range(1,30):
 previous={k:len(witness[k]) for k in ks if k not in found}
 for A in range(1<<(n-1),1<<n):
  p=A>>1;h=H[p] if A&1 else H[p]+p+(p&1);H.append(h)
  for k in ks:
   if k not in found:witness[k].setdefault(h% (1<<k),A)
 for k in ks:
  if k not in found and len(witness[k])==1<<(k-1):
   found[k]=n
   require(previous[k]<(1<<(k-1)),'minimal word length','automaton_minimal')
   table.append({'k':k,'minimal_free_word_length':n,'insulated_middle_length':2*k+n,'previous_length_residues':previous[k],'odd_residues':1<<(k-1)})
 if len(found)==len(list(ks)):break
require(found=={2:2,3:5,4:6,5:7,6:9,7:10,8:12,9:13,10:14,11:15,12:16,13:17,14:19},'finite table regression','automaton_table')
for k,n in found.items():
 M=1<<k;zeros=(0,)*k;ones=(1,)*k
 for u,v in [((),()),((1,),(0,)),((1,0,1),(0,1,0))]:
  base=count(u+zeros+ones+v)
  for r,A in witness[k].items():
   t=bits(A,n)
   require(count(t)%M==r,'nonlinear witness','automaton_witness')
   require((count(u+zeros+t+ones+v)-base)%M==(r-1)%M,'universal insulation interface','automaton_insulation')
 # Simultaneous finite selector using the same two-ended construction.
 L=2*k+n;total=L+12
 for h in (1<<(total-2),3*(1<<(total-2)),(1<<total)-1):
  u,v,E=ends(total,L,h);base=count(u+zeros+ones+v)
  for r in range(1,M,2):
   needed=(r-base+1)%M;t=bits(witness[k][needed],n);w=u+zeros+t+ones+v
   require(len(w)==total and w[0]==1 and w[-1]==0 and count(w)%M==r and abs(count(w)-h)<=E,'nonlinear joint selector','automaton_joint')

receipt={'status':'PASS; proved scoped architecture theorem plus finite nonlinear bank census','checks':CHECKS,'sharp_banks':banks,'simultaneous_tests':joint,'finite_nonlinear_banks':table,'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
(HERE/'receipt.json').write_text(json.dumps(receipt,indent=2,sort_keys=True)+'\n')
print(json.dumps({'status':'PASS','checks':CHECKS,'nonlinear_lengths':found,'receipt_sha256':hashlib.sha256((HERE/'receipt.json').read_bytes()).hexdigest()},sort_keys=True))
