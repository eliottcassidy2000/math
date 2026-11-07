from collections import defaultdict,Counter
from fractions import Fraction as Q
from itertools import product
from math import comb

def need(ok,msg):
 if not ok: raise RuntimeError(msg)

def move(s,c,m=3):
 D,E,U,d=s
 return tuple(x%m for x in ((D,0,2*U,d) if c=='O' else (D-E+U,d,2*U-d,d) if c=='Z' else (D-E+U,D-E,D+U,D)))
def F(s):return (s[0]-s[1]+s[2])%3
def unit(a):return { (D,e*a,U,a):Q(1,18) for D,U,e in product(range(3),range(3),range(2)) }
def add(out,a,c):
 for s,p in a.items():out[s]+=c*p

def project(s):
 D,E,U,d=s
 if d:return unit(d)
 if U:return {(x,0,u,0):Q(1,6) for x in range(3) for u in (1,2)}
 z=defaultdict(Q);z[(D,0,0,0)]+=Q(1,2);z[((D-E)%3,0,0,0)]+=Q(1,2)
 return z

def Jproj(a):
 out=defaultdict(Q)
 for s,p in a.items():add(out,project(move(s,'J')),p)
 return out

def residue(a):
 out=[Q(0)]*3
 for s,p in a.items():out[F(s)]+=p
 return out

# Independent finite graph control for the claimed O/Z recurrent laws.
for a in (1,2):
 S=unit(a); out=defaultdict(Q)
 for s,p in S.items():
  for c in 'OZ':out[move(s,c)]+=p/2
 need(dict(out)==S,"stationary measure")
 for root in S:
  seen={root};todo=[root]
  for s in todo:
   for c in 'OZ':
    v=move(s,c)
    if v not in seen:seen.add(v);todo.append(v)
  need(set(S)==seen,"strong connectivity")
 need(any(move(s,'O')==s for s in S),'aperiodicity')
pi=unit(1)
for k in range(3):
 print('PROVED projector k=',k,'probabilities=',residue(pi),'local factors=',[3*x for x in residue(pi)])
 pi=Jproj(pi)
need(residue(Jproj(Jproj(unit(1))))==[Q(34,81),Q(47,162),Q(47,162)],"two-J law")
# Entire-length exact integer recurrence, no stationary formula.
a=[{(1,0,2,1):1},{},{}]
for L in range(1,1001):
 b=[defaultdict(int) for _ in range(3)]
 for j in range(3):
  for s,v in a[j].items():
   for c in 'OZ':b[j][move(s,c)]+=v
   if j<2:b[j+1][move(s,'J')]+=v
 a=b
 if L<=9:
  lit=Counter()
  for w in product('OZJ',repeat=L):
   if w.count('J')!=2:continue
   s=(1,0,2,1)
   for c in w:s=move(s,c,10**12)
   lit[F(s)]+=1
  need([lit[r] for r in range(3)]==[sum(v for s,v in a[2].items() if F(s)==r) for r in range(3)],"literal enumeration")
 if L in (20,80,200,1000):
  den=comb(L,2)*2**(L-2)
  print('FINITE-EXACT length',L,'factors',[float(Q(3*sum(v for s,v in a[2].items() if F(s)==r),den)) for r in range(3)])
print('PASS: stationary invariance, strong connectivity, self-loops, literal-word control through length 9, exact recurrence through length 1000')
