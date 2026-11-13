from fractions import Fraction
from math import lcm
from itertools import combinations
def first(P,n):
 for f in P:
  if n%f.denominator==0:return n//f.denominator*f.numerator
 return n
def compile_program(code,registers,states,clock):
 # instructions ('inc',reg,dst), ('dec',reg,positive,zero), ('halt',)
 p={}
 for j,I in enumerate(code):
  for z in (0,1):
   src=states[j][z]
   if I[0]=='inc':
    _,r,d=I;p[src]=[Fraction(clock*states[d][1-z]*registers[r],src)]
   elif I[0]=='dec':
    _,r,d,e=I;p[src]=[Fraction(clock*states[d][1-z],src*registers[r]),Fraction(clock*states[e][1-z],src)]
   else:p[src]=[Fraction(1,src)]
 return tuple(f for L in p.values() for f in L)
# Nonterminal transition replay, including self-loop through the other layer.
C=[('inc',0,1),('dec',0,1,2),('halt',)]
R=(2,); P=((5,7),(11,13),(17,19)); clock=3
f=compile_program(C,R,P,clock)
mod=lcm(*(g.denominator for g in f))
for j,I in enumerate(C):
 for z in (0,1):
  for a in range(6):
   for t in range(5):
    n=P[j][z]*2**a*3**t
    if j==0:w=P[1][1-z]*2**(a+1)*3**(t+1)
    elif j==1 and a:w=P[1][1-z]*2**(a-1)*3**(t+1)
    elif j==1:w=P[2][1-z]*3**(t+1)
    else:w=2**a*3**t
    if first(f,n)!=w:raise ValueError((j,z,a,t))
    slope=next((q for q in f if (n%mod)%q.denominator==0),1)
    if slope*n!=w:raise ValueError('residue rule')
# A separate cleanup: token 23 switches to 29 while draining 2,
# then two different tokens 31,37 drain the clock and remove state.
cleanup=(Fraction(29,46),Fraction(31,23),Fraction(23,58),Fraction(31,29),
         Fraction(37,93),Fraction(1,31),Fraction(31,111),Fraction(1,37))
for a in range(7):
 for t in range(7):
  n=23*2**a*3**t
  for _ in range(a+t+4):
   n=first(cleanup,n)
  if n!=1:raise ValueError(('cleanup',a,t,n))
print('prime-state compiler, fixed residue branches, and cleanup pass')
def count(A):
 n=len(A);F=[[0]*n for _ in range(1<<n)]
 for i in range(n):F[1<<i][i]=1
 for S in range(1,1<<n):
  for i,c in enumerate(F[S]):
   if not c:continue
   r=A[i]&~S
   while r:
    b=r&-r;r^=b;F[S|b][b.bit_length()-1]+=c
 return sum(F[-1])
def through(maxn):
 H=set()
 for n in range(1,maxn+1):
  pairs=list(combinations(range(n),2))
  for mask in range(1<<len(pairs)):
   A=[0]*n
   for i,(u,v) in enumerate(pairs):
    if mask>>i&1:A[u]|=1<<v
    else:A[v]|=1<<u
   H.add(count(A))
 return H
if through(4)!={1,3,5}:raise ValueError('small spectrum')
if {x for x in through(5) if x<=11}!={1,3,5,9,11}:raise ValueError('small cutoff')
print('literal finite spectrum checker controls pass')
def odd_reset(n):
 if not n%2:n+=1
 return (n+1)//2 if n%4==1 else (n-1)//2
for n in range(1,2000):
 q=n if n%2 else n+1
 while q!=1:
  j=odd_reset(q)
  if not (0<j<q and j%2):raise ValueError('reset lane')
  q=j
print('positive-slope odd reset lane checked')
