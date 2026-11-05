import json,math,sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).parent/'source'))
from converse import diagonal_inverse
from bilinear import source_bank

def need(b,s):
 if not b:raise RuntimeError(s)
# Recompute the gadget constants from literal subset Hamilton words, not g/Q recurrences.
def gadget_constants(L):
 bits=[]if L==1 else[1-int(c)for c in bin(L-1)[2:][::-1]];n=len(bits)+1;a=[0]*n
 for i,b in enumerate(bits,1):a[0 if b else i]|=1<<(i if b else 0)
 for i in range(1,n):
  for j in range(i+1,n):a[i]|=1<<j
 def count(S,root=None):
  total=0
  def dfs(v,left):
   nonlocal total
   if not left:total+=1;return
   z=a[v]&left
   while z:
    b=z&-z;z-=b;dfs(b.bit_length()-1,left^b)
  if not S:return 1 if root is None else 0
  for v in range(n):
   if S>>v&1 and(root is None or v==root):dfs(v,S^(1<<v))
  return total
 allset=(1<<n)-1;h=[count(s)for s in range(1<<n)];R=sum(count(s,0)*h[allset^s]for s in range(1<<n));return h[-1],R,1<<(n-1)
cap=57;constants={L:gadget_constants(L)for L in range(1,cap+1)}
# Layer by the strictly increasing U coordinate rather than the source BFS implementation.
levels=[set()for _ in range(cap+1)];levels[1].add((1,1,1,1));bank=set()
for U in range(1,cap+1):
 for h,e,u,d in levels[U]:
  bank.add((h,e,u,d));moves=[(h,0,2*u,d)]
  for L in range(1,cap//u+1):
   g,Q,B=constants[L];moves.append((g*h-(L%2)*e+L*u,h-(e if L==1 else 0),Q*h+B*u,h))
  for s in moves:
   need(s[2]>U,'strict source height')
   if s[2]<=cap:levels[s[2]].add(s)
need(bank==set(source_bank(cap)),'independent physical-constant state bank')
def val(v,kind):
 D,E,U,d=v
 return (D-E)**2+2*D*U+3*U*U-2*d*(E+U)if kind=='N'else D*D+4*D*U+6*U*U-2*E*U-2*d*(E+U)
image={15}
for v in bank-{(1,1,1,1)}:
 D,E,U,d=v
 for kind in('N','M'):
  h=val(v,kind);low=3*U*U+1 if kind=='N'else 6*U*U+2*U+1
  need(h>=low and h%4==1,'diagonal lower bound and parity')
  if h<=10000:image.add(h)
need(len(image)==194,'complete support through10000')
for target in range(1,1001):need(bool(diagonal_inverse(target))==(target in image),'target inverse against independent bank')
# Formal GF identities to order100, with literal Fibonacci recurrence.
f=[0,1]
for _ in range(200):f.append(sum(f[-2:]))
w=[0]*102;w[1]=1
for n in range(2,102):w[n]=f[2*n-3]
a=[0]*102;a[5]=1
for k in range(2,50):a[2*k+1]+=w[k];a[2*k+2]+=w[k]
# Multiply SC2 by 1-3z²+z⁴; expected z5(2+z-4z²-z³+z⁴).
for n in range(101):
 lhs=a[n]-(3*a[n-2]if n>=2 else 0)+(a[n-4]if n>=4 else 0)
 rhs={5:2,6:1,7:-4,8:-1,9:1}.get(n,0)
 need(lhs==rhs,'exact rational GF')
print(json.dumps({'status':'PASS','literal_gadget_constants':cap,'independent_source_states':len(bank),'complete_diagonal_support_through10000':len(image),'exact_inverse_targets':1000,'GF_coefficients':101},sort_keys=True,indent=2))
