"""Check the multiplicative scope and literal first-hole repair."""
from pathlib import Path
from math import isqrt
from hashlib import sha256
from json import loads
def ck(t,why):
 if not t:raise AssertionError(why)
def prime(n):return n>1 and all(n%d for d in range(2,isqrt(n)+1))
def hull(A,B,lim=59999):
 C=set(A)
 buff=sorted(B-{1})
 for k in range(1,lim+1,2):
  if k not in C:continue
  for b in buff:
   if k*b>lim:break
   C.add(k*b)
 return C
def starts(A):
 n=len(A);F=(1<<n)-1
 D=[[0]*n for _ in range(F+1)]; H=[0]*(F+1);H[0]=1
 for s in range(1,F+1):
  for i in range(n):
   if not s>>i&1:continue
   other=s^(1<<i)
   D[s][i]=1 if not other else sum(D[other][j] for j in range(n) if (other&A[i])>>j&1)
   H[s]+=D[s][i]
 return D,H
def strong(A):
 n=len(A);F=(1<<n)-1
 for R in (A,tuple(F^x^(1<<i) for i,x in enumerate(A))):
  s=1
  while True:
   a=s
   for i,r in enumerate(R):
    if s>>i&1:a|=r
   if a==s:break
   s=a
  if s!=F:return False
 return True
def run():
 f=Path(__file__).with_name('pinned_counts_9.json').read_bytes()
 ck(sha256(f).hexdigest()=='bd44780743699a2c65dc8b355c9b7f1435b8cfb2ae464da1a624f10b91ef5de1','count snapshot hash')
 data=loads(f)
 ck(data['census_sha256']=='1e0c721d111f51a6a1c5e81032b557a5123920d6daf33993fc15a123db50001d','historical bank identity')
 A=set(data['counts']);ck(len(A)==1520 and max(A)==3357 and 7 not in A,'bank')
 B={2**j+1 for j in range(1,16)}
 sets=(hull(A,{3,5,17,257}),hull(A,B),hull({1},A),hull({1},A|B))
 targets=set(range(69,60000,2))
 misses=[sorted(targets-C) for C in sets]
 ck([len(x) for x in misses]==[22013,21209,10849,10849],'miss totals')
 ck(all(x[0]==2887 for x in misses) and prime(2887),'first miss')
 p={k for k in targets if prime(k)}
 p7={k for k in targets if k%7==0 and prime(k//7)}
 ck(len(p-A)==5622 and len(p7-A)==983,'lane miss totals')
 ck(all(not ((p|p7)&C)-A for C in sets),'unchanged lanes')
 print('multiplicative misses:',[len(x) for x in misses],'; first 2887')
 w=data["H911_upper_bits"];root=8;R=[0]*9;j=0
 ck(w=='111110101111111111100111111101111111','literal source')
 for x in range(9):
  for y in range(x+1,9):
   if w[j]=='1': R[x]|=1<<y
   else:R[y]|=1<<x
   j+=1
 dp,H=starts(R);F=511^(1<<root)
 L=0
 for s in range(F+1):
  if s&~F:continue
  L+=H[F^s]*(sum(dp[s][i] for i in range(9) if (s&R[root])>>i&1) if s else 1)
 ck(H[511]==911 and L==988 and strong(R),'rooted source')
 for i in range(8):R[i]|=(1<<9)|(1<<10)
 R[8]|=1<<9
 R.extend((1<<10,1<<8))
 ck(starts(R)[1][-1]==2887 and strong(R),'literal glue')
 print('first multiplicative hole repaired by one triangle overlap:',H[511],L,2887)
if __name__=='__main__':run()
