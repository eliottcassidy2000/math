import json,pathlib,hashlib,random
from fractions import Fraction as Q
D=pathlib.Path(__file__).parent
checks={'physical_joint_terms':0,'physical_masks':0,'attachment_counts':0,'rational_identities':0,'escape_bounds':0,'clock_shifts':0}
def require(x,why):
 if not x:raise RuntimeError(why)
def joint(f):
 # Scan high endpoint; the accumulator counts old low endpoint choices.
 rev=0;j=0;prev=0
 for k in range(f.bit_length()):
  bit=(f>>k)&1
  if bit:j+=rev+prev
  rev=2*rev+bit;prev=bit
 return j
# Covers are literal disjoint path slots; transitive slot orders are forced.
for n in range(1,9):
 N=1<<n;J=[[0]*n for _ in range(n)];F=[0]*n
 for slot1 in range(1,N):
  u=slot1.bit_length()-1;F[u]+=1;slot2=(N-1)^slot1
  if slot2:J[u][(slot2&-slot2).bit_length()-1]+=1
 for i in range(n):
  for k in range(n):
   expected=(1<<(i-k-1))if i>k else (1 if k==i+1 else 0)
   require(J[i][k]==expected,'directed joint formula');checks['physical_joint_terms']+=1
 for f in range(1,N):
  ports=[i for i in range(n)if f>>i&1]
  require(sum(F[i]for i in ports)==f,'physical F mask')
  require(sum(J[i][k]for i in ports for k in ports)==joint(f),'physical joint mask');checks['physical_masks']+=1
# Full path DP independently verifies strong TT1/root/TT_m chart.
for n in range(1,7):
 for f in range(1<<(n-1),1<<n):
  v=n+2;rows=[0]*v
  rows[0]=sum(1<<i for i in range(2,v));rows[1]=1
  for i in range(n):
   rows[i+2]=sum(1<<(k+2)for k in range(i+1,n))
   if f>>i&1:rows[i+2]|=2
   else:rows[1]|=1<<(i+2)
  E=[[0]*v for _ in range(1<<v)]
  for i in range(v):E[1<<i][i]=1
  for S in range(1,1<<v):
   for i in range(v):
    if S>>i&1:
     for k in range(v):
      if not S>>k&1 and rows[i]>>k&1:E[S|1<<k][k]+=E[S][i]
  require(sum(E[-1])==1-joint(f)+2*f,'full graph count');checks['attachment_counts']+=1
  for start in range(v):
   reach=1<<start
   while True:
    nxt=reach
    for i in range(v):
     if reach>>i&1:nxt|=rows[i]
    if nxt==reach:break
    reach=nxt
   require(reach==(1<<v)-1,'full graph strongness')
rng=random.Random(183)
samples=list(range(1,257))+[rng.getrandbits(k)|1|(1<<(k-1))for k in range(1,257)]
for u in samples:
 bits=[i for i in range(u.bit_length())if u>>i&1];w=len(bits);a=sum(i+1 in bits for i in bits)
 R=sum((Q(1,1<<i)for i in bits),Q());B=sum((Q(1,1<<(k-i))for i in bits for k in bits if i<k),Q())
 require(Q(joint(u))-u*R/2==a-Q(w,2)-B/2,'exact rational identity')
 require(abs(Q(joint(u),u)-R/2)<=Q(u.bit_length(),u),'rational error bound');checks['rational_identities']+=1
 if u%2:
  R3=sum((Q(1,1<<i)for i in range((3*u).bit_length())if (3*u)>>i&1),Q())
  require(abs(R3/2-R/2)>=Q(1,8),'mod8 gap')
for u in range(1,200000,2):
 require(8*abs(joint(3*u)-3*joint(u))>=3*u-32*u.bit_length()-16,'escape bound');checks['escape_bounds']+=1
for f in range(1,128):
 for k in range(-3,4):
  fp=Q(3*f,1<<k)if k>=0 else Q(3*f*(1<<(-k)))
  if fp.denominator!=1:continue
  fp=int(fp);d=3*joint(f)-joint(3*f)-2
  require(joint(fp)==joint(3*f),'shift invariance')
  for t in [4,5]:require(1-joint(fp)+(1<<(t+k))*fp==3*(1-joint(f)+(1<<t)*f)+d,'clock relation')
  odd=lambda x:x//(x&-x)
  require(odd(fp)==3*odd(f),'odd part transition');checks['clock_shifts']+=1
hits={str(d):[u for u in range(1,128,2)if 3*joint(u)-joint(3*u)-2==d]for d in [2,4]}
require(hits=={'2':[],'4':[7]},'bounded alphabet classification')
require(Q(3*128,8)-4*8-2>6,'cutoff')
for m in range(8,1000):require(Q(3*(1<<(m-1)),8)-4>0,'dyadic monotonicity')
receipt={'status':'PASS','checks':checks,'small_digit_hits':hits,'scope':'Independent exact finite controls; all-order symbolic theorem reviewed separately.','source_sha256':hashlib.sha256((D.parent/'induction/PROOF.md').read_bytes()).hexdigest()}
(D/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt,sort_keys=True))
