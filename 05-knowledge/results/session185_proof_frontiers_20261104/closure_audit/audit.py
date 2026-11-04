"""Independent permutation path counts for two literal carriers; standard library."""
import itertools,json
from pathlib import Path

def need(ok,msg):
 if not ok:raise RuntimeError(msg)
def sig(rows,mask):
 n=len(rows);N=(1<<n)-1
 h=[0]*(N+1);e=[0]*(N+1);s=[0]*(N+1);h[0]=1
 for k in range(1,n+1):
  for path in itertools.permutations(range(n),k):
   if all(rows[a]>>b&1 for a,b in zip(path,path[1:])):
    S=sum(1<<v for v in path);h[S]+=1
    e[S]+=(mask>>path[-1])&1;s[S]+=(mask>>path[0])&1
 return (h[N],sum(h[A]*h[N^A] for A in range(N+1)),sum(e[A]*h[N^A] for A in range(N+1)),sum(h[A]*s[N^A] for A in range(N+1)),sum(e[A]*s[N^A] for A in range(N+1)),e[N],s[N])
def mul(x,y):
 h,p,F,G,J,e,s=x;H,P,f,g,j,E,S=y
 return(h*H,p*P,F*H+p*f,G*P+h*g,J*H+h*j+e*S+G*f,h*E,s*H)

def strong(r):
 n=len(r)
 reach=[[i==j or bool(r[i]>>j&1) for j in range(n)] for i in range(n)]
 for k in range(n):
  for i in range(n):
   for j in range(n):reach[i][j]=reach[i][j] or (reach[i][k] and reach[k][j])
 return all(map(all,reach))
def attach(r,mask,t):
 n=len(r);out=[0]*(n+t+1);root=n+t
 for i in range(n+t+1):
  for j in range(i+1,n+t+1):
   if j==root:forward=(i>=t and bool(mask>>(i-t)&1))
   elif i<t:forward=True
   else:forward=bool(r[i-t]>>(j-t)&1)
   if forward:out[i]|=1<<j
   else:out[j]|=1<<i
 return out
def count(r):
 n=len(r);full=(1<<n)-1;total=0
 def walk(mask,last):
  nonlocal total
  if mask==full:total+=1;return
  remaining=r[last]&~mask
  while remaining:
   b=remaining&-remaining;remaining-=b;walk(mask|b,b.bit_length()-1)
 for v in range(n):walk(1<<v,v)
 return total
ws=json.loads((Path(__file__).parent.parent/'closure'/'witnesses.json').read_text())
states={'':([6,4,0],7),'2':([30,28,8,16,4],7),'4':([30,28,24,16,0],21)}
for x in ws:states[x['address']]=(x['rows'],x['mask'])
sigs={a:sig(r,m) for a,(r,m) in states.items()}
for x in ws:need(sigs[x['address']]==tuple(x['sig']),'witness signature')
shifts={'':0};edges=[];fullcounts={}
from fractions import Fraction
for a,(r,m) in states.items():
 n=len(r)
 need(all(0<=x<1<<n and not(x>>i&1) for i,x in enumerate(r)),'row range')
 need(all(((r[i]>>j)&1)+((r[j]>>i)&1)==1 for i in range(n) for j in range(i)),'tournament')
 if len(a)>=2:need(strong(r),'strong suffix')
 z=sigs[a];c=z[0]-z[4]
 if a:
  p=sigs[a[:-1]];d=int(a[-1]);qs=[q for q in range(-20,21) if Fraction(2)**q*z[2]==3*p[2]]
  need(len(qs)==1,'unique clock shift');q=qs[0]
  need(c==3*(p[0]-p[4])+d,'constant');shifts[a]=shifts[a[:-1]]+q
  edges.append((a[:-1],a,q))
 for t in (1,2,4+shifts[a]):
  need(t>=1,'legal common clock');A=attach(r,m,t);need(strong(A),'strong full attachment')
  H=count(A);need(H==c+2**t*z[2],'independent full Hamilton count')
  if t==4+shifts[a]:fullcounts[a]=H
for a in states:
 if a:need(fullcounts[a]==3*fullcounts[a[:-1]]+int(a[-1]),'common-clock tree')
# Stronger author's final valuation family, independently evaluate the transitive block.
r=(24,1,3,6,14);seed=sig(r,14)
need(strong(r) and seed==(11,100,64,49,28,9,6),'final valuation seed')
A=671088;n=20;ports=[i for i in range(n) if A>>i&1]
# In a transitive two-slot partition with max(first)=i,min(second)=j:
# i>j gives 2^(i-j-1) free middle choices; i<j only adjacent endpoints work.
J=sum(2**(i-j-1) if i>j else int(j==i+1) for i in ports for j in ports if i!=j)
V=(1,2**n,sum(2**i for i in ports),sum(2**(n-1-i) for i in ports),J,int(n-1 in ports),int(0 in ports))
need(V==(1,1048576,671088,60357,38626,1,0),'independent boundary partition block')
need(25*A==16*(2**20-1),'geometric identity')
z=seed
for k in range(1,101):
 z=mul(z,V)
 need(z[2]==2**(6+20*k) and (z[0]-z[4])%3==(1+k)%3,'valuation residue formula')
need(tuple(x%3 for x in V)==(1,1,0,0,1,1,0),'exact residue recurrence')
result={'status':'PASS','permutation_signature_states':len(states),'directed_DFS_full_attachment_counts':fullcounts,'edges':edges,'valuation_seed_signature':seed,'independent_transitive_block_signature':V,'valuation_instances':100,'scope':'Finite verification plus written audit of all-order claims.'}
Path(__file__).with_name('receipt.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,sort_keys=True))
