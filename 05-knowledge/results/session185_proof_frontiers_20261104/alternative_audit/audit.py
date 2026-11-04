import itertools,math,json
from pathlib import Path

def need(ok,msg):
 if not ok:raise RuntimeError(msg)
Q=[(16,5,33,39,46,3),(16,5,49,39,42,3),(18,4,49,39,42,3)]
expected=[[(4,4,1),(5,4,1),(3,2,0),(12,10,2),(3,1,0)],[(4,4,1),(7,5,1),(4,2,0),(11,9,2),(3,1,0)],[(6,5,1),(5,4,1),(3,1,0),(12,9,2),(5,2,0)]]
def count(R):
 N=(1<<len(R))-1;ans=0
 def walk(S,v):
  nonlocal ans
  if S==N:ans+=1;return
  M=R[v]&~S
  while M:
   b=M&-M;M-=b;walk(S|b,b.bit_length()-1)
 for i in range(len(R)):walk(1<<i,i)
 return ans
def strong(R):
 n=len(R);d=[[i==j or bool(R[i]>>j&1) for j in range(n)] for i in range(n)]
 for k in range(n):
  for i in range(n):
   for j in range(n):d[i][j]|=d[i][k] and d[k][j]
 return all(map(all,d))
def sub(q,A):
 n=len(A);R=[x&31 for x in q[:5]]+[x<<5 for x in A]
 for i in range(5):
  for j in range(n):
   if q[i]>>5&1:R[i]|=1<<(5+j)
   else:R[5+j]|=1<<i
 return R
def covers(A,k):
 n=len(A);s=0
 for p in itertools.permutations(range(n)):
  for cuts in itertools.combinations(range(1,n),k-1):
   if all(j in cuts or A[p[j-1]]>>p[j]&1 for j in range(1,n)):s+=1
 return s if k<=n else 0
profiles=[]
for z,q in enumerate(Q):
 need(strong(q),'seed strong')
 cyc=(0,4,3,1,2,5,0);need(all(q[a]>>b&1 for a,b in zip(cyc,cyc[1:])),'common cycle')
 table=[[0]*3 for _ in range(5)]
 for p in itertools.permutations(range(5)):
  G={i for i in range(6) if (i==0 or q[p[i-1]]>>5&1) and (i==5 or q[5]>>p[i]&1)}
  B={i for i in range(1,5) if not(q[p[i-1]]>>p[i]&1)}
  if not B<=G:continue
  for k in range(1,4):
   if 0<=k-len(B)<=len(G)-len(B):table[p[0]][k-1]+=math.comb(len(G)-len(B),k-len(B))
 need(table==[list(t) for t in expected[z]],'full grouped profile')
 profiles.append([sum(row[k] for row in table) for k in range(3)])
for a,b,pair in [(Q[0],Q[1],(2,4)),(Q[1],Q[2],(0,1))]:
 diff=[(i,j) for i in range(6) for j in range(i+1,6) if ((a[i]^b[i])>>j)&1]
 need(diff==[pair],'single literal edge')
records=[]
for r in range(1,8):
 A=[sum(1<<j for j in range(i+1,r)) for i in range(r)]
 vals=[]
 for j,q in enumerate(Q):
  T=sub(q,A);need(strong(T),'inflated strong');h=count(T)
  need(h==4*3**r+9*2**r-3+2*j,'all-r finite direct count');vals.append(h)
 records.append(vals)
controls=0
for n in range(1,5):
 pairs=list(itertools.combinations(range(n),2))
 for m in range(1<<len(pairs)):
  A=[0]*n
  for bit,(u,v) in enumerate(pairs):
   if m>>bit&1:A[u]|=1<<v
   else:A[v]|=1<<u
  h=count(A);q2=covers(A,2);q3=covers(A,3)
  for j,q in enumerate(Q):
   T=sub(q,A);need(strong(T),'arbitrary module strong');need(count(T)==(27+2*j)*h+21*q2+4*q3,'arbitrary module count');controls+=1
need(pow(2,42,49)==pow(3,42,49)==1,'period42')
for r,res,p in [(4,28,67),(6,14,499)]:
 H=4*3**r+9*2**r+1;need(H%49==res and H==7*p,'exact residue/cofactor')
 need(all(p%d for d in range(2,math.isqrt(p)+1)),'cofactor prime')
for r in range(1,100):
 for j,l in itertools.product(range(3),repeat=2):
  need(-9*2**r+6+2*l-6*j<=-8,'uniform step inequality')
result={'status':'PASS','profiles':profiles,'grouped_profiles_verified':True,'direct_DFS_inflations':records,'permutation_cut_arbitrary_module_checks':controls,'mod49_period':42,'prime_cofactors':[67,499]}
Path(__file__).with_name('receipt.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,sort_keys=True))
