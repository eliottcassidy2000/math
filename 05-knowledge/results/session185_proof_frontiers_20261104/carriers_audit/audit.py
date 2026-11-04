import json,itertools,hashlib
from pathlib import Path
P=Path('/workspace/shared/session185/carriers')
def need(x,m):
 if not x:raise RuntimeError(m)
def strong(r):
 n=len(r);a=[[i==j or bool(r[i]>>j&1) for j in range(n)] for i in range(n)]
 for k in range(n):
  for i in range(n):
   for j in range(n):a[i][j]|=a[i][k] and a[k][j]
 return all(map(all,a))
def signature(r,W):
 n=len(r);N=(1<<n)-1;h=[0]*(N+1);h[0]=1;e=[0]*(N+1);s=[0]*(N+1)
 def rec(mask,first,last):
  h[mask]+=1;e[mask]+=last in W;s[mask]+=first in W
  for v in range(n):
   if not(mask>>v&1) and r[last]>>v&1:rec(mask|1<<v,first,v)
 for v in range(n):rec(1<<v,v,v)
 return h[N],sum(h[x]*h[N^x] for x in range(N+1)),sum(e[x]*h[N^x] for x in range(N+1)),sum(h[x]*s[N^x] for x in range(N+1)),sum(e[x]*s[N^x] for x in range(N+1)),e[N],s[N]
def tailjoint(u):
 n=u.bit_length();N=(1<<n)-1
 return sum(bool(u>>(S.bit_length()-1)&1) and bool(u>>(((N^S)&-(N^S)).bit_length()-1)&1) for S in range(1,N))
def forward(r):
 n=len(r);dp=[{} for _ in range(1<<n)]
 for v in range(n):dp[1<<v][v]=1
 for mask in range(1,1<<n):
  for end,c in dp[mask].items():
   avail=r[end]&~mask
   while avail:
    b=avail&-avail;avail-=b;v=b.bit_length()-1;d=dp[mask|b];d[v]=d.get(v,0)+c
 return sum(dp[-1].values())
data=json.loads((P/'bank.json').read_text());reps=data['rows'];bank=data['bank'];counts={}
need(len(reps)==9,'nine strong types')
for n in (1,3,4,5):
 rr=[r for r in reps if len(r)==n];pairs=list(itertools.combinations(range(n),2));count=0
 for bits in range(1<<len(pairs)):
  r=[0]*n
  for k,(u,v) in enumerate(pairs):
   if bits>>k&1:r[u]|=1<<v
   else:r[v]|=1<<u
  if not strong(r):continue
  count+=1
  need(any(any(all(bool(r[i]>>j&1)==bool(q[p[i]]>>p[j]&1) for i,j in pairs) for p in itertools.permutations(range(n))) for q in rr),'bank coverage')
 counts[n]=count
need(counts=={1:1,3:2,4:24,5:544},'strong census')
need({(i,m) for i,m,z in bank}=={(i,m) for i,r in enumerate(reps) for m in range(1<<len(r))},'all masks')
independent=[]
for i,m,z in bank:
 r=reps[i];need(strong(r),'bank strong');zz=signature(r,{k for k in range(len(r)) if m>>k&1});need(tuple(z)==zz,'bank signature');independent.append((i,m,zz))
checks=0;hits=[]
for i,m,z in independent:
 if not m:continue
 h,p,F,G,J,e,s=z
 for C in (-35,-29):
  need(G>0,'positive G');bound=(h-J-C)//G
  for u in [0]+list(range(1,bound+1)):
   checks+=1
   if u:f=F+p*u;c=h-J-h*tailjoint(u)-e*(u&1)-G*u
   else:f=F;c=h-J
   if c==C and f//(f&-f)==63:hits.append((i,m,u))
need(not hits,'no targets')
need(checks==json.loads((P/'receipt.json').read_text())['fixed_front_tail_candidates_checked'],'candidate count')
r=(2,12,17,21,3);need(signature(r,{3})==(11,100,12,27,0,1,4),'seed')
for k in range(12):
 u=(63*2**(20*k+8)-3)//25;need(25*u==63*2**(20*k+8)-3,'integer');need(12+100*u==63*2**(20*k+10),'ray')
r=list(r);W={3};u=645
for k in range(u.bit_length()):
 n=len(r);r=[x|(1<<n) for x in r]+[0]
 if u>>k&1:W.add(n)
n=len(r);rows=[sum(1<<j for j in range(1,n+1))]+[(x<<1)|((1<<(n+1)) if i in W else 0) for i,x in enumerate(r)]+[1|sum(1<<(i+1) for i in range(n) if i not in W)]
need(strong(rows),'family full strong');j=tailjoint(u);c=11-11*j-(u&1)-27*u;H=forward(rows);need(H==c+2*(12+100*u),'full literal Hamilton count')
out={'status':'PASS','labelled_strong_census':counts,'actual_port_signatures':len(bank),'fixed_front_candidates':checks,'hits':hits,'first_family':{'u':u,'J_tail':j,'c':c,'F':12+100*u,'H':H,'order':len(rows)},'sha256':{f:hashlib.sha256((P/f).read_bytes()).hexdigest() for f in ['PROOF.md','verify.py','inherited.py','bank.json','escape_witnesses.json','receipt.json']}}
Path(__file__).with_name('receipt.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,sort_keys=True))
