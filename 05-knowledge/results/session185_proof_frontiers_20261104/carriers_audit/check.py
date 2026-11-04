import itertools,json,hashlib
from pathlib import Path
P=Path('/workspace/shared/session185/carriers')
def need(b,m):
 if not b:raise RuntimeError(m)
def table(a):
 n=len(a);N=1<<n;dp=[[0]*n for _ in range(N)];h=[0]*N;h[0]=1
 for i in range(n):dp[1<<i][i]=1
 for s in range(1,N):
  h[s]=sum(dp[s])
  for i,c in enumerate(dp[s]):
   if c:
    avail=a[i]&~s
    while avail:
     b=avail&-avail;avail-=b;dp[s|b][b.bit_length()-1]+=c
 return h,dp
def sig(a,m):
 h,e=table(a);rev=[sum(1<<j for j in range(len(a)) if a[j]>>i&1) for i in range(len(a))];_,st=table(rev);N=(1<<len(a))-1
 E=[sum(c for i,c in enumerate(row) if m>>i&1) for row in e];S=[sum(c for i,c in enumerate(row) if m>>i&1) for row in st]
 return [h[N],sum(h[t]*h[N^t] for t in range(N+1)),sum(E[t]*h[N^t] for t in range(N+1)),sum(h[t]*S[N^t] for t in range(N+1)),sum(E[t]*S[N^t] for t in range(N+1)),E[N],S[N]]
def strong(a):
 n=len(a)
 for root in range(n):
  seen={root};todo=[root]
  while todo:
   u=todo.pop()
   for v in range(n):
    if a[u]>>v&1 and v not in seen:seen.add(v);todo.append(v)
  if len(seen)<n:return False
 return True
def canon(a):
 n=len(a);return min(tuple(sum(1<<j for j in range(n) if a[p[i]]>>p[j]&1) for i in range(n)) for p in itertools.permutations(range(n)))
def append(a,m,u):
 n=len(a);k=u.bit_length();full=((1<<k)-1)<<n
 return [v|full for v in a]+[sum(1<<j for j in range(n+i+1,n+k)) for i in range(k)],m|(u<<n)
def attach(a,m,t):
 n=len(a);b=[sum(1<<j for j in range(i+1,n+t)) for i in range(t)]+[v<<t for v in a]+[(1<<t)-1]
 for i in range(n):
  if m>>i&1:b[t+i]|=1<<(n+t)
  else:b[-1]|=1<<(t+i)
 return b
def odd(f):return f//(f&-f)
D=json.loads((P/'bank.json').read_text());reps=set();labelled={}
for n in range(1,6):
 count=0;pairs=list(itertools.combinations(range(n),2))
 for mask in range(1<<len(pairs)):
  a=[0]*n
  for b,(i,j) in enumerate(pairs):
   if mask>>b&1:a[i]|=1<<j
   else:a[j]|=1<<i
  if strong(a):reps.add(canon(a));count+=1
 labelled[n]=count
need(reps=={tuple(r) for r in D['rows']},'complete canonical bank')
need(len(D['bank'])==sum(1<<len(r) for r in D['rows']),'bank size')
seen=set();checked=0;hits=[]
for i,m,z in D['bank']:
 need((i,m) not in seen,'bank duplicate');seen.add((i,m));a=D['rows'][i];need(sig(a,m)==z,'bank signature')
 if not m:continue
 h,p,F,G,J,e,s=z;need(G>0,'G positivity')
 for C in (-35,-29):
  limit=(h-J-C)//G
  for u in range(max(0,limit)+1):
   aa,mm=append(a,m,u);ss=sig(aa,mm);checked+=1
   need(strong(attach(aa,mm,1)),'physical candidate strong')
   if ss[0]-ss[4]==C and odd(ss[2])==63:hits.append((i,m,u,C))
need(checked==867 and hits==[],'inverse receipt')
a=[2,12,17,21,3];m=8;need(strong(a),'example strong');need(sig(a,m)==[11,100,12,27,0,1,4],'example signature')
aa,mm=append(a,m,645);ss=sig(aa,mm);need(ss[2]==64512 and ss[0]-ss[4]==-21849,'first family')
full=attach(aa,mm,1);need(strong(full) and table(full)[0][-1]==107175,'physical first attachment')
need(pow(2,20,25)==1 and 63*2**8%25==3,'all-order integrality')
wits=json.loads((P/'escape_witnesses.json').read_text());ww=[]
for w in wits:
 a=w['rows'];m=w['mask'];ss=sig(a,m);need(ss==w['sig'] and strong(a),'witness sig and SCC')
 counts=[]
 for t in (1,2,3):
  b=attach(a,m,t);H=table(b)[0][-1];need(strong(b) and H==ss[0]-ss[4]+2**t*ss[2],'witness full graph');counts.append(H)
 ww.append([w['address'],counts])
print(json.dumps({'status':'PASS','labelled_strong_counts':labelled,'strong_types':len(reps),'bank_rows':len(D['bank']),'physically_recomputed_inverse_candidates':checked,'hits':hits,'first_attachment_H':107175,'witnesses':ww,'source_sha256':{n:hashlib.sha256((P/n).read_bytes()).hexdigest() for n in ['PROOF.md','verify.py','bank.json','inherited.py','escape_witnesses.json']}},indent=2,sort_keys=True))
