import itertools, random, math, json
rng=random.Random(172)
def req(ok,msg):
 if not ok:raise RuntimeError(msg)
def randgraph(n):
 a=[0]*n
 for i in range(n):
  for j in range(i):
   if rng.randrange(2):a[i]|=1<<j
   else:a[j]|=1<<i
 return a
def wrap(parts):
 groups=[(g,i) for g,p in enumerate(parts) for i in range(len(p))]; a=[0]*(len(groups)+1); last=len(groups)
 for i,(g,v) in enumerate(groups):
  for j,(h,w) in enumerate(groups):
   if g<h or (g==h and parts[g][v]>>w&1):a[i]|=1<<j
  if g==0:a[last]|=1<<i
  else:a[i]|=1<<last
 return a
def covers(a):
 n=len(a);p=[0]*(n+1)
 for seq in itertools.permutations(range(n)):
  for cut in range(1<<(n-1)):
   if all(cut>>j&1 or a[seq[j]]>>seq[j+1]&1 for j in range(n-1)):p[cut.bit_count()+1]+=1
 return [c//math.factorial(j) for j,c in enumerate(p)]
def val(p,u):return sum(c*math.prod(u-t for t in range(j)) for j,c in enumerate(p))
def form(ps,u):
 f=[val(p,u) for p in ps]; l=[val(p,u-1) for p in ps]; g=[val(p,u+1) for p in ps]
 return u*((l[0]-f[0])*math.prod(l[1:])+(3*f[0]-g[0]-l[0])*math.prod(f[1:])+(g[0]-f[0])*math.prod(g[1:]))
count=0
for sizes in [(1,1),(1,2),(2,1),(2,2),(1,3),(3,1),(1,1,1),(2,1,1),(1,2,1),(1,1,2),(2,2,1),(1,1,1,1),(3,2),(2,3)]:
 for repeat in range(3):
  parts=[randgraph(n) for n in sizes];p=covers(wrap(parts));ps=list(map(covers,parts))
  for u in range(-2,len(p)+1):req(val(p,u)==form(ps,u),'arbitrary block full profile')
  count+=1
B=[(16,1,3,7,14),(20,1,2,7,14)]
def query(p):return tuple(v|(32 if i else 0) for i,v in enumerate(p))+(1,)
def H(a):
 n=len(a);dp=[[0]*n for _ in range(1<<n)]
 for i in range(n):dp[1<<i][i]=1
 for s in range(1,1<<n):
  for j in range(n):
   ways=dp[s][j]
   if not ways:continue
   nxt=a[j]&~s
   while nxt:
    bit=nxt&-nxt;nxt-=bit;t=bit.bit_length()-1;dp[s|bit][t]+=ways
 return sum(dp[-1])
records=[]
for bits,address in [((0,0,0),None),((1,0,1),None),((0,1,0),0),((1,0,1),0),((0,0,1),2),((1,1,0),1)]:
 parts=[B[b] for b in bits]
 if address is not None:parts[address]=query(parts[address])
 actual=H(wrap(parts)); M=3
 if address is None:expected=75*84**2-57*9**2
 elif address==0:expected=(207+30*bits[0])*84**2-(173+26*bits[0])*9**2
 else:expected=(16800+2400*bits[address])*84-(969+114*bits[address])*9
 req(actual==expected,'three block readout')
 records.append({'bits':bits,'address':address,'H':actual})
print(json.dumps({'status':'PASS','arbitrary_profile_examples':count,'three_block_examples':records},indent=2))
