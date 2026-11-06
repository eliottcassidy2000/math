import itertools,json,pathlib
P=pathlib.Path(__file__).resolve().parent

def compositions(t):
 if t==0:yield ();return
 for w in range(2,t+1):
  for rest in compositions(t-w):yield(w,)+rest
physical=0
for t in range(8):
 for parts in compositions(t):
  n=t+1;a=[[i<j for j in range(n)] for i in range(n)];s=0
  for w in parts:a[s][s+w]=False;a[s+w][s]=True;s+=w
  h=sum(all(a[u][v] for u,v in zip(p,p[1:])) for p in itertools.permutations(range(n)))
  target=1
  for w in parts:target*=2**(w-1)+1
  assert h==target
  physical+=1
N=35;vals=[set() for _ in range(N+1)];vals[0]={1};cnt=[0]*(N+1);cnt[0]=1
for w in range(2,N+1):
 if w==4:continue
 q=2**(w-1)+1
 for t in range(w,N+1):vals[t]|={q*v for v in vals[t-w]};cnt[t]+=cnt[t-w]
p=[0]*(N+1);p[0]=1
for w in range(1,N+1):
 for t in range(w,N+1):p[t]+=p[t-w]
def pp(t):return p[t] if t>=0 else 0
seen=set()
for t in range(N+1):
 assert len(vals[t])==cnt[t]==pp(t)-pp(t-1)-pp(t-4)+pp(t-5)
 assert not seen.intersection(vals[t]);seen.update(vals[t])
 assert len(seen)==pp(t)-pp(t-4)
# Independent priority closure of numerical monoid, carrying its weighted order.
import heapq
cap=50000;gens=[(2**m+1,m+1) for m in range(1,cap.bit_length()) if m!=3 and 2**m+1<=cap]
known={1:0};todo=[1]
while todo:
 x=heapq.heappop(todo)
 for q,w in gens:
  y=x*q
  if y>cap:continue
  cost=known[x]+w
  if y in known:assert known[y]==cost
  else:known[y]=cost;heapq.heappush(todo,y)
for x in range(1,cap+1):
 rem=x;weight=0
 for m in reversed(range(1,x.bit_length())):
  if m==3:continue
  q=2**m+1
  while rem%q==0:rem//=q;weight+=m+1
 assert (rem==1)==(x in known)
 if rem==1:assert weight==known[x]
print(json.dumps({'status':'PASS','direct_permutation_chains':physical,'product_partition_max_weight':N,'greedy_membership_inputs':cap,'monoid_members':len(known)}))
