#!/usr/bin/env python3
"""Recipe checks and independent canonical strong-unit grammar control."""
import itertools,json
import compiler
# Direct counters copied explicitly here to avoid importing verifier side effects.
def need(ok,msg):
 if not ok:raise RuntimeError(msg)
def profile(a):
 n=len(a);N=1<<n;F=N-1;E=[[0]*n for _ in range(N)];S=[[0]*n for _ in range(N)];H=[0]*N;H[0]=1
 for v in range(n):E[1<<v][v]=S[1<<v][v]=1
 for mask in range(1,N):
  for v in range(n):
   rem=mask^(1<<v)
   if not(mask>>v&1)or not rem:continue
   for u in range(n):
    if rem>>u&1:
     if a[u]>>v&1:E[mask][v]+=E[rem][u]
     else:S[mask][v]+=S[rem][u]
  H[mask]=sum(E[mask])
 C=[0]*N
 for mask in range(N):
  sub=mask
  while True:
   C[mask]+=H[sub]*H[mask^sub]
   if not sub:break
   sub=(sub-1)&mask
 L=sum(E[x][0]*H[F^x]for x in range(N));R=sum(S[x][0]*H[F^x]for x in range(N));T=sum(S[x][0]*C[F^x]for x in range(N))
 return(H[F],S[F][0],E[F][0],L,R,H[F^1],C[F],C[F^1],T)
def canonical_child(U,w,bits):
 # New root first, then old unique-rooted Hamilton order, then increasing tail.
 m=len(U);n=m+w+1;a=[0]*n
 for i,row in enumerate(U):
  a[i+1]=row<<1
  for j in range(m+1,n):a[i+1]|=1<<j
  if i==0:a[0]|=1<<(i+1)
  else:a[i+1]|=1
 for k in range(w):
  i=m+1+k
  for j in range(i+1,n):a[i]|=1<<j
  if bits>>k&1:a[0]|=1<<i
  else:a[i]|=1
 return tuple(a)
trials=0
for k in range(6):
 for w in itertools.product('OZJ',repeat=k):
  p=compiler.SEED
  for op in w:p=compiler.move(p,op)
  need(profile(compiler.literal(w))==p,('literal word',w));trials+=1
report,w=compiler.inventory(137);need(report['complete'],'inventory incomplete')
rootone={1:{(0,)},2:set()};support={1};count=0;order_counts={}
for n in range(3,10):
 bank=set()
 for m in range(1,n):
  wlen=n-m-1
  for U in rootone[m]:
   for bits in range(1<<max(0,wlen-1)):bank.add(canonical_child(U,wlen,bits))
 rootone[n]=bank;order_counts[n]=len(bank)
 for a in bank:
  p=profile(a);need(p[1]==1,'canonical source not unit');B=compiler.base(p)
  if B<=137:support.add(B)
  count+=1
need(set(w)==support and 137 not in w,'independent support mismatch')
report139,w139=compiler.inventory(139);need(139 in w139,'positive target missing');p=profile(compiler.literal(w139[139]));need(compiler.base(p)==139 and p[1]==1,'positive target wrong')
limited,_=compiler.inventory(137,max_states=1);need(not limited['complete']and limited['reason']=='RESOURCE_LIMIT','resource guard')
r997,w997=compiler.inventory(997);need(r997['complete']and len(w997)==330,'inventory997 changed')
print(json.dumps(dict(status='PASS',literal_source_words_through_length5=trials,independent_canonical_unit_grammar_order_counts=order_counts,independent_strong_graphs=count,exact_unit_B_support_through137=sorted(support),target137_absent_only_at_start1=True,target139_literal_witness=list(compiler.literal(w139[139])),resource_limit_rejects_absence_claim=True,inventory997_base_count=len(w997),scope='All-order unit-start membership uses the proved order cap and inherited complete grammar; no other starts or F8 occupancy.'),sort_keys=True,indent=2))
