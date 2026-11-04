import itertools,json,sys,math,hashlib
from pathlib import Path
from functools import lru_cache
P=Path('/workspace/shared/session187/gate');raw=json.loads((P/'full_response.json').read_text());receipt=json.loads((P/'receipt.json').read_text())
def need(x,m):
 if not x:raise RuntimeError(m)
def hp(a):
 n=len(a);pred=[sum(1<<u for u in range(n) if a[u]>>v&1) for v in range(n)]
 @lru_cache(None)
 def e(mask,v):
  rem=mask^(1<<v)
  if not rem:return 1
  b=rem&pred[v];s=0
  while b:
   bit=b&-b;b-=bit;s+=e(rem,bit.bit_length()-1)
  return s
 return sum(e((1<<n)-1,v) for v in range(n))
def extend(a,x,y):
 n=len(a);out=list(a)+[1<<(n+1),0]
 for u in range(n):
  for v,m in [(n,x),(n+1,y)]:
   if m>>u&1:out[u]|=1<<v
   else:out[v]|=1<<u
 return out
def inflate(seed,r,order=None):
 n=r+5;order=order or list(range(5,n));out=[0]*n
 for u,v in itertools.combinations(range(n),2):
  if u<5 and v<5:forward=seed[u]>>v&1
  elif u<5:forward=seed[u]>>5&1
  else:forward=order.index(u)<order.index(v)
  if forward:out[u]|=1<<v
  else:out[v]|=1<<u
 return out
def reverse(a):return [((1<<len(a))-1)^r^(1<<u) for u,r in enumerate(a)]
seeds=[[16,5,33,39,46,3],[16,5,49,39,42,3],[18,4,49,39,42,3]]
def skeleton(a,root,ignore=None,deleted=(),boundary=None):
 names=[i for i in range(len(a)) if i!=root and i not in deleted];hist=[0]*(len(names)+2)
 for p in itertools.permutations(names):
  for mask in range(1<<(len(names)+1)):
   word=[]
   for pos in range(len(names)+1):
    if mask>>pos&1:word.append(root)
    if pos<len(names):word.append(p[pos])
   if boundary and not any((u,v)==boundary for u,v in zip(word,word[1:])):continue
   good=all(u==ignore or (boundary and (u,v)==boundary) or a[u]>>v&1 for u,v in zip(word,word[1:]))
   if good:hist[mask.bit_count()]+=1
 return hist
profiles=0
for j,seed in enumerate(seeds):
 row=receipt['profiles'][j];a=inflate(seed,2)
 checks=[('first_score_profile',skeleton(a,6,5)),('exterior0_score_profile',skeleton(seed,5,0)),('op_exterior0_score_profile',skeleton(reverse(seed),5,0)),('quartic_forward_profile',skeleton(a,6,deleted=(0,),boundary=(5,3))),('quartic_backward_profile',skeleton(a,6,deleted=(0,),boundary=(3,5)))]
 for key,v in checks:need(v==row[key],('independent explicit gap subsets',j,key,v));profiles+=1
seen=set();count=0;checks=0
for case in raw['cases']:
 j,r,kind=case['state'],case['r'],case['kind'];need((j,r,kind) not in seen,'unique packet');seen.add((j,r,kind));a=inflate(seeds[j],r);b=reverse(a)
 if kind=='pair_reversal':b=list(a);b[5]^=1<<6;b[6]^=1<<5
 need(a==case['rows'] and b==case['other_rows'],'independent physical rows');ms=[m for m in range(1<<len(a)) if m.bit_count()<=2];need(ms==case['battery_masks'],'full masks');count+=2*len(case['battery']);N=len(ms)
 # Diverse named queries: zeros, singleton, quartic-support, outer boundary.
 for x,y in [(0,0),(1,0),(32,0),(9,65 if r>=2 else 33),(ms[-1],ms[-2])]:
  pair=case['battery'][ms.index(x)*N+ms.index(y)]
  need([hp(extend(a,x,y)),hp(extend(b,x,y))]==pair,('fresh actual count',j,r,kind,x,y));checks+=2
 for x,y,d in case['dense_controls']:
  need(hp(extend(a,x,y))-hp(extend(b,x,y))==d,'all high-weight controls');checks+=2
 need(len(case['coefficients'])==len({(x,y) for x,y,z in case['coefficients']}),'unique coefficient names')
expect={(j,r,k) for j in range(3) for r in range(1,8) for k in ['global_converse','pair_reversal'] if r>=2 or k!='pair_reversal'}
need(seen==expect and count==221208,'packet universe and retained count')
# Every physical module order for r=3, all six external seeds, tested directly.
family=0
for j,seed in enumerate(seeds):
 for orient in range(2):
  root=seed if not orient else reverse(seed)
  for order in itertools.permutations(range(5,8)):
   a=inflate(root,3,list(order));vals={v:hp(extend(a,1<<v,0)) for v in order};need(tuple(sorted(order,key=vals.get))==order,'rank reconstruction');family+=1
   H=hp(a);need(H==4*27+9*8-3+2*j,'state recovery')
   expected=[8*27+15*8-9,14*27+21*8-9,14*27+9*8-6][j] if not orient else [4*27+9*8-3,4*27+9*8,4*27+12*8+3][j]
   need(hp(extend(a,1,0))==H+2*expected,'external orientation recovery');checks+=5
# Literal Hamilton-word expansion, no endpoint tensors or interpolation.
@lru_cache(None)
def wordpoly(a):
 n=len(a);poly={}
 for word in itertools.permutations(range(n+2)):
  factors=[];good=True
  for u,v in zip(word,word[1:]):
   if u<n and v<n:
    if not a[u]>>v&1:good=False;break
   elif u>=n and v>=n:
    if u!=n:good=False;break
   elif u<n:factors.append((u+(n if v==n+1 else 0),True))
   else:factors.append((v+(n if u==n+1 else 0),False))
  if not good:continue
  terms={0:1}
  for bit,positive in factors:
   new={}
   for mask,c in terms.items():
    if not positive:new[mask]=new.get(mask,0)+c
    key=mask|(1<<bit);new[key]=new.get(key,0)+(c if positive else -c)
   terms=new
  for key,c in terms.items():poly[key]=poly.get(key,0)+c
 return {key:c for key,c in poly.items() if c}
word_cases=0
for case in raw['cases']:
 if case['r']>2:continue
 p=wordpoly(tuple(case['rows']));q=wordpoly(tuple(case['other_rows']));n=len(case['rows'])
 diff={(k&((1<<n)-1),k>>n):p.get(k,0)-q.get(k,0) for k in p.keys()|q.keys() if p.get(k,0)!=q.get(k,0)}
 need(diff=={(x,y):v for x,y,v in case['coefficients']},'literal entire polynomial');word_cases+=1
out={'status':'PASS','literal_full_polynomial_cases_r1_r2':word_cases,'independently_enumerated_gap_profiles':profiles,'verified_unique_packets':len(seen),'retained_paired_battery_counts':count,'fresh_predecessor_DP_counts':checks,'named_reordering_family_cases_r3':family,'full_response_sha256':hashlib.sha256((P/'full_response.json').read_bytes()).hexdigest(),'scope':'Independent direct gap-subset enumeration and predecessor DP; full coefficient inversion separately author replayed, not independently every large battery value.'}
print(json.dumps(out,sort_keys=True,indent=2))
