"""Exact highest-changed-quotient reduction controls."""
from itertools import combinations,permutations
from math import comb,prod
from functools import lru_cache
import json

def need(ok,msg):
 if not ok:raise RuntimeError(msg)

def rows(n,code):
 a=[0]*n
 for bit,(i,j) in enumerate(combinations(range(n),2)):
  u,v=(j,i) if code>>bit&1 else(i,j);a[u]|=1<<v
 return tuple(a)

def induced(a,mask):
 names=[i for i in range(len(a)) if mask>>i&1]
 return tuple(sum(1<<j for j,v in enumerate(names) if a[u]>>v&1) for u in names)

def op(a):return tuple(((1<<len(a))-1)^r^(1<<i) for i,r in enumerate(a))

def module(a,M):return all(r&M in(0,M) for i,r in enumerate(a) if not M>>i&1)

def flip(a,M):return tuple(r^(M^(1<<i)) if M>>i&1 else r for i,r in enumerate(a))

@lru_cache(None)
def q_words(a):
 n=len(a)
 if not n:return (1,)
 q=[0]*(n+1)
 for w in permutations(range(n)):
  bad=sum(not(a[u]>>v&1) for u,v in zip(w,w[1:]));good=n-1-bad
  for k in range(bad+1,n+1):q[k]+=comb(good,k-1-bad)
 return tuple(q)

@lru_cache(None)
def hp(a):
 n=len(a);full=(1<<n)-1
 @lru_cache(None)
 def finish(mask,last):
  if mask==full:return 1
  out=0;rem=a[last]&~mask
  while rem:
   b=rem&-rem;rem-=b;out+=finish(mask|b,b.bit_length()-1)
  return out
 return sum(finish(1<<i,i) for i in range(n)) if n else 1

def blow(Q,blocks):
 offsets=[0]
 for b in blocks:offsets.append(offsets[-1]+len(b))
 a=[0]*offsets[-1]
 for i,b in enumerate(blocks):
  for u,r in enumerate(b):
   a[offsets[i]+u]=r<<offsets[i]
   for j,c in enumerate(blocks):
    if Q[i]>>j&1:a[offsets[i]+u]|=((1<<len(c))-1)<<offsets[j]
 return tuple(a)

def quotient(a,blocks):
 reps=[(M&-M).bit_length()-1 for M in blocks]
 return tuple(sum(1<<j for j,v in enumerate(reps) if a[u]>>v&1) for u in reps)

def extend(a,x,c):
 n=len(a);b=list(a)+[1<<(n+1),0]
 for v in range(n):
  if x>>v&1:b[v]|=1<<n
  else:b[n]|=1<<v
  if c>>v&1:b[v]|=1<<(n+1)
  else:b[n+1]|=1<<v
 return tuple(b)


def linear(a):return sorted(r.bit_count() for r in a)==list(range(len(a)))

def tree(a):
 mods=[M for M in range(1,1<<len(a)) if module(a,M)]
 strong=[M for M in mods if all(M&N in(0,M,N) for N in mods)]
 nodes={}
 for M in strong:
  if M.bit_count()<2:continue
  inside=[N for N in strong if N!=M and N&M==N]
  children=[N for N in inside if not any(N!=K and N&K==N for K in inside)]
  nodes[M]=tuple(children)
 return mods,nodes

def highest(T,U,nodes):
 changed=[M for M,ch in nodes.items() if quotient(T,ch)!=quotient(U,ch)]
 return [M for M in changed if not any(M!=N and M&N==M for N in changed)]

def check_pair(T,U,nodes):
 tops=highest(T,U,nodes)
 need(bool(tops)==(T!=U),'highest existence')
 need(all(A&B==0 for A,B in combinations(tops,2)),'highest antichain')
 covered=0
 for M in tops:covered|=M
 for u,v in combinations(range(len(T)),2):
  if (T[u]>>v&1)!=(U[u]>>v&1):
   need(any(M>>u&1 and M>>v&1 for M in tops),'changed edge outside top interiors')
 records=[]
 for M in tops:
  ch=nodes[M];Q,R=quotient(T,ch),quotient(U,ch)
  need(all(q_words(induced(T,A))==q_words(induced(U,A)) for A in ch),'hereditary child profile equality')
  if not linear(Q):
   need(R==op(Q),'changed prime not converse');records.append({'module':M,'kind':'prime','quotient_order':len(Q)})
   continue
  need(linear(R),'scalar equal linear quotient became cyclic')
  first=sorted(range(len(ch)),key=lambda i:-Q[i].bit_count())
  second=sorted(range(len(ch)),key=lambda i:-R[i].bit_count())
  t=next(t for t in range(1,len(ch)) if set(first[:t])!=set(second[:t]))
  prefix=sum(ch[i] for i in first[:t]);b=(M&-M).bit_length()-1
  outside=sum(1<<v for v in range(len(T)) if not M>>v&1 and T[v]>>b&1)
  a,c=outside,outside|prefix
  ambient=[M]+[1<<v for v in range(len(T)) if not M>>v&1]
  lower=2*hp(quotient(T,ambient))*prod(hp(induced(T,A)) for A in ch)
  ft,fu=hp(extend(T,a,c)),hp(extend(U,a,c))
  need(fu-ft>=lower,'highest linear physical separator')
  records.append({'module':M,'kind':'linear','quotient_order':len(Q),'query':[a,c],'values':[ft,fu],'lower_bound':lower})
 return records

def main():
 stats={'families':0,'equal_families':0,'equal_nonidentity':0,'linear_highest_nodes':0,'prime_highest_nodes':0,'multiple_highest_families':0,'prime_only_families':0};census=[]
 for n in range(1,6):
  local={k:0 for k in stats}
  for code in range(1<<(n*(n-1)//2)):
   T=rows(n,code);mods,nodes=tree(T);mods=[M for M in mods if M.bit_count()>1]
   for bits in range(1<<len(mods)):
    U=T
    for j,M in enumerate(mods):
     if bits>>j&1:U=flip(U,M)
    local['families']+=1
    if hp(T)!=hp(U):continue
    local['equal_families']+=1
    records=check_pair(T,U,nodes)
    if not records:continue
    local['equal_nonidentity']+=1
    lin=sum(r['kind']=='linear' for r in records);pr=len(records)-lin
    local['linear_highest_nodes']+=lin;local['prime_highest_nodes']+=pr
    local['multiple_highest_families']+=len(records)>1
    local['prime_only_families']+=not lin
  for k,v in local.items():stats[k]+=v
  census.append({'n':n,**local})
 # Two disjoint highest changes, one linear and one prime.
 T=blow((2,4,1),(rows(2,0),(2,4,1),(0,)));U=blow((2,4,1),(op(rows(2,0)),op((2,4,1)),(0,)))
 _,nodes=tree(T);multi=check_pair(T,U,nodes)
 need(len(multi)==2 and {r['kind'] for r in multi}=={'linear','prime'},'multiple control')
 # A changed C3 prime root and changed linear descendants: only root is highest.
 T3=blow((2,4,1),(rows(2,0),)*3);U3=op(T3);_,nodes3=tree(T3);prime3=check_pair(T3,U3,nodes3)
 need(len(prime3)==1 and prime3[0]['kind']=='prime' and prime3[0]['quotient_order']==3,'prime-three root control')
 changed_children=[M for M,ch in nodes3.items() if M!=(1<<len(T3))-1 and quotient(T3,ch)!=quotient(U3,ch)]
 need(len(changed_children)==3,'changed linear descendants control')
 return {'status':'PASS','totals':stats,'by_order':census,
         'multiple_highest_control':{'T':T,'U':U,'nodes':multi},
         'prime_three_control':{'T':T3,'U':U3,'highest':prime3,'changed_descendants':changed_children},
         'scope':'Necessary highest-prime quotient reduction, not full injectivity.'}

if __name__=='__main__':print(json.dumps(main(),indent=2,sort_keys=True))
