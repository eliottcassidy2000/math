"""Independent literal injection, quantified skew and actual quartic batteries."""
from pathlib import Path
from itertools import combinations,permutations,product
from functools import lru_cache
import json,hashlib,subprocess,sys
HERE=Path(__file__).resolve().parent;SOURCE=HERE.parent/'ordered_module_skew'
def need(x,msg):
 if not x:raise RuntimeError(msg)
def rows(n,code):
 R=[0]*n
 for k,(u,v) in enumerate(combinations(range(n),2)):
  if code>>k&1:R[u]|=1<<v
  else:R[v]|=1<<u
 return tuple(R)
def op(R):return tuple(((1<<len(R))-1)^r^(1<<u) for u,r in enumerate(R))
def induce(R,vs):return tuple(sum(1<<j for j,v in enumerate(vs) if R[u]>>v&1) for u in vs)
def mod(R,M):return all(r&M in(0,M) for u,r in enumerate(R) if not M>>u&1)
def literalH(R):return sum(all(R[u]>>v&1 for u,v in zip(p,p[1:])) for p in permutations(range(len(R))))
def covers(R):
 out=[]
 for p in permutations(range(len(R))):
  for k in range(1,len(R)):
   A,B=p[:k],p[k:]
   if all(R[u]>>v&1 for path in(A,B) for u,v in zip(path,path[1:])):out.append((A,B))
 return out
@lru_cache(None)
def endpoint(R):
 n=len(R);d=[[0]*n for _ in range(1<<n)]
 for u in range(n):d[1<<u][u]=1
 for used in range(1,1<<n):
  for u,z in enumerate(d[used]):
   if not z:continue
   rem=R[u]&~used
   while rem:
    bit=rem&-rem;rem-=bit;d[used|bit][bit.bit_length()-1]+=z
 return d
def hp(R):return sum(endpoint(tuple(R))[-1])
def attach(R,a,c):
 n=len(R);T=list(R)+[1<<(n+1),0]
 for u in range(n):
  for mask,v in ((a,n),(c,n+1)):
   if mask>>u&1:T[u]|=1<<v
   else:T[v]|=1<<u
 return tuple(T)
def main():
 paths=[SOURCE/'PROOF.md',SOURCE/'verify.py',HERE.parent/'marked_pair_quotient'/'PROOF.md',HERE.parent/'marked_pair_quotient'/'verify.py'];pins={str(p.relative_to(HERE.parent)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths};data=json.loads((SOURCE/'receipt.json').read_text())
 count=injections=0
 for n in range(2,5):
  for code in range(1<<(n*(n-1)//2)):
   R=rows(n,code);cs=covers(R)
   for choices in product(range(3),repeat=n):
    A={u for u in range(n) if choices[u]==1};B={u for u in range(n) if choices[u]==2};M=A|B;mask=sum(1<<u for u in M)
    if not A or not B or not mod(R,mask) or not all(R[u]>>v&1 for u in A for v in B):continue
    left=[(P,Q) for P,Q in cs if P[-1] in A and Q[0] in B];right={(P,Q) for P,Q in cs if P[-1] in B and Q[0] in A};image=set()
    for P,Q in left:
     i=len(P)
     while i and P[i-1] in M:i-=1
     j=0
     while j<len(Q) and Q[j] in M:j+=1
     need(set(P[i:])<=A and set(Q[:j])<=B,'pure endpoint runs');v=(P[:i]+Q[:j],P[i:]+Q[j:]);need(v in right and v not in image,'literal injective switch');image.add(v);injections+=1
    AA,BB=induce(R,sorted(A)),induce(R,sorted(B));K=induce(R,[min(M)]+[u for u in range(n) if u not in M]);lower=literalH(K)*(len(covers(AA))*literalH(BB)+literalH(AA)*len(covers(BB))+len(covers(AA))*len(covers(BB)))
    need(len(right)-len(left)>=lower and ((len(right)>len(left))==(len(M)>=3)),'quantified bound and pair boundary');count+=1
 # Every actual quartic battery and its independently computed endpoint-cover skew.
 batteries=queries=0
 for rec in data['literal_prime_five_batteries']+data['larger_cyclic_prime_batteries']:
  T,U=tuple(rec['rows_T']),tuple(rec['rows_U']);z=rec['pair_index'];sizes=list(map(len,rec['children']));offs=[sum(sizes[:i]) for i in range(len(sizes)+1)];removed=set(range(offs[z],offs[z+1]));keep=[u for u in range(len(T)) if u not in removed];R=induce(T,keep);index={u:i for i,u in enumerate(keep)}
  A,B=[{index[u] for child in group for u in range(offs[child],offs[child+1])} for group in rec['marked_child_sets']];D=endpoint(R);E=endpoint(op(R));N=(1<<len(R))-1;skew=0
  for used in range(1,N):skew+=sum(D[used][u] for u in B)*sum(E[N^used][v] for v in A)-sum(D[used][u] for u in A)*sum(E[N^used][v] for v in B)
  need(skew==rec['cover_skew']==rec['quartic_coefficient'] and skew>0,'independent full cover skew');total=0
  for q in rec['literal_battery']:
   a,b=hp(attach(T,q['a'],q['c'])),hp(attach(U,q['a'],q['c']));need(a==q['T'] and b==q['U'],'actual named finite-difference query');total+=q['sign']*(a-b);queries+=1
  need(total==skew,'exact grouped quartic');batteries+=1
 cyclic=0
 for rec in data['cyclic_prime_structure']:
  R=tuple(rec['quotient']);n=len(R);r=(n-1)//2;need(R==tuple(sum(1<<((u+d)%n) for d in range(1,r+1)) for u in range(n)),'cyclic definition')
  need(all(not mod(R,M) for M in range(1,(1<<n)-1) if M.bit_count()>1),'independent prime check')
  for z in range(n):
   names=[u for u in range(n) if u!=z];V=induce(R,names);i,j=names.index((z+r)%n),names.index((z+r+1)%n);need(mod(V,(1<<i)|(1<<j)) and V[i]>>j&1,'deleted opposite pair');cyclic+=1
 normal=subprocess.check_output([sys.executable,str(SOURCE/'verify.py')]);optimized=subprocess.check_output([sys.executable,'-O',str(SOURCE/'verify.py')]);need(normal==optimized and json.loads(normal)==data,'source receipt replay')
 need(pins=={str(p.relative_to(HERE.parent)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},'source freeze')
 result={'status':'PASS','independent_ordered_module_cases_through4':count,'independent_literal_cover_injections':injections,'independent_grouped_quartic_batteries':batteries,'independent_actual_query_values_per_graph':queries,'independent_cyclic_deleted_pair_checks':cyclic,'producer_normal_optimized_equal':True,'producer_receipt_sha256':hashlib.sha256(normal).hexdigest(),'source_sha256':pins,'audit_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
 (HERE/'receipt.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,sort_keys=True))
if __name__=='__main__':main()
