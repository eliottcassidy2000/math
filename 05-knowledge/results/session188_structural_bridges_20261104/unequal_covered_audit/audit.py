"""Independent run-word kernel and full named coefficient audit."""
from pathlib import Path
from itertools import product,combinations
from functools import lru_cache
import json,hashlib,subprocess,sys
HERE=Path(__file__).resolve().parent;SOURCE=HERE.parent/'unequal_covered_response'
def need(ok,msg):
 if not ok:raise RuntimeError(msg)
def extend(R,a,c):
 n=len(R);U=list(R)+[1<<(n+1),0]
 for u in range(n):
  for mask,port in ((a,n),(c,n+1)):
   if mask>>u&1:U[u]|=1<<port
   else:U[port]|=1<<u
 return tuple(U)
def hp(R):
 # Backward recurrence on unvisited names, independent of forward producer DP.
 N=(1<<len(R))-1
 @lru_cache(None)
 def f(remaining,last):
  if not remaining:return 1
  available=remaining&R[last];total=0
  while available:
   bit=available&-available;available-=bit
   total+=f(remaining^bit,bit.bit_length()-1)
  return total
 return sum(f(N^(1<<u),u) for u in range(len(R)))
def runwords(Q,a,c,r):
 R=extend(Q,a,c)
 @lru_cache(None)
 def count(left,last):
  if not any(left):return 1
  ans=0
  for v,num in enumerate(left):
   if num and (last==5 or R[last]>>v&1):
    new=list(left);new[v]-=1;ans+=count(tuple(new),v)
  return ans
 return count(tuple(r)+(1,1),5)
def kernel_coeff(r):
 # S numerator listed independently as exponent triples.
 S={(2,2,0):1,(0,2,2):1,(2,0,2):1,
    (2,1,1):-1,(1,2,1):-1,(1,1,2):-1,
    (2,1,0):1,(2,0,1):1,(1,2,0):1,(0,2,1):1,(1,0,2):1,(0,1,2):1,(1,1,1):-6}
 return sum(v*(r[0]-p[0]+1) for p,v in S.items() if r[0]-p[0]==r[1]-p[1]==r[2]-p[2]>=0)
def main():
 pins={n:hashlib.sha256((SOURCE/n).read_bytes()).hexdigest() for n in ('PROOF.md','verify.py')}
 kernels=0
 for r in product(range(4),repeat=3):
  lhs=sum(runwords((2,4,1),1<<i,1<<j,r)-runwords((4,1,2),1<<i,1<<j,r) for i,j in ((0,1),(1,2),(2,0)))
  need(lhs==kernel_coeff(r),'literal run-word rational numerator certificate');kernels+=1
 record=json.loads((SOURCE/'normal.json').read_text());case=record['mixed7'];T=tuple(case['T']);U=tuple(case['U']);n=len(T)
 need(n==7 and T[:2]==(14,12),'literal physical names 0,1 and outgoing block A')
 # Derive XOR independently from the original whole and child sets.
 target=list(T)
 for a,b in combinations(range(n),2):
  if (1+int(a//2==b//2 and a<4 and b<4)+int(a>=4 and b>=4))%2:target[a]^=1<<b;target[b]^=1<<a
 need(tuple(target)==U,'whole-plus-three-child simultaneous XOR')
 need(hp(T)==hp(U)==159,'literal old scalar')
 values={};coeff={};checks=0
 for degree in range(5):
  for variables in combinations(range(2*n),degree):
   mask=sum(1<<v for v in variables);a=mask&((1<<n)-1);c=mask>>n
   value=hp(extend(T,a,c))-hp(extend(U,a,c));values[mask]=value
   proper=(mask-1)&mask;coefficient=value
   if mask:
    while True:
     coefficient-=coeff.get(proper,0)
     if not proper:break
     proper=(proper-1)&mask
   if coefficient:coeff[mask]=coefficient
   checks+=1
 saved={x['mask']:x['coefficient'] for x in case['full_named_contrast']}
 need(coeff==saved,'all mixed7 named coefficients including cubic and quartic')
 for q in case['all_64_uniform_queries']:
  need(hp(extend(T,q['a'],q['c']))==q['T'] and hp(extend(U,q['a'],q['c']))==q['U'],'all actual coarse queries')
 need(sum(q['T']-q['U'] for q in case['three_queries'])==-324,'mixed7 cyclic sum')
 # Replay producer unchanged, including its independent full word expansion.
 normal=subprocess.check_output([sys.executable,str(SOURCE/'verify.py')]);optimized=subprocess.check_output([sys.executable,'-O',str(SOURCE/'verify.py')]);need(normal==optimized,'producer normal optimized agreement');need(json.loads(normal)==record,'producer saved receipt reproducible')
 for name,data in [('producer.normal.stdout',normal),('producer.optimized.stdout',optimized)]: (HERE/name).write_bytes(data)
 need(pins=={n:hashlib.sha256((SOURCE/n).read_bytes()).hexdigest() for n in pins},'source stayed frozen')
 out={'status':'PASS','independent_run_word_coefficients':kernels,'kernel_certificate':'After multiplication by (1-abc)^3, both numerators have each variable degree at most three, so the 4^3 series coefficients determine the identity.','independent_named_queries':checks,'independent_uniform_queries':64,'mixed7_degree_counts':{d:sum(m.bit_count()==d for m in coeff) for d in range(5)},'source_sha256':pins,'audit_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'producer_normal_optimized_equal':True}
 (HERE/'receipt.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,sort_keys=True))
if __name__=='__main__':main()
