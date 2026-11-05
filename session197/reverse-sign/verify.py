import itertools,json
from graphs import SEED,birth,clone,run_dp,hp,strong,require
P=[2,12,145,149,99,399,431,19,159]
N=[2,12,145,149,611,911,431,19,671,207]
def dfs(rows):
 n=len(rows);full=(1<<n)-1
 def rec(last,used):
  if used==full:return 1
  todo=rows[last]&~used;ans=0
  while todo:
   b=todo&-todo;todo-=b;ans+=rec(b.bit_length()-1,used|b)
  return ans
 return sum(rec(v,1<<v)for v in range(n))
def same(A,B,f):return len(A)==len(B)==len(f)and sorted(f)==list(range(len(A)))and all(bool(A[i]>>j&1)==bool(B[f[i]]>>f[j]&1)for i in range(len(A))for j in range(len(A)))
rows=SEED[:];ports=[0,1,4,5,6]
rows,ports=birth(rows,ports,'L',24,4);require(rows==P and ports==[0,1,4,5,8],'first gate')
rows,ports=birth(rows,ports,'L',28,5);require(rows==N and ports==[0,1,4,5,8],'second gate')
cores=[]
for T,a0,A,beta in [(SEED,59,60,-1),(P,121,120,1),(N,419,420,-1)]:
 a=run_dp(T,7);require(a[:4]==[a0,3*A+beta,3*A,A]and not any(a[4:]),'complete run vector');require(dfs(T)==hp(T)==a[1],'independent source count');require(strong(T),'core strongness');O=[i for i in range(len(T))if T[7]>>i&1];I=[i for i in range(len(T))if T[i]>>7&1];require(all(i in {1,4}or j in {1,4}for i in O for j in I if T[i]>>j&1),'crossing cover');cores.append({'rows':T,'runs':a[:4]})
controls=[]
for r0,r1,r2 in itertools.product(range(4),repeat=3):
 if r0+r1+r2>3:continue
 rows=SEED[:];ports=[0,1,4,5,6];roots=[7];gates=[]
 for phase,r in enumerate([r0,r1,r2]):
  for _ in range(r):roots.append(len(rows));rows,ports=birth(rows,ports,'R',7,5)
  if phase<2:
   gates.append(len(rows));rows,ports=birth(rows,ports,'L',24 if phase==0 else 28,4 if phase==0 else 5)
  A,beta=[(120,1),(420,-1),(420,-1)][phase];h=hp(rows);require(h==A*3**len(roots)+beta,'phase count');require(strong(rows),'phase strongness')
 m=len(roots);canonical=clone(N,7,m);mapping=list(range(7))+gates+roots;require(same(canonical,rows,mapping),'all-interleaving clone identity');require(ports==[0,1,4,5,gates[0]],'final exposed gate');controls.append({'root_births':[r0,r1,r2],'H':hp(rows),'ports':ports})
rows=N[:];ports=[0,1,4,5,8]
rows,ports=birth(rows,ports,'L',24,4);rows,ports=birth(rows,ports,'L',28,5);a=run_dp(rows,7);require(a[:4]==[3607,11071,11364,3900]and not any(a[4:]),'repeat obstruction');require(hp(rows)==11071!=8819,'repeat scalar obstruction')
trans=[30,28,24,16,0];ports=list(range(5));trans,ports=birth(trans,ports,'L',24,4);trans,ports=birth(trans,ports,'L',28,5);require(trans==[30,28,88,112,32,71,19]and hp(trans)==dfs(trans)==121,'arbitrary-state hostile')
print(json.dumps({'status':'PASS','cores':cores,'interleaving_controls':controls,'repeated_word_count':11071,'repeated_word_modes':[3900,-336,43],'transitive_hostile_count':121},sort_keys=True,indent=2))
