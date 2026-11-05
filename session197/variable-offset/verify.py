import itertools,json
from graphs import require,hp,run_dp,strong,birth,SEED
from core import *
def same(A,B,f):return len(A)==len(B)==len(f)and sorted(f)==list(range(len(A)))and all(bool(A[i]>>j&1)==bool(B[f[i]]>>f[j]&1)for i in range(len(A))for j in range(len(A)))
def dfs(rows):
 full=(1<<len(rows))-1
 def rec(v,M):
  if M==full:return 1
  todo=rows[v]&~M;ans=0
  while todo:b=todo&-todo;todo-=b;ans+=rec(b.bit_length()-1,M|b)
  return ans
 return sum(rec(v,1<<v)for v in range(len(rows)))
require(canonical(1,0)==SEED_B,'literal initial gate');require(dfs(SEED_B)==537,'independent source count')
a=run_dp(SEED_B,7);require(a[:4]==[177,537,540,180]and not any(a[4:]),'boundary runs')
grid=[];vectors=[]
for t in range(1,4):
 rows=canonical(1,t);a=run_dp(rows,7);A,C=coefficient(t);require(a[1:4]==[3*A+C,3*A,A]and not any(a[4:]),'complete root vectors');vectors.append({'t':t,'runs':a[:4]})
 for m in range(1,4):
  rows=canonical(m,t);h=hp(rows);require(h==value(m,t),'tensor grid');require(strong(rows),'grid strongness')
  for r,module,cover in [(7,set(range(7,7+m)),{1,4}),(8+m,set(range(8+m,len(rows))),{4,7+m})]:
   outside=[v for v in range(len(rows))if v not in module];O=[v for v in outside if rows[r]>>v&1];I=[v for v in outside if rows[v]>>r&1];require(all(x in cover or y in cover for x in O for y in I if rows[x]>>y&1),'uniform cover')
  grid.append({'m':m,'t':t,'H':h})
words=[]
for length in range(5):
 for tup in itertools.product('RL',repeat=length):
  word=''.join(tup);m=1+word.count('R');t=word.count('L');rows,ports,f=physical(word);require(ports==PORTS and same(canonical(m,t),rows,f),'interleaving graph');require(hp(rows)==value(m,t),'interleaving count');require(strong(rows),'interleaving strongness')
  state=(537,3,3,2)
  for c in word:
   h,x,z,d=state;state=(2*h+40*x+20*z+5,3*x,z,2*d)if c=='L'else(3*h-2*d+10,3*x,3*z,d)
  require(state==(value(m,t),3**(m+t),3**m,2**(t+1)),'four-coordinate transfer');words.append({'word':word,'H':state[0]})
N=[2,12,145,149,611,911,431,19,671,207];require(same(N,canonical(1,1),list(range(8))+[9,8]),'reverse-core gate swap')
rows,_=birth(SEED_B,PORTS,'L',12,5);require(hp(rows)==1483!=1611,'nonreturning initial gate')
for m,t in itertools.product(range(1,7),range(7)):require((m,t)in inverse(value(m,t)),'inverse completeness control')
solutions=[]
for n in range(1,10001):
 for m,t in inverse(n):require(value(m,t)==n,'inverse soundness');solutions.append([n,m,t])
print(json.dumps({'status':'PASS','grid':grid,'root_vectors':vectors,'word_checks':words,'inverse_solutions_to10000':solutions,'repeat_B_hostile':1483},sort_keys=True,indent=2))
