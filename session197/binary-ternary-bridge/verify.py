import itertools,json
from graphs import require,hp,run_dp,strong,SEED
from bridge import graph,physical,value,inverse,PORTS

def same(A,B,f):return len(A)==len(B)==len(f)and sorted(f)==list(range(len(A)))and all(bool(A[i]>>j&1)==bool(B[f[i]]>>f[j]&1)for i in range(len(A))for j in range(len(A)))
def cover(rows,root,module,C):
 outside=[i for i in range(len(rows))if i not in module];O=[i for i in outside if rows[root]>>i&1];I=[i for i in outside if rows[i]>>root&1]
 require(all(i in C or j in C for i in O for j in I if rows[i]>>j&1),'crossing cover')
P=[2,12,145,149,99,399,431,19,159];require(graph(1,0)==SEED,'boundary source');require(graph(1,1)==P,'positive source')
require(sum(all(P[a]>>b&1 for a,b in zip(w,w[1:]))for w in itertools.permutations(range(9)))==361,'literal positive count')
require(run_dp(SEED,7)[:4]==[59,179,180,60]and not any(run_dp(SEED,7)[4:]),'negative run vector')
require(run_dp(P,7)[:4]==[121,361,360,120]and not any(run_dp(P,7)[4:]),'positive root vector')
require(run_dp(P,8)[:3]==[179,361,180]and not any(run_dp(P,8)[3:]),'bridge root vector')
grid=[]
for m,t in itertools.product(range(1,4),range(1,3)):
 rows=graph(m,t);h=hp(rows);require(h==value(m,t),'tensor grid');require(strong(rows),'grid strongness');cover(rows,7,set(range(7,7+m)),{1,4});cover(rows,7+m,set(range(7+m,len(rows))),{4});grid.append({'m':m,'t':t,'H':h})
words=[]
for L in range(5):
 for w in itertools.product('RL',repeat=L):
  word=''.join(w);rows,ports,f=physical(word);m=1+word.count('R');t=word.count('L');require(ports==PORTS,'persistent ports');require(same(graph(m,t),rows,f),'physical relabeling');h=hp(rows);require(h==value(m,t),'physical count');require(strong(rows),'physical strongness');words.append({'word':word,'H':h})
extra=[(4,0),(4,1),(1,4),(3,3),(5,2)]
for m,t in extra:require(hp(graph(m,t))==value(m,t),'off-grid count')
inverses=[]
for N in range(1,100001):
 pair=inverse(N)
 if pair:
  require(value(*pair)==N,'inverse soundness');inverses.append({'N':N,'m':pair[0],'t':pair[1]})
for m,t in itertools.product(range(1,9),repeat=2):require(inverse(value(m,t))==(m,t),'inverse completeness control')
print(json.dumps({'status':'PASS','grid':grid,'word_checks':words,'off_grid':extra,'inverse_targets_to100000':inverses,'literal_permutations':362880},sort_keys=True,indent=2))
