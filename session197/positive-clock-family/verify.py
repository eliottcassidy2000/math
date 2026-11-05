import itertools,json
from graphs import hp,run_dp,literal_runs,strong,require
from family import Q,graph,coefficients,value,returning

def same(A,B,f):
 return len(A)==len(B)==len(f) and sorted(f)==list(range(len(A))) and all(bool(A[i]>>j&1)==bool(B[f[i]]>>f[j]&1)for i in range(len(A))for j in range(len(A)))

def covers(m,p,q,u,v):
 lengths=[m,p,u,1,q,1,v];starts=[sum(lengths[:i])for i in range(7)];rows=graph(m,p,q,u,v);cut=[{starts[3],starts[5]},{starts[5]},{starts[5]},None,{starts[3]},None,{starts[3]}]
 for b in [0,1,2,4,6]:
  r=starts[b];outside=[x for x in range(len(rows))if not starts[b]<=x<starts[b]+lengths[b]]
  O=[x for x in outside if rows[r]>>x&1];I=[x for x in outside if rows[x]>>r&1]
  require(all(x in cut[b]or y in cut[b]for x in O for y in I if rows[x]>>y&1),'uniform crossing cover')

require(all(Q[(0,2,1,5,3,4,6)[i]]>>(0,2,1,5,3,4,6)[(i+1)%7]&1 for i in range(7)),'quotient cycle')
source=[30,64,10,66,174,335,21,367,95]
require(same(graph(1,2,2,1,1),source,[0,3,1,2,4,7,5,6,8]),'source relabeling')
require(sum(all(source[a]>>b&1 for a,b in zip(w,w[1:]))for w in itertools.permutations(range(9)))==217,'literal source count')
runs=[];grid=[]
for p,q,u,v in itertools.product([1,2],repeat=4):
 a=run_dp(graph(1,p,q,u,v),0);A,B,C=coefficients(p,q,u,v)
 require(not any(a[4:]),'root run cutoff');require((a[3],a[2]-3*a[3],a[1]-2*a[2]+3*a[3])==(A,B,C),'run coefficient identity')
 runs.append({'parameters':[p,q,u,v],'runs':a[:4],'coefficients':[A,B,C]})
 for m in [1,2,3]:
  rows=graph(m,p,q,u,v);h=hp(rows);require(h==value(m,p,q,u,v),'grid formula');require(strong(rows),'grid strongness');covers(m,p,q,u,v)
  grid.append({'parameters':[m,p,q,u,v],'H':h})
extra=[(4,1,2,1,1),(2,3,1,1,2),(1,2,3,2,1),(2,1,1,3,2),(3,2,2,1,2),(1,3,3,2,2),(2,2,1,2,3)]
for t in extra:require(hp(graph(*t))==value(*t),'off-grid formula');covers(*t)
returns=[]
for v,m in itertools.product(range(1,4),range(1,5)):
 rows,ports,f=returning(m,v);require(ports==[1,2,3,4],'persistent ports');require(same(graph(m,2,2,1,v),rows,f),'return graph identity');h=hp(rows);require(h==72*2**(v-1)*3**m+1,'return count');require(strong(rows),'return strongness')
 returns.append({'m':m,'v':v,'order':len(rows),'H':h})
for p,q,u,v in itertools.product(range(1,9),repeat=4):require((coefficients(p,q,u,v)[1]==0)==(p==q==2),'cancellation control')
print(json.dumps({'status':'PASS','grid':grid,'runs':runs,'off_grid':len(extra),'literal_source_permutations':362880,'returns':returns,'cancellation_controls':4096},sort_keys=True,indent=2))
