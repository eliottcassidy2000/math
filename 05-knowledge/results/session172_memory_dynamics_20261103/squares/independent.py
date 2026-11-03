import itertools, math, json
checks=0
def require(x,label):
 global checks
 checks+=1
 if not x: raise RuntimeError(label)
def edge(i,j):return tuple(sorted((i,j)))
def hp(V,E):
 return [p for p in itertools.permutations(V) if all(edge(a,b) in E for a,b in zip(p,p[1:]))]
def coefficients(n,E):
 V=tuple(range(n)); full=(1<<n)-1
 banks={m:hp(tuple(i for i in V if m>>i&1),E) for m in range(1,full+1)}
 ends=[sum(p[-1]==i for p in banks[full]) for i in V]
 K=[[0]*n for i in V]
 for m in range(1,full):
  for p in banks[m]:
   for q in banks[full^m]:K[p[-1]][q[0]]+=1
 return ends,K
cases=0
for n in range(1,5):
 pairs=list(itertools.combinations(range(n),2))
 for mask in range(1<<len(pairs)):
  E={e for k,e in enumerate(pairs) if mask>>k&1}; ends,K=coefficients(n,E)
  for i,j in pairs:
   toggled=E^{(i,j)}
   delta=abs(len(hp(tuple(range(n)),toggled))-len(hp(tuple(range(n)),E)))
   require(delta==2*K[i][j],'edge influence')
  for sm in range(1<<n):
   S=[i for i in range(n) if sm>>i&1]
   observed=len(hp(tuple(range(n+1)),E|{(i,n) for i in S}))
   predicted=2*sum(ends[i] for i in S)+sum(K[i][j] for i in S for j in S)
   require(observed==predicted,'all-context polynomial')
   cases+=1
seeds=[{(0,1),(0,2),(0,3),(0,4),(1,2),(1,3),(1,4)}, {(0,1),(0,2),(0,3),(0,4),(1,2),(1,3),(2,3)}]
seed_data=[]
for E in seeds:
 n=5; oriented_component_order_counts=[0]*6
 for p in itertools.permutations(range(n)):
  for cuts in range(1<<(n-1)):
   if all((cuts>>i&1) or edge(p[i],p[i+1]) in E for i in range(n-1)):
    oriented_component_order_counts[cuts.bit_count()+1]+=1
 P=[oriented_component_order_counts[j]//math.factorial(j) for j in range(6)]
 require(all(oriented_component_order_counts[j]%math.factorial(j)==0 for j in range(6)),'component order quotient')
 require(P==[0,12,60,54,14,1],'profile')
 ends,K=coefficients(n,E)
 leaf=len(hp(tuple(range(6)), E|{(1,5)}))
 seed_data.append({'E':ends,'profile':P,'leaf':leaf})
 for u in range(0,10):
  F=sum(P[j]*math.prod(range(u-j+1,u+1)) for j in range(1,6))
  require(F==u*u*(u+1)**2*(u+2),'slot polynomial')
require([x['E'] for x in seed_data]==[[0,0,4,4,4],[0,2,2,2,6]],'endpoint vectors')
require([x['leaf'] for x in seed_data]==[0,4],'read-head')
for E,expected in [(set(),0),({(0,1),(2,3)},48)]:
 ends,K=coefficients(5,E)
 require(not any(ends) and not any(map(any,K)),'one-step blind')
 ext=E|{(i,j) for i in range(5) for j in (5,6)}
 require(len(hp(tuple(range(7)),ext))==expected,'two-step count')
# Every listed arithmetic witness uses precisely the declared vertices, ports, edges.
cert=json.load(open('/tmp/s172_squares/certificate.json'))
for rec in cert['resurrections']:
 n=rec['base_n'];i,j=rec['new_ports']
 require(rec['K']=={22:3,24:10}[n],'two-piece count')
 seen=set()
 for p,q in rec['covers']:
  require(p[-1]==i and q[0]==j,'ports')
  require(sorted(p+q)==list(range(1,n+1)),'spanning disjoint pieces')
  full=p+[n+1]+q
  require(all(math.isqrt(a+b)**2==a+b for a,b in zip(full,full[1:])),'literal square sums')
  seen.add(tuple(full))
 require(len(seen)==rec['K'],'distinct certificates')
require({x['n']:x['oriented_H'] for x in cert['square_prefixes']}=={15:2,16:2,17:2,18:0,22:0,23:6,24:0,25:20},'count transcription')
print(json.dumps({'status':'pass','checks':checks,'exhaustive_contexts_n_le_4':cases,'seed_data':seed_data,'note':'independent subset-pair coefficients versus permutation extension counts; independent permutation-cut profile enumeration; arithmetic witness checks'},sort_keys=True,indent=2))
