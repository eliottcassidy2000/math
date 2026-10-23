"""Independent fixed-start Hamiltonian audit of the n=15 E-signature hostile.
Recursive completion counts are multiplied by the transitive cyclic action.
"""
import functools,itertools,json,pathlib
P=pathlib.Path(__file__).resolve().parent
sets=((14,13,12,11,5,6,8),(14,13,12,4,10,9,8));result=[]
def need(q,m):
 if not q:raise ValueError(m)
def measurements(A):
 n=len(A);c=[0]*n;e=[[0]*n for _ in range(n)];v4=[0]*n
 for i,j,k in itertools.combinations(range(n),3):
  if A[i][j]==A[j][k]==A[k][i]:
   for v in (i,j,k):c[v]+=1
   for u,v in itertools.combinations((i,j,k),2):e[u][v]+=1;e[v][u]+=1
 for S in itertools.combinations(range(n),4):
  if sorted(sum(A[u][v] for v in S) for u in S)==[1,1,2,2]:
   for v in S:v4[v]+=1
 E=tuple(sorted(tuple(sorted((c[v],e[u][v]) for v in range(n) if A[u][v])) for u in range(n)))
 tags=[(sum(A[v]),c[v],v4[v]) for v in range(n)]
 return c,e,v4,E,tags
signatures=[]
for diffs in sets:
 A=[[int((v-u)%15 in diffs) for v in range(15)] for u in range(15)]
 need(all(A[u][v]+A[v][u]==1 for u in range(15) for v in range(u+1,15)),'tournament')
 out=[sum(A[u][v]<<v for v in range(15)) for u in range(15)]
 @functools.lru_cache(None)
 def complete(last,remaining):
  if not remaining:return 1
  m=out[last]&remaining;s=0
  while m:
   b=m&-m;m-=b;s+=complete(b.bit_length()-1,remaining^b)
  return s
 fixed=complete(0,((1<<15)-1)^1);H=15*fixed
 subset_H=[0]*(1<<15)
 for S in range(1,1<<15):
  subset_H[S]=sum(complete(u,S^(1<<u)) for u in range(15) if S>>u&1)
 need(subset_H[-1]==H,'full-start versus fixed-start DP')
 F2=sum(subset_H[S]*subset_H[((1<<15)-1)^S] for S in range(1,(1<<15)-1))
 c,e,v4,E,tags=measurements(A)
 need(set(c)=={28} and set(v4)=={184},'uniform incidence')
 need(all(sorted(e[u][v] for v in range(15) if A[u][v])==[3,3,3,4,4,5,6] for u in range(15)),'uniform edge valency')
 need(all(A[(u+1)%15][(v+1)%15]==A[u][v] for u in range(15) for v in range(15)),'translation automorphism')
 # Same initial vertex tag and same multiset weighted neighbors implies stable
 # equality at every ordinary or edge-weighted one-dimensional refinement round.
 signatures.append(E)
 result.append({'differences':list(diffs),'fixed_start_paths':fixed,'total_paths':H,'ordered_two_block_path_covers':F2,'memo_states':complete.cache_info().currsize,'vertex_tags':sorted(set(tags)),'edge_incidence_multiset':[3,3,3,4,4,5,6]})
need(signatures[0]==signatures[1],'exact E collision')
need([r['total_paths'] for r in result]==[196264905,196277235],'independent H values')
# Explicit padding controls at n16..18; proof covers every nonnegative padding.
padding=[]
for k in range(1,4):
 ss=[]
 for diffs in sets:
  n=15+k;A=[[int((v-u)%15 in diffs) if u<15 and v<15 else int(u<v) for v in range(n)] for u in range(n)]
  ss.append(measurements(A)[3])
 need(ss[0]==ss[1],'sink-padding E collision');padding.append(n)
O={'status':'PASS','method':'fixed-start recursive completions times cyclic transitivity, independent of full-start forward producer DP','circulants':result,'E_equal':True,'nonisomorphic_by_unequal_H':True,'all_order_refinement_obstruction':'Identical uniform rooted tags and outgoing edge-weight multiplicities remain identical at every one-dimensional refinement round. Adjoin an ordered chain of universal sinks to preserve E equality and H inequality at all n>=15.','explicit_padding_tests':padding}
(P/'circulant_audit_receipt.json').write_text(json.dumps(O,indent=2)+'\n');print(json.dumps(O,indent=2))
