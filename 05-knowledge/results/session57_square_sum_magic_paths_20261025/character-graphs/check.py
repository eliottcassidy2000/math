"""Independent tuple-field, forward subset-DP, trace/Newton spectrum checks."""
import sys,json,itertools,collections,functools,pathlib
BASE=pathlib.Path(__file__).parent;sys.path.insert(0,str(BASE/'sources'))
from field_arithmetic import K
LINES=((0,1,2),(3,4,5),(6,7,8),(0,3,6),(1,4,7),(2,5,8),(0,4,8),(2,4,6))
def decode(x,p,n):return tuple((x//p**j)%p for j in range(n))
def matrix(F,c,u,v):return [F.add(c,F.add(F.scale(s,u),F.scale(t,v))) for s,t in ((1,0),(-1,-1),(0,1),(-1,1),(0,0),(1,-1),(0,-1),(1,1),(-1,0))]
def valid(F,m):return len(set(m))==9 and set(m)<=F.sq and len({F.add(F.add(m[a],m[b]),m[c]) for a,b,c in LINES})==1
@functools.lru_cache(None)
def stats(A):
 n=len(A);dp=[[0]*n for _ in range(1<<n)]
 for i in range(n):dp[1<<i][i]=1
 for mask in range(1,1<<n):
  for i in range(n):
   if not dp[mask][i]:continue
   for j in range(n):
    if A[i][j] and not mask>>j&1:dp[mask|1<<j][j]+=dp[mask][i]
 powers=[list(r) for r in A];traces=[]
 for k in range(1,n+1):
  traces.append(sum(powers[i][i] for i in range(n)))
  powers=[[sum(powers[i][h]*A[h][j] for h in range(n)) for j in range(n)] for i in range(n)]
 co=[1]
 for k in range(1,n+1):
  v=-sum(co[k-i]*traces[i-1] for i in range(1,k+1));assert v%k==0;co.append(v//k)
 return sum(map(sum,A))//2,sum(dp[-1])//2,tuple(co)
def adjacency(F,V,z):return tuple(tuple(int(i!=j and F.add(x,y) in F.sq and (z or F.add(x,y)!=F.z)) for j,y in enumerate(V)) for i,x in enumerate(V))
def signature(F,V):return tuple((z,*stats(adjacency(F,V,z))) for z in (True,False))
def expected(g):return tuple((z,g['graphs'][str(z)]['edges'],g['graphs'][str(z)]['hp_unoriented'],tuple(g['graphs'][str(z)]['charpoly_coefficients'])) for z in (True,False))
def check_record(F,g):
 V=[decode(x,F.p,F.n) for x in g['vertices']];assert signature(F,V)==expected(g)
 if 'matrix' in g:assert valid(F,[decode(x,F.p,F.n) for x in g['matrix']])
 for z in (True,False):
  s=g['graphs'][str(z)];A=adjacency(F,V,z);assert sorted(map(sum,A))==s['degrees'];path=s['hp_witness_indices']
  if path:assert sorted(path)==list(range(9)) and all(A[a][b] for a,b in zip(path,path[1:]))
D=json.loads((BASE/'results.json').read_text());total=0
for r in D['fields']:
 F=K(r['p'],r['modulus']);counts=collections.Counter()
 # numeric order independently reconstructed from tuples
 enc=lambda x:sum(v*F.p**j for j,v in enumerate(x))
 for c in (F.z,F.one):
  for u in ([F.one] if c==F.z else F.E):
   for v in F.E:
    m=matrix(F,c,u,v)
    if valid(F,m):counts[signature(F,sorted(m,key=enc))]+=1
 want=collections.Counter({expected(g):g['count'] for g in r['spectral_hp_groups']})
 assert counts==want,(r['q'],counts,want)
 assert sum(counts.values())==r['normalized_magic_count'];total+=sum(counts.values())
 for g in r['spectral_hp_groups']:check_record(F,g)
 if 'first_nine_square_graph' in r:check_record(F,r['first_nine_square_graph'])
 print('PASS',r['q'],flush=True)
for q in (729,2187):
 r=json.loads((BASE/f'witness-{q}.json').read_text());F=K(r['p'],r['modulus']);check_record(F,r)
# independent brute-permutation check on both F17 conventions
r=next(r for r in D['fields'] if r['q']==17);F=K(17,r['modulus']);V=[(x,) for x in r['first_nine_square_graph']['vertices']]
for z in (True,False):
 A=adjacency(F,V,z);brute=sum(all(A[a][b] for a,b in zip(P,P[1:])) for P in itertools.permutations(range(9)))//2
 assert brute==stats(A)[1]
result={'status':'PASS','fields':len(D['fields']),'normalized_magic_grids':total,'large_witnesses':[729,2187],'independent_checks':['tuple finite-field arithmetic','exhaustive normalized grids','forward subset Hamiltonian path DP','Newton trace characteristic polynomial','F17 brute all permutations both zero conventions']}
(BASE/'check-result.json').write_text(json.dumps(result,indent=2)+'\n');print(result)
# Fixed affine-plane set: all 9! arrangements, eight prescribed lines.
plane_results=[]
for q in (81,2187):
 r=next(r for r in D['fields'] if r['q']==81) if q==81 else json.loads((BASE/f'witness-{q}.json').read_text())
 g=r['spectral_hp_groups'][0] if q==81 else r;F=K(r['p'],r['modulus']);V=[decode(x,F.p,F.n) for x in g['vertices']]
 H={tuple(t) for t in itertools.combinations(range(9),3) if F.add(F.add(V[t[0]],V[t[1]]),V[t[2]])==F.z};assert len(H)==12
 count=sum(all(tuple(sorted(P[i] for i in line)) in H for line in LINES) for P in itertools.permutations(range(9)))
 assert count==432;plane_results.append({'q':q,'hyperedges':len(H),'magic_arrangements_of_fixed_set':count,'all_arrangements':362880})
result['affine_plane_permutation_checks']=plane_results
(BASE/'check-result.json').write_text(json.dumps(result,indent=2)+'\n');print(plane_results)
