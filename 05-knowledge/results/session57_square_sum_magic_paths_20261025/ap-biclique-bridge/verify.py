import itertools,json
import sympy as s
V=tuple(range(6)); ed=lambda p:frozenset(tuple(sorted((u,v))) for u,v in zip(p,p[1:]))
E=frozenset((i,j) for i in range(3) for j in range(3,6))
paths=[p for p in itertools.permutations(V) if p<p[::-1] and ed(p)<=E]
pairs=[(p,q) for p,q in itertools.combinations(paths,2) if ed(p)|ed(q)==E]
cycles=set()
for p in itertools.permutations(V):
 if all(tuple(sorted((p[i],p[(i+1)%6]))) in E for i in range(6)):
  cycles.add(frozenset(tuple(sorted((p[i],p[(i+1)%6]))) for i in range(6)))
M=s.zeros(9,6)
for i in range(3):
 for j in range(3):M[3*i+j,i]=M[3*i+j,3+j]=1
R=[]
for i,j in itertools.product((1,2),repeat=2):
 r=[0]*9
 for a,b,v in [(i,j,1),(i,0,-1),(0,j,-1),(0,0,1)]:r[3*a+b]+=v
 R.append(r)
for cells in [[(0,0,1),(1,0,-2),(2,0,1)],[(0,0,1),(0,1,-2),(0,2,1)]]:
 r=[0]*9
 for i,j,v in cells:r[3*i+j]=v
 R.append(r)
a1,a2,a3,b1,b2,b3=s.symbols('a1 a2 a3 b1 b2 b3')
G=s.Matrix([[a3+b3,a2+b2,a1+b1],[a1+b2,a3+b1,a2+b3],[a2+b1,a1+b3,a3+b2]])
T=a1+a2+a3+b1+b2+b3
assert all(s.expand(sum(G.row(i))-T)==0 and s.expand(sum(G.col(i))-T)==0 for i in range(3))
D=[s.expand(sum(G[i,i] for i in range(3))-T),s.expand(sum(G[i,2-i] for i in range(3))-T)]
P=(0,3,1,4,2,5);Q=(4,0,5,1,3,2)
assert ed(P)|ed(Q)==E and len(ed(P)&ed(Q))==1
out={'hamilton_paths_unoriented':len(paths),'hamilton_cycles_unoriented':len(cycles),'unordered_two_path_covers':len(pairs),'cover_partners_per_path':sorted(set(sum(p in pq for pq in pairs) for p in paths)),'cover_intersections':sorted(set(len(ed(p)&ed(q)) for p,q in pairs)),'unsigned_incidence_rank':M.rank(),'cycle_constraints_rank':s.Matrix(R[:4]).rank(),'cycle_plus_AP_rank':s.Matrix(R).rank(),'diagonal_defects':[str(d) for d in D],'P':P,'Q':Q,'shared_edge':list(ed(P)&ed(Q)),'cover_pairs':pairs}
assert out['hamilton_paths_unoriented']==36 and out['hamilton_cycles_unoriented']==6
assert M.rank()==5 and s.Matrix(R).rank()==6
adj={p:set() for p in paths}
for p,q in pairs:adj[p].add(q);adj[q].add(p)
seen=set();sizes=[]
for p in paths:
 if p in seen:continue
 stack=[p];comp=set()
 while stack:
  z=stack.pop()
  if z in comp:continue
  comp.add(z);stack.extend(adj[z]-comp)
 seen|=comp;sizes.append(len(comp))
assert sorted(sizes)==[6]*6
out['cover_partner_graph_component_sizes']=sorted(sizes)
print(json.dumps(out,indent=2));open(__import__('pathlib').Path(__file__).with_name('results.json'),'w').write(json.dumps(out,indent=2))
