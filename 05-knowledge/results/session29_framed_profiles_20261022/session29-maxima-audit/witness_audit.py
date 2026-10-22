import itertools,json,pathlib,sys
src=pathlib.Path(sys.argv[1]) if len(sys.argv)>1 else pathlib.Path(__file__).resolve().parent.parent/'session29-halfpath-asymptotics'
results=[]
for w in json.loads((src/'witnesses.json').read_text()):
 m=w['m']; T=w['outmasks']; n=2*m
 assert all(((T[i]>>j)&1)+((T[j]>>i)&1)==1 for i in range(n) for j in range(i))
 assert all(((T[i]>>j)&1)==((T[j^1]>>(i^1))&1) for i in range(n) for j in range(n) if i!=j)
 assert all((T[2*i]>>(2*i+1))&1 for i in range(m))
 D=[[0]*n for _ in range(1<<m)]
 for v in range(n):D[1<<(v//2)][v]=1
 for mask in range(1,1<<m):
  for v in range(n):
   for z in range(n):
    if not(mask>>(z//2)&1) and (T[v]>>z&1):D[mask|1<<(z//2)][z]+=D[mask][v]
 f=sum(D[-1][::2]);assert f==w['fixed_hp']
 autos=[]
 for p in itertools.permutations(range(m)):
  a=[2*p[v//2]+v%2 for v in range(n)]
  if all(((T[i]>>j)&1)==((T[a[i]]>>a[j])&1) for i in range(n) for j in range(i)):
   autos.append(list(p))
 # Any automorphism commuting with sigma preserves its pair orbits. Internal
 # plus->minus arcs force it to preserve sign, so these exhaust the centralizer.
 assert f%len(autos)==0
 results.append(dict(m=m,fixed_paths=f,even_centralizer_order=len(autos),even_blue_tiles=f//len(autos),centralizer_pair_permutations=autos))
print(json.dumps(results,indent=2))
