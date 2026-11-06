import itertools,json,math
from pathlib import Path
R=Path(__file__).resolve().parent
selected={9:{1875,2025,2187},10:{10125,9375},11:{54675,46875},12:{151875,253125}}
out=[]
for n,need in selected.items():
 for line in (R/f'search{n}.out').read_text().splitlines():
  nn,h,a,b,code,t=map(int,line.split())
  if h not in need:continue
  assert nn==n and h==3**a*5**b
  incoming=[0]*n;reach=[1<<i for i in range(n)]
  for k,(i,j) in enumerate(itertools.combinations(range(n),2)):
   u,v=(i,j) if code>>k&1 else(j,i);incoming[v]|=1<<u;reach[u]|=1<<v
  for k in range(n):
   for i in range(n):
    if reach[i]>>k&1:reach[i]|=reach[k]
  assert all(x==(1<<n)-1 for x in reach)
  f=[[0]*n for _ in range(1<<n)]
  for i in range(n):f[1<<i][i]=1
  for mask in range(1,1<<n):
   if mask&(mask-1)==0:continue
   for v in range(n):
    if mask>>v&1:f[mask][v]=sum(f[mask^(1<<v)][u] for u in range(n) if (mask&incoming[v])>>u&1)
  assert sum(f[-1])==h
  out.append(dict(n=n,h=h,exponents=[a,b],code=str(code),iteration=t,strong=True,verified_count=sum(f[-1])))
assert len(out)==sum(map(len,selected.values()))
assert 1875**2*2025**3==15**14
assert 10125**4*9375==15**17
assert 54675*46875==15**8
assert 151875*253125==15**9
assert 16+24*4>11*9 # Old dual fails at order nine.
# Elementary bound on all permutation words, checked on witness controls.
for row in out:assert row['h']<=math.factorial(row['n'])//2**(row['n']//2)
(R/'WITNESSES.json').write_text(json.dumps(out,indent=2)+'\n')
result={'witnesses_checked':len(out),'global_upper_rate':'8/3','block_exponent':9,'block_order':24,'blocks':[[14,45],[17,50],[8,22],[9,24]],'old_dual_extends_globally':False,'search_exhaustive':False,'checks':'pass'}
(R/'RESULT.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
