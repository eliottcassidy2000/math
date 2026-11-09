import csv,json,math,itertools
from pathlib import Path
root=Path(__file__).parent
scan=root.parent/'census'
def hist(b,n):
 return {int(r['H']):int(r['count']) for r in csv.DictReader(open(root/f'hist-base{b}-n{n}.csv'))}
def classes(n):
 out={}
 for line in open(root/f'semicomplete-{n}.txt'):
  if line.startswith('class '):
   d=dict(t.split('=') for t in line.split()[1:]);out[int(d['bits'])]=(int(d['H']),int(d['strong']),int(d['reciprocal']))
 return out
raw=json.load(open(scan/'results.json'))
for b,orders in [(3,range(1,7)),(2,range(3,8))]:
 for n in orders:
  h=hist(b,n); m=n*(n-1)//2
  assert sum(h.values())==b**m
  expected=math.factorial(n)*(2**(n-1)*3**(m-n+1) if b==3 else 2**(m-n+1))
  assert sum(k*v for k,v in h.items())==expected
  if b==3:
   if n<6: ref={int(k):v for k,v in raw[str(n)]['histogram'].items()}
   else: ref={int(r['hamiltonian_paths']):int(r['labelled_count']) for r in csv.DictReader(open(scan/'results6.csv'))}
   assert h==ref
  print(f'PASS base={b} n={n}: graph total, path-incidence total'+(', complete scanner histogram equality' if b==3 else ''))
a={k:v for k,v in classes(6).items() if v[0] in [7,21]}
b={int(r['canonical_binary_integer']):(int(r['hamiltonian_paths']),int(r['strong']),int(r['reciprocal_pairs'])) for r in csv.DictReader(open(scan/'classes6.csv'))}
assert a==b
print('PASS all n6 target canonical representatives, strongness and reciprocal counts')
for n,H,expected in [(3,3,12),(4,7,48),(5,21,480),(6,21,1024)]:
 total=0;autstats=[]
 for bits,(h,s,r) in classes(n).items():
  if h!=H or not s:continue
  rows=[[int(c) for c in format(bits,f'0{n*n}b')[i*n:(i+1)*n]] for i in range(n)]
  autos=[p for p in itertools.permutations(range(n)) if all(rows[i][j]==rows[p[i]][p[j]] for i in range(n) for j in range(n))]
  masks=set()
  for mask in range(1<<n):
   masks.add(min(sum(((mask>>p[i])&1)<<i for i in range(n)) for p in autos))
  total+=len(masks);autstats.append((len(autos),len(masks)))
 assert total==expected
 print(f'PASS diagonal-decorated strong cores n={n} H={H}: {total}; (automorphism order, loop orbits)={autstats}')
print('ALL CHECKS PASS')
