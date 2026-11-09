import json,csv,itertools,collections

def matrix(n,code):
 a=[[0]*n for _ in range(n)]
 for i in range(n):
  for j in range(i+1,n):
   s=code%3;code//=3;a[i][j]=int(s!=1);a[j][i]=int(s!=0)
 return a

def paths(a):return [list(p) for p in itertools.permutations(range(len(a))) if all(a[p[i]][p[i+1]] for i in range(len(a)-1))]
records=[]
for row in csv.DictReader(open('classes6.csv')):
 r={k:int(v) for k,v in row.items()};a=matrix(6,r['representative_code']);r['matrix']=a;r['paths']=paths(a)
 assert len(r['paths'])==r['hamiltonian_paths']
 pairs=[(i,j) for i in range(6) for j in range(i+1,6) if a[i][j] and a[j][i]]
 r['reciprocal_pair_locations']=pairs;spec=collections.Counter();resolutions=[]
 for bits in itertools.product((0,1),repeat=len(pairs)):
  b=[row[:] for row in a]
  for (i,j),bit in zip(pairs,bits):b[i][j]=bit;b[j][i]=1-bit
  h=len(paths(b));spec[h]+=1;resolutions.append({'matrix':b,'hamiltonian_paths':h})
 r['tournament_resolution_histogram']=dict(sorted(spec.items()));r['tournament_resolutions']=resolutions
 records.append(r)
json.dump(records,open('classes6.json','w'),indent=2)
for r in records:
 if r['strong'] and r['reciprocal_pairs']==1:print(r['representative_code'],' / '.join(''.join(map(str,x)) for x in r['matrix']),r['tournament_resolution_histogram'],r['paths'])
