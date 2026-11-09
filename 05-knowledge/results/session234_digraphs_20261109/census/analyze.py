import json,itertools,collections
# Avoid importing scan.py which runs enumeration.
r=json.load(open('results.json'))
def direct(a):return sum(all(a[p[i]][p[i+1]] for i in range(len(a)-1)) for p in itertools.permutations(range(len(a))))
for n in ('4','5'):
 for h,cs in r[n]['target_classes'].items():
  for c in cs.values():
   a=c['matrix'];d=[(i,j) for i in range(int(n)) for j in range(i+1,int(n)) if a[i][j] and a[j][i]]
   spectrum=collections.Counter()
   for bits in itertools.product((0,1),repeat=len(d)):
    b=[row[:] for row in a]
    for (i,j),bit in zip(d,bits):b[i][j]=bit;b[j][i]=1-bit
    spectrum[direct(b)]+=1
   c['reciprocal_pairs']=d;c['tournament_resolution_histogram']=dict(sorted(spectrum.items()))
   c['hamiltonian_paths']=[list(p) for p in itertools.permutations(range(int(n))) if all(a[p[i]][p[i+1]] for i in range(int(n)-1))]
json.dump(r,open('results.json','w'),indent=2)
for n in ('4','5'):
 for h,cs in r[n]['target_classes'].items():
  print(n,h,'reciprocal counts',dict(collections.Counter(len(c['reciprocal_pairs']) for c in cs.values())))
  for c in cs.values():print(c['code'],c['tournament_resolution_histogram'])
