import contextlib,io,itertools,json
with contextlib.redirect_stdout(io.StringIO()):
 from independent import subset_count
 from sandwich import fuse

def context(bits):
 n=len(bits)+1;G=[[False]*n for _ in range(n)]
 for i,b in enumerate(bits,1):
  G[0][i]=bool(b);G[i][0]=not b
  for j in range(i+1,n):G[i][j]=True
 return G
Gs=[[[bool(r>>j&1) for j in range(6)] for r in rows] for rows in [(56,13,25,48,34,6),(52,13,24,49,34,6),(52,13,24,33,42,6)]]
checks=0;shifts=set()
for n in range(6):
 for bits in itertools.product(range(2),repeat=n):
  C=context(bits);left=[subset_count(fuse(C,G)) for G in Gs]
  right=[subset_count(fuse(G,C)) for G in Gs]
  assert (left[2]-left[1]==2)==(sum(bits)<=1)
  assert (right[2]-right[1]==2)==all(bits)
  assert left[1]-left[0]==2 and right[1]-right[0]==2
  if sum(bits)<=1:
   t=(n-1-bits.index(1)) if sum(bits) else None
   shift=0 if t is None else 48 if t==0 else 57*2**t
   assert left==[41+shift,43+shift,45+shift];shifts.add(shift)
  checks+=1
print(json.dumps({'status':'PASS','all_unit_deletion_contexts_spine_lengths':'0..5','contexts':checks,'translated_shifts':sorted(shifts)},indent=2))
