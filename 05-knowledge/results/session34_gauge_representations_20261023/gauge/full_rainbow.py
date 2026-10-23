from itertools import combinations, permutations
from collections import Counter
import json
E=list(combinations(range(7),2));idx={e:i for i,e in enumerate(E)}
L=[{i,(i+1)%7,(i+3)%7} for i in range(7)]
el={e:next(i for i,l in enumerate(L) if set(e)<=l) for e in E}
cycles=[]
for vs in combinations(range(7),4):
 for p in permutations(vs[1:]):
  c=(vs[0],)+p
  if c[1]>c[-1]:continue
  edges=[tuple(sorted((c[i],c[(i+1)%4]))) for i in range(4)]
  if len({el[e] for e in edges})==4:cycles.append(sum(1<<idx[e] for e in edges))
basis={}
for c in cycles:
 x=c
 while x:
  j=x.bit_length()-1
  if j in basis:x^=basis[j]
  else:basis[j]=x;break
free=[i for i,e in enumerate(E) if 0 not in e]
hist=Counter()
for m in range(1<<len(free)):
 x=sum(((m>>j)&1)<<k for j,k in enumerate(free))
 R=sum((-1)**((x&c).bit_count()%2) for c in cycles)
 hist[R]+=64

# Certificate: 11 independent original rainbow rows, one XOR echelon basis.
original=[];red={}
for c in cycles:
 x=c
 while x:
  j=x.bit_length()-1
  if j in red:x^=red[j]
  else:red[j]=x;original.append(c);break
if len(basis)!=11 or len(original)!=11 or sum(hist.values())!=2**21 or hist[21]!=1024:
 raise RuntimeError('rank or census failed')
A=[4,1,3,7,6,5,2]
line=[sum(1<<idx[e] for e in E if el[e]==i) for i in range(7)]
code=[]
for m in range(128):
 syndrome=0;em=0
 for i in range(7):
  if m>>i&1:syndrome^=A[i];em^=line[i]
 if syndrome==0:code.append(em)
for em in code:
 for c in cycles:
  if (em&c).bit_count()%2:raise RuntimeError('code does not annihilate')
if len(code)!=16:raise RuntimeError('code dimension failed')
print(json.dumps({'rainbow_rank':len(basis),'full_field_count':sum(hist.values()),'R_histogram':dict(sorted(hist.items())),'independent_original_rows':original,'echelon_rows':sorted(red.items(),reverse=True),'kernel_contains_code_words':len(code)},indent=2))
