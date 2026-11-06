from itertools import permutations
import json

def matrix(c,n):
 a=[[False]*n for _ in range(n)];k=0
 for i in range(n):
  for j in range(i+1,n):
   a[i][j]=bool(c>>k&1);a[j][i]=not a[i][j];k+=1
 return a

def words(a,vertices):
 return [p for p in permutations(vertices) if all(a[u][v] for u,v in zip(p,p[1:]))]
results={}
for c in (2,5):
 a=matrix(c,5);assert len(words(a,range(5)))==3
 out=set()
 for p in range(32):
  b=[r+[bool(p>>i&1)] for i,r in enumerate(a)]+[[not bool(p>>i&1) for i in range(5)]+[False]]
  h=len(words(b,range(6)));out.add(h)
  calc=0
  for s in range(32):
   aa=[i for i in range(5) if s>>i&1];bb=[i for i in range(5) if not(s>>i&1)]
   left=sum(not w or b[w[-1]][5] for w in words(a,aa))
   right=sum(not w or b[5][w[0]] for w in words(a,bb))
   calc+=left*right
  assert calc==h
 results[c]=sorted(out)
print(json.dumps(results,sort_keys=True))
assert results[2]==[3,5,9,13,15,17,23,25,27,33]
assert results[5]==[3,5,9,13,15,27,29]
print('PASS independent direct permutation enumeration, all 64 extensions, and cut formula')
