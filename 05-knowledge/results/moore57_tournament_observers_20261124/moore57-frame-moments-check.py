from itertools import combinations

A8=[19,8,2,25,17,10,8,2,5,13,26,16,13,26,19,15,21,8,9,3,10,2,18,19,14,4,7,5]
A9=[14,9,9,8,0,6,7,9,0,1,8,27,10,18,24,14,9,25,26,25,11,0,22,18,14,17,8,7,22,11,25,18,5,25,18,18]

def check(k, offsets, components):
 assert len(offsets)==k*(k-1)//2 and sum(components)==56
 a=[[0]*k for _ in range(k)]
 for (i,j),v in zip(combinations(range(k),2),offsets):
  a[i][j]=v
  a[j][i]=-v
 for s in components:
  for i,j in combinations(range(k),2):
   # Direct arc and every possible two-matching word have distinct gains.
   routes=[a[i][j]]+[a[i][t]+a[t][j] for t in range(k) if t not in (i,j)]
   assert len({v%s for v in routes})==k-1
 return k-1

assert check(8,A8,[27,29])==7 and (8*27)%56!=0
assert check(9,A9,[28,28])==8
print('two explicit locally triangle/square-free but globally rejected block systems')
