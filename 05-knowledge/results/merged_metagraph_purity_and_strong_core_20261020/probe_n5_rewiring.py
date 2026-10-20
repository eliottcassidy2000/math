from itertools import permutations,combinations
from collections import Counter
n=5
pairs=list(combinations(range(n),2)); tiles=[(x,y) for x in range(n) for y in range(x-1)]
perms=list(permutations(range(n))); full=(1<<len(tiles))-1

def matrix(t):
 A=[[0]*n for _ in range(n)]
 for x in range(1,n): A[x][x-1]=1
 for k,(x,y) in enumerate(tiles):
  if t>>k&1: A[x][y]=1
  else: A[y][x]=1
 return A

def canonical(A):
 return min(sum(A[p[i]][p[j]]<<k for k,(i,j) in enumerate(pairs)) for p in perms)

def key(t):
 A=matrix(t); return min(canonical(A),canonical([list(x) for x in zip(*A)]))

def refl(t): return sum(((t>>tiles.index((n-1-y,n-1-x)))&1)<<k for k,(x,y) in enumerate(tiles))
keys={t:key(t) for t in range(full+1)}
rep={v:min(t for t in keys if keys[t]==v) for v in set(keys.values())}
meta={}
for v,t in rep.items():
 A=matrix(t);s=tuple(sorted(map(sum,A)));c3=10-sum(x*(x-1)//2 for x in s)
 N=sum(k==v for k in keys.values());g=sum(keys[u]==v and refl(u)==u for u in keys)
 meta[v]=(t,s,c3,N,g)
lines=[]
for t in range(full+1):
 if t<(t^full):
  blue=refl(t)==t;a,b=keys[t],keys[t^full]
  lines.append((t,t^full,a,b,blue))
print('Coordinates:',tiles)
print('Merged node: representative, sorted scores, c3, fibre N, blue g')
for v,row in sorted(meta.items()):print(v,row)
print('Blue lines: mask pair, merged node codes, C3 pair')
for t,u,a,b,blue in lines:
 if blue:print(t,u,a,b,(meta[a][2],meta[b][2]))
for L1,L2 in combinations([l for l in lines if l[-1]],2):
 t,u,a,b,_=L1;v,w,c,d,_=L2
 if len({a,b,c,d})<4: continue
 if meta[a][2]%2==meta[c][2]%2:
  old=((a,b),(c,d));new=((a,c),(b,d))
 else:old=((a,b),(c,d));new=((a,d),(b,c))
 assert Counter(sum(old,()))==Counter(sum(new,()))
 assert all((meta[x][2]-meta[y][2])%2==0 for x,y in new)
 print('HOSTILE',{'old_lines':(L1,L2),'old_pairs':old,'new_pairs':new,'old_c3':[[meta[x][2],meta[y][2]] for x,y in old],'new_c3':[[meta[x][2],meta[y][2]] for x,y in new]})
 break
