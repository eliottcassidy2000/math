from itertools import permutations
import json,pathlib
P=pathlib.Path(__file__).resolve().parent

def graph(n,a=None):
 out=[0]*n
 for i in range(n):
  for j in range(i+1,n):
   x,y=(i,j) if j==i+1 else (j,i)
   if a is not None and (i,j)==(a,a+2):x,y=i,j
   out[x]|=1<<y
 return out

def cycles(out):
 n=len(out);found=set()
 def walk(start,path,mask):
  v=path[-1]
  if len(path)>=3 and out[v]>>start&1:found.add(tuple(path))
  for w in range(start+1,n):
   if not(mask>>w&1) and out[v]>>w&1:walk(start,path+[w],mask|1<<w)
 for v in range(n):walk(v,[v],1<<v)
 return found

def hp(out):
 n=len(out);dp=[[0]*n for _ in range(1<<n)]
 for v in range(n):dp[1<<v][v]=1
 for mask in range(1,1<<n):
  for v,c in enumerate(dp[mask]):
   if not c:continue
   z=out[v]&~mask
   while z:
    one=z&-z;z-=one;w=one.bit_length()-1;dp[mask|one][w]+=c
 return sum(dp[-1])
rows=json.loads((P/'RESULT.json').read_text())['rows'];count=0
for n,a,b,h in rows:
 if n>10:continue
 out=graph(n,a);assert hp(out)==h
 if n<=9:
  old=cycles(graph(n));new=cycles(out)
  add=set()
  for i in range(a+1):
   for j in range(a+2,n):
    if j-i>=3:add.add(tuple(range(i,a+1))+tuple(range(a+2,j+1)))
  for i in range(a):
   for j in range(a+3,n):add.add(tuple(range(i,a+1))+tuple(range(a+2,j+1))+(a+1,))
  assert old-new=={(a,a+1,a+2)} and new-old==add
  count+=1
print('PASS independent forward path DP through n10 and DFS directed-cycle inventories through n9:',count,'cycle cases')
