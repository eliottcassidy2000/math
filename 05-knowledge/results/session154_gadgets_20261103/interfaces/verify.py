from itertools import permutations,product

def original(code=41):
 a=[[0]*5 for _ in range(5)];k=0
 for i in range(5):
  for j in range(i+1,5):
   if (code>>k)&1:a[i][j]=1
   else:a[j][i]=1
   k+=1
 return a

def literal(a):
 return sum(all(a[p[i]][p[i+1]] for i in range(len(p)-1)) for p in permutations(range(len(a))))

def strong(a):
 n=len(a);r=[row[:] for row in a]
 for k in range(n):
  for i in range(n):
   for j in range(n):r[i][j]|=r[i][k] and r[k][j]
 return all(i==j or r[i][j] for i in range(n) for j in range(n))

def build(v,masks):
 a=original(); old=[u for u in range(5) if u!=v];b=[[0]*7 for _ in range(7)]
 for i,j in product(range(4),repeat=2):b[i][j]=a[old[i]][old[j]]
 for t,m in enumerate(masks):
  for j in range(4):
   if (m>>j)&1:b[4+t][j]=1
   else:b[j][4+t]=1
 b[4][5]=b[5][6]=b[6][4]=1
 return b

def covers(a,v):
 U=[u for u in range(len(a)) if u!=v]; ret=[0,0,0,0]
 # enumerate a permutation and all weak cuts into r+1 directed paths.
 for r in (2,3):
  for p in permutations(U):
   for cuts in product(range(len(U)+1),repeat=r):
    if tuple(sorted(cuts))!=cuts:continue
    q=(0,)+cuts+(len(U),); seg=[p[q[i]:q[i+1]] for i in range(r+1)]
    if any(not seg[i] for i in range(1,r)):continue
    if any(any(not a[s[j]][s[j+1]] for j in range(len(s)-1)) for s in seg):continue
    if seg[0] and not a[seg[0][-1]][v]:continue
    if seg[-1] and not a[v][seg[-1][0]]:continue
    if any(not a[v][seg[i][0]] or not a[seg[i][-1]][v] for i in range(1,r)):continue
    ret[r]+=1
 return ret

a=original(); print('seed H',literal(a),'strong',strong(a))
for masks in [(6,1,1),(3,4,1),(15,7,1),(15,15,9)]:
 b=build(0,masks); print('ports',masks,'H',literal(b),'strong',strong(b))
for v in range(5):
 old=[u for u in range(5) if u!=v];m=sum(a[v][u]<<i for i,u in enumerate(old));b=build(v,(m,m,m)); p=covers(a,v); x=literal(b);expected=3*literal(a)+6*(p[2]+p[3]);
 if x!=expected:raise RuntimeError((v,p,x,expected))
 print('uniform root',v,'P2',p[2],'P3',p[3],'H',x)
