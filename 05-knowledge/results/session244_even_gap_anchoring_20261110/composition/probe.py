from itertools import product
O=((1,0,0,0),(0,0,0,0),(0,0,2,0),(0,0,0,1))
Z=((1,-1,1,0),(0,0,0,1),(0,0,2,-1),(0,0,0,1))
J=((1,-1,1,0),(1,-1,0,0),(1,0,1,0),(1,0,0,0))
M=dict(O=O,Z=Z,J=J)
def mul(a,b):return tuple(tuple(sum(a[i][k]*b[k][j] for k in range(4)) for j in range(4)) for i in range(4))
I=tuple(tuple(int(i==j) for j in range(4)) for i in range(4))
states={'':I}
for n in range(1,10):
 states={w+c:mul(M[c],a) for w,a in states.items() for c in M}
 rows={a[0]:w for w,a in states.items() if w[-1]=='Z'}
 found=[]
 for r,w in rows.items():
  s=(r[0],r[1]+2,r[2],r[3])
  if s in rows:found.append((r,w,rows[s]))
 print(n,len(states),len(rows),len(found), sorted(found,key=lambda x:(sum(abs(a) for a in x[0]),x))[:8],flush=True)
