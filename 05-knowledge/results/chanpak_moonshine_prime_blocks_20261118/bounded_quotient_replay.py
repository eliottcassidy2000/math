from itertools import product, combinations
def H(A):
 n=len(A); D=[[0]*n for _ in range(1<<n)]
 for i in range(n): D[1<<i][i]=1
 for S,row in enumerate(D):
  for i,k in enumerate(row):
   if not k: continue
   V=A[i]&~S
   while V:
    b=V&-V;V^=b;D[S|b][b.bit_length()-1]+=k
 return sum(D[-1])
def inflate(Q,m):
 lab=[i for i,k in enumerate(m) for _ in range(k)]
 return [sum(1<<j for j,b in enumerate(lab) if (i<j if a==b else (Q[a]>>b)&1)) for i,a in enumerate(lab)]
Q=[2,4,1]
for a,b in product(range(1,5),repeat=2):
 z=H(inflate(Q,(a,b,1)))
 if z!=((1<<a)-1)*((1<<b)-1)+2: raise ValueError((a,b,z))
tests=0
E=list(combinations(range(4),2))
for code in range(64):
 Q=[0]*4
 for k,(u,v) in enumerate(E):
  if not code>>k&1:u,v=v,u
  Q[u]|=1<<v
 active={i for i in range(4) if any(Q[i]>>j&1 and Q[j]>>k&1 and Q[k]>>i&1 for j,k in product(range(4), repeat=2))}
 for m in ((1,2,3,2),(3,1,2,2),(2,3,1,2),(2,2,3,1)):
  z=H(inflate(Q,m)); normal=tuple(a if i in active else 1 for i,a in enumerate(m))
  if z!=H(inflate(Q,normal)) or any(z<(1<<m[i])+1 for i in active):raise ValueError((Q,m,z))
  tests+=1
print("PASS two-knob identity, inactive normalization and active lower bound:",tests,"four-vertex controls")
