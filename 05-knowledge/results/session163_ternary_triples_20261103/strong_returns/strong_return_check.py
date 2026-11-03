"""Literal source-state controls for the two offset-then-zero branches."""

def need(p,msg):
 if not p: raise AssertionError(msg)

def H(rows):
 n=len(rows); D=[[0]*n for _ in range(1<<n)]
 for i in range(n): D[1<<i][i]=1
 for mask in range(1,1<<n):
  for i in range(n):
   a=D[mask][i]
   if a:
    z=rows[i]&~mask
    while z:
     b=z&-z;j=b.bit_length()-1;D[mask|b][j]+=a;z-=b
 return sum(D[-1])

def a(rows):
 n=len(rows)
 return tuple(r|((1<<n) if i<n-1 else 0) for i,r in enumerate(rows))+(1<<(n-1),)

def state(rows):
 h=H(rows);ha=H(a(rows));haa=H(a(a(rows)))
 need((ha-h)%2==0 and (haa-ha)%2==0,'parity')
 return h,(ha-h)//2,(haa-ha)//2

def strong(rows):
 n=len(rows)
 for root in range(n):
  todo=[root];seen={root}
  while todo:
   u=todo.pop()
   for v in range(n):
    if rows[u]>>v&1 and v not in seen:seen.add(v);todo.append(v)
  if len(seen)!=n:return False
 return True

def join(rows,block):
 n=len(rows); r=list(rows); allnew=(1<<len(block))-1
 for i in range(n):r[i] |= ((allnew ^ (1 if i==n-1 else 0))<<n)
 return tuple(r)+tuple((x<<n) | ((1<<(n-1)) if i==0 else 0) for i,x in enumerate(block))

def mul(v,P):return tuple(sum(v[i]*P[i][j] for i in range(3)) for j in range(3))
P=((5,3,4),(8,7,8),(6,6,7));B=((1,2,2),(4,3,4),(4,2,3));t=(2,-1,0)
need(mul(t,P)==t and mul(t,B)==tuple(-v for v in t),'t gains')
seed=(112,97,115,103,74,16,32)
out1=(252,249,194,228,204,148,32,64)
out2=(496,497,499,487,456,400,288,64,128)
ret1=(1008,1017,1011,997,776,912,816,592,128,256)
ret2=(996,993,994,999,975,912,800,576,128,256)
expect=((31,14,25),(95,41,75),(97,40,75),(285,123,225),(291,120,225))
for rows,v in zip((seed,out1,out2,ret1,ret2),expect):
 need(all(not (rows[i]>>i&1) and all(((rows[i]>>j)&1)+((rows[j]>>i)&1)==1 for j in range(len(rows)) if i!=j) for i in range(len(rows))), 'tournament orientation')
 need(strong(rows),'source strong');need(state(rows)==v,'literal state')
 # Test q2 transfer at full source state on actual tournament via 3 H DPs.
 extended=join(rows,(8,1,3,6));need(strong(extended),'q2 extension strong')
 need(state(extended)==mul(v,P),'q2 literal state transfer')
 extended=join(rows,(0,1,3));need(strong(extended),'b extension strong')
 need(state(extended)==mul(v,B),'b literal state transfer')
 print('n',len(rows),'state',v,'q2 state',mul(v,P))
U,U2,U4,R2,R4=expect
need(U2==tuple(3*U[i]+t[i] for i in range(3)),'offset 2')
need(U4==tuple(3*U[i]+2*t[i] for i in range(3)),'offset 4')
need(R2==tuple(3*x for x in U2) and R4==tuple(3*x for x in U4),'both return')
print('PASS: five strong literal states, b/q2 transfers, and two offset-then-zero returns')
