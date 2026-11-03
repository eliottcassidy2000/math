"""Exact all-word cone certificate; source atom checks by Hamilton subset DP."""
from collections import deque
from itertools import product
M={'a':((1,0,1),(2,0,1),(0,1,0)), 'b':((1,2,2),(4,3,4),(4,2,3)), 'q2':((5,3,4),(8,7,8),(6,6,7)), 'r21':((13,8,10),(20,14,17),(16,11,14)), 'r24':((5,16,16),(18,28,30),(18,22,24)), 'q4':((3,6,6),(6,8,9),(6,5,6)), 'r19':((11,7,10),(28,12,19),(26,11,16)), 'r32':((3,12,12),(12,22,23),(12,19,20)), 'r48':((5,16,16),(10,20,21),(10,13,14)), 'r18':((13,12,14),(28,19,24),(24,16,19)), 'r5':((9,5,6),(22,16,19),(20,15,18)), 'q1':((3,3,3),(12,6,9),(12,5,8)), 'q0':((1,4,4),(8,9,10),(8,8,9))}
G={'a':(0,), 'b':(2,4,1),'q2':(10,4,9,2),'r21':(22,20,8,19,4),'r24':(2,28,9,17,5),'q4':(2,12,9,1),'r19':(14,20,8,18,5),'r32':(2,4,25,19,3),'r48':(2,20,25,19,1),'r18':(10,20,9,18,5),'r5':(22,4,8,19,6),'q1':(6,4,8,3),'q0':(2,4,9,3)}
W={'L2':('r24','a','a','a'),'L6':('r24','q4'),'L45':('r19','r32','r48','q4'),'L13':('r19','r19','a','a','r18'),'b':('b',),'q2':('q2',),'Lp2':('r24','a','a','r21')}
I3=((1,0,0),(0,1,0),(0,0,1)); B=((0,1),(0,2),(1,-1)); t=(2,-1,0)
def check(ok,msg):
 if not ok: raise RuntimeError(msg)
def mm(A,B):return tuple(tuple(sum(x*y for x,y in zip(row,col)) for col in zip(*B)) for row in A)
def vm(v,A):return mm((v,),A)[0]
def mat(word):
 p=I3
 for x in word:p=mm(p,M[x])
 return p
def graph(word):
 rows=[]
 for l in reversed(word):
  src=G[l]; n=len(rows); m=len(src)
  rows=[r|(((1<<m)-1)<<n) for r in rows]+[r<<n for r in src]
 n=len(rows)
 for i in range(n):
  for j in range(i+2,n):rows[i]^=1<<j;rows[j]^=1<<i
 return tuple(sum(1<<(n-1-j) for j in range(n) if rows[n-1-i]>>j&1) for i in range(n))
def hp(rows):
 n=len(rows); dp=[[0]*n for _ in range(1<<n)]
 for i in range(n):dp[1<<i][i]=1
 for mask in range(1,1<<n):
  for i in range(n):
   if dp[mask][i]:
    for j in range(n):
     if not mask>>j&1 and rows[i]>>j&1:dp[mask|1<<j][j]+=dp[mask][i]
 return sum(dp[-1])
def add_a(rows):
 n=len(rows);return tuple(r | ((1<<n) if i<n-1 else 0) for i,r in enumerate(rows))+(1<<(n-1),)
def observed(word):
 rows=graph(word);h=hp(rows);ha=hp(add_a(rows));haa=hp(add_a(add_a(rows)))
 check((ha-h)%2==0 and (haa-h)%2==0,'observer parity')
 x=(ha-h)//2;y=(haa-h)//2-x
 return (h,x,y)
for k,P in M.items():
 v0=observed((k,));v1=observed(('a',k));v2=observed(('a','a',k))
 rows=(v0,tuple(a-b for a,b in zip(v1,v0)),tuple(a-b for a,b in zip(v2,v1)))
 # rows above are row0,row2,row1
 check((rows[0],rows[2],rows[1])==P,'independent atom '+k)
for u,expected in [(('r5','a','a'),(31,14,25)),(('a','a','q1','a','a'),(95,41,75))]:
 check(observed(u)==expected,'source observer')
print('13 atom matrices and two actual seed observers checked by literal graph DP')
U=(31,14,25); U2=(95,41,75)
check(tuple(U2[i]-3*U[i] for i in range(3))==t,'seed relation')
macro={k:mat(v) for k,v in W.items()}; quot={}
for k,p in macro.items():
 z=mm(p,B); q=((z[0][0]+z[2][0],z[0][1]+z[2][1]),z[0])
 check(mm(B,q)==z,'quotient '+k)
 gt=vm(t,p);check(gt[2]==0 and gt[0]==-2*gt[1], 'return '+k)
 quot[k]=q
 if k not in ('b',):
  low=all(11*y<17*x for x,y in q)
  high=all(5*y>8*x for x,y in q)
  check(low or high,'slope '+k)
  print(k,'low' if low else 'high','gain',-gt[1])
check(quot['b']==((5,8),(2,3)),'b quotient')
D=mm(mm(macro['L2'],macro['b']),macro['L2']);QD=mm(mm(quot['L2'],quot['b']),quot['L2'])
check(vm(t,D)==tuple(-4*z for z in t),'D gain')
check(QD[0][0]*QD[1][1]-QD[0][1]*QD[1][0]!=0,'D invertible')
check(11*34<17*25,'seed ratio low')
# An explicit nonempty mixed prefix: A'=Lp2 L6, A=A' b.
P=mm(macro['Lp2'],macro['L6']);v=vm(vm(vm(U,P),macro['b']),D);w=vm(U2,P)
lhs=vm(vm(w,macro['b']),D);g=-vm(t,P)[1];diff=tuple(lhs[i]-3*v[i] for i in range(3))
check(diff==tuple(4*g*z for z in t),'mixed physical identity')
print('mixed witness A prime = Lp2 L6, carry =',4*g)
print('w =',w,'v =',v)
print('PASS')
