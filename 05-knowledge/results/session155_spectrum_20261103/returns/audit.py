from math import gcd
from itertools import permutations

P=((336359,155087,268050),(672726,310170,536100),(560922,258618,447000))
B=((0,1),(0,2),(1,-1))
Q=((715050,1009641),(268050,378483))
R=((238350,336547),(89350,126161))
t=(2,-1,0)
seeds=((31,14,25),(95,41,75),(97,40,75))
rows=((112,97,115,103,74,16,32),(252,249,194,228,204,148,32,64),(496,497,499,487,456,400,288,64,128))
def mm(A,B):return tuple(tuple(sum(x*y for x,y in zip(row,col)) for col in zip(*B)) for row in A)
def mv(v,A):return mm((v,),A)[0]
def eye(n):return tuple(tuple(int(i==j) for j in range(n)) for i in range(n))
def pw(A,n):
 out=eye(len(A))
 for _ in range(n):out=mm(out,A)
 return out
def beta(v):return mv(v,B)
def hp(rows):
 n=len(rows); dp=[[0]*n for _ in range(1<<n)]
 for i in range(n):dp[1<<i][i]=1
 for mask in range(1,1<<n):
  for i in range(n):
   if not dp[mask][i]:continue
   for j in range(n):
    if not (mask>>j)&1 and (rows[i]>>j)&1:dp[mask|1<<j][j]+=dp[mask][i]
 return sum(dp[-1])
if not (mm(P,B)==mm(B,Q)): raise AssertionError('mm(P,B)==mm(B,Q)')
if not (all(Q[i][j]==3*R[i][j] for i in range(2) for j in range(2))): raise AssertionError('all(Q[i][j]==3*R[i][j] for i in range(2) for j in range(2))')
if not (mv(t,P)==tuple(-4*x for x in t)): raise AssertionError('mv(t,P)==tuple(-4*x for x in t)')
if not (R[0][0]*R[1][1]-R[0][1]*R[1][0]==-100): raise AssertionError('R[0][0]*R[1][1]-R[0][1]*R[1][0]==-100')
if not (all((a+50*b)%100==0 for a,b in R)): raise AssertionError('all((a+50*b)%100==0 for a,b in R)')
if not ([beta(v) for v in seeds]==[(25,34),(75,102),(75,102)]): raise AssertionError('[beta(v) for v in seeds]==[(25,34),(75,102),(75,102)]')
if not ([hp(x) for x in rows]==[31,95,97]): raise AssertionError('[hp(x) for x in rows]==[31,95,97]')
print('seed totals / bulk:', [(hp(x),beta(v)) for x,v in zip(rows,seeds)])
for k in range(1,9):
 M=pw(R,k); A,Bb=M[0]; C,D=M[1]; N=100**k
 if not (D%10==1 and A*D-Bb*C==(-100)**k): raise AssertionError('D%10==1 and A*D-Bb*C==(-100)**k')
 for a in range(-20,21):
  for b in range(-20,21):
   x=a*A+b*C; y=a*Bb+b*D
   if not ((x*D-y*C)%N==0): raise AssertionError('(x*D-y*C)%N==0')
 print('power',k,'index',N,'second-row-second-entry mod10',D%10)
for v in seeds:
 for j in range(5):
  for k in range(1,j+1):
   c,d=mv(beta(v),pw(R,j-k));c*=3**j;d*=3**j
   w=(c+d,0,c) # beta(w)=(c,d)
   lhs=mv(w,pw(P,k));rhs=tuple(3**k*x for x in mv(v,pw(P,j)))
   diff=tuple(x-y for x,y in zip(lhs,rhs));m=-diff[1]
   if not (diff==tuple(m*x for x in t)): raise AssertionError('diff==tuple(m*x for x in t)')
print('formal integer witnesses, all seeds, 0<=j<=4 and 1<=k<=j: PASS')
print('all target bulks odd first coordinate excluded at k=1: PASS')
def add_a(rows):
 n=len(rows);return tuple(r | ((1<<n) if i<n-1 else 0) for i,r in enumerate(rows)) + (1<<(n-1),)
for v,g in zip(seeds,rows):
 h,x,y=v
 if not (hp(add_a(g))==h+2*x): raise AssertionError('hp(add_a(g))==h+2*x')
 if not (hp(add_a(add_a(g)))==h+2*(x+y)): raise AssertionError('hp(add_a(add_a(g)))==h+2*(x+y)')
 print('direct graph h,h_a,h_aa:',h,hp(add_a(g)),hp(add_a(add_a(g))))
M={'a':((1,0,1),(2,0,1),(0,1,0)),'b':((1,2,2),(4,3,4),(4,2,3)), 'r24':((5,16,16),(18,28,30),(18,22,24)),'r5':((9,5,6),(22,16,19),(20,15,18)), 'q1':((3,3,3),(12,6,9),(12,5,8)), 'q0':((1,4,4),(8,9,10),(8,8,9))}
G={'a':(0,), 'b':(2,4,1),'r24':(2,28,9,17,5),'r5':(22,4,8,19,6),'q1':(6,4,8,3),'q0':(2,4,9,3)}
def wordmat(word):
 m=eye(3)
 for s in word:m=mm(m,M[s])
 return m
def graph(word):
 rows=[]
 for l in reversed(word):
  src=G[l];n=len(rows);m=len(src)
  rows=[r|(((1<<m)-1)<<n) for r in rows]+[r<<n for r in src]
 n=len(rows)
 for i in range(n):
  for j in range(i+2,n): rows[i]^=1<<j;rows[j]^=1<<i
 return tuple(sum(1<<(n-1-j) for j in range(n) if rows[n-1-i]>>j&1) for i in range(n))
words=(('r5','a','a'),('a','a','q1','a','a'),('q0','a','a','a','a','a'))
D=('r24','a','a','a','b','r24','a','a','a')
if not (wordmat(D)==P): raise AssertionError('wordmat(D)==P')
for word,v,expected in zip(words,seeds,rows):
 if not (mv((1,0,0),wordmat(word))==v and graph(word)==expected): raise AssertionError('mv((1,0,0),wordmat(word))==v and graph(word)==expected')
print('independent atom matrix product, seed reconstruction: PASS')
# General theorem: construct targets whose bulk is in L_(k-j), for k>j.
for j in range(4):
 for k in range(j+1,6):
  for q in ((1,0),(0,1),(3,-2)):
   a,b=mv(q,pw(R,k-j));v=(a+b,0,a)
   a,b=tuple(3**j*z for z in q);w=(a+b,0,a)
   left=mv(w,pw(P,k));right=tuple(3**k*z for z in mv(v,pw(P,j)))
   diff=tuple(x-y for x,y in zip(left,right));m=-diff[1]
   if not (diff==tuple(m*z for z in t)): raise AssertionError('general criterion constructive test')
print('general criterion constructive tests for k>j: PASS')
