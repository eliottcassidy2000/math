from itertools import product,permutations
def check(p,x):
 if not p:raise ValueError(x)
def move(s,b):
 h,e,t=s
 return (h+t-e,1,2*t-1) if b else (h,0,2*t)
def add(T,b):
 n=len(T)
 return [r|((1<<n) if i or b==0 else 0) for i,r in enumerate(T)]+[b]
def W(word):
 A=[2,0]
 for b in word:A=add(A,int(b))
 return A
def h(A):
 n=len(A);D=[[0]*n for _ in range(1<<n)]
 for i in range(n):D[1<<i][i]=1
 for S,row in enumerate(D):
  for i,v in enumerate(row):
   if v:
    w=A[i]&~S
    while w:
     z=w&-w;w-=z;D[S|z][z.bit_length()-1]+=v
 return sum(D[-1])
def slot(A,r):
 n=len(A);C=[0]*n
 for P in permutations(i for i in range(n) if i!=r):
  for s in range(1,1<<n):
   if s&(s<<1):continue
   V=[]
   for i in range(n):
    if s>>i&1:V.append(r)
    if i<n-1:V.append(P[i])
   if all(A[i]>>j&1 for i,j in zip(V,V[1:])):C[s.bit_count()]+=1
 return C[1:]
def inflate(A,r,m):
 labs=sum(([i]*(m if i==r else 1) for i in range(len(A))),[])
 return [sum(1<<j for j,y in enumerate(labs) if (i<j if x==y else A[x]>>y&1)) for i,x in enumerate(labs)]
for n in range(9):
 for w in product((0,1),repeat=n):
  s=(1,0,2)
  for b in w:s=move(s,b)
  word=''.join(map(str,w))
  expect=word.startswith('01') and all(x=='1' for x in word[1:]) or word.startswith('101') and all(x=='1' for x in word[2:])
  check((s[2]-s[0]==2)==expect,(word,s))
  check((s[2]-s[0]==4)==(word.startswith('1101') and all(x=='1' for x in word[3:])),(word,s))
  check((s[2]-s[0]==6)==(word.startswith('001') and all(x=='1' for x in word[2:])),(word,s))
A=W('011')
B=A.copy();B[1]^=1<<3;B[3]^=1<<1
check(slot(A,4)==[11,6,0,0],A)
check(slot(B,4)==[13,6,0,0],B)
C=[110,72,2,4,15,30,60]
check(slot(C,0)==[41,20,0,0,0,0],C)
check(h([24,1,3,6,14])==11,'count 11')
for m in range(1,6):
 X=inflate(A,4,m);Y=inflate(B,4,m)
 check((h(X),h(Y))==(6*2**m-1,6*2**m+1),m)
 check(h(inflate(C,0,m))==20*2**m+1,m)
 w=W('0'+'1'*(m+1));check(w==X,(w,X))
print('PASS exact profiles, first-family fixed arc reversal, and defect-two word classification through length8')
D=[30,104,2,4,14,29,61]
check(slot(D,0)==[37,18,0,0,0,0],D)
check(h([40,1,3,6,15,30])==19,'count 19')
for m in range(1,6):check(h(inflate(D,0,m))==18*2**m+1,m)
print('PASS defects four and six word classification and literal defect-four readers')
