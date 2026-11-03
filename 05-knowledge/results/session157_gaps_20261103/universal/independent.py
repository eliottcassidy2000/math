from math import gcd,lcm

def factors(n):
 p=3;z=[]
 while p*p<=n:
  if n%p==0:
   z.append(p)
   while n%p==0:n//=p
  p+=2
 if n>1:z.append(n)
 return z

def Lvalue(m):
 L=1
 for p in factors(m):
  z=1
  while pow(2,z,p)!=1:z+=1
  L=lcm(L,z)
 return L

def hword(w):
 H=2;B=0;prev=None
 for x in w:
  if not x:
   H+=B
  elif prev==0:H+=1
  B=2*B+x;prev=x
 return H

def graph_count(w):
 N=len(w)+1;O=[0]*N;r=N-1
 for i,b in enumerate(w):
  for j in range(i+1,len(w)):O[i]|=1<<j
  if b:O[r]|=1<<i
  else:O[i]|=1<<r
 D=[[0]*N for _ in range(1<<N)]
 for i in range(N):D[1<<i][i]=1
 for s in range(1,1<<N):
  for i,c in enumerate(D[s]):
   if c:
    z=O[i]&~s
    while z:
     t=z&-z;z-=t;D[s|t][t.bit_length()-1]+=c
 return sum(D[-1])
for k in (2,3,4):
 for n in range(3):
  for a in range(1,3):
   w=[1]+([0]*(k-1)+[1])*n+[1,0]*a
   if len(w)<=12 and hword(w)!=graph_count(w): raise AssertionError((k,n,a))
print('Independent graph DP agrees with transfer for 17 small motifs')
for e in range(5):
 for m in range(1,46,2):
  M=2**e*m;L=Lvalue(m);k=1
  while k<e+1:k+=L
  q=pow(2,k,M) if M>1 else 0;values=set()
  for a in range(1,2**max(e,1)+e+3):
   H=2%M;B=1%M
   for n in range(m*L):
    G=H;Z=B;prev=1
    for x in [1,0]*a:
     if x==0:G=(G+Z)%M
     elif prev==0:G=(G+1)%M
     Z=(2*Z+x)%M;prev=x
    values.add(G)
    H=(H+(pow(2,k-1,M)-1)*B+1)%M;B=(q*B+1)%M
  want={r for r in range(M) if not(M%2==0 and r%2==0)}
  if values!=want: raise AssertionError((e,m,k,len(values),len(want),sorted(want-values)[:5]))
print('All parity-eligible residues confirmed for e=0..4 and odd m<=45')
