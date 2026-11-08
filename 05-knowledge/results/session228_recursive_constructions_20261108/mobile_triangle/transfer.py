from itertools import permutations,product
ST=list(product(range(3),range(2))) # P0empty P1in P2out; R0empty R1good

def local(block,k):
 n=len(block); full=(1<<n)-1; out=set(range(k)); a=[]
 for s in range(1<<n):
  v=[i for i in range(n) if s>>i&1];w=[i for i in range(n) if not(s>>i&1)]
  for p in permutations(v):
   if any(not(block[u]>>v&1) for u,v in zip(p,p[1:])):continue
   for r in permutations(w):
    if any(not(block[u]>>v&1) for u,v in zip(r,r[1:])):continue
    a.append((0 if not p else (1 if p[-1] not in out else 2),0 if not r else (1 if r[0] in out else 2)))
 return a

def matrix(block,k):
 M=[[0]*6 for _ in ST]
 for j,(p,r) in enumerate(ST):
  for pp,rr in local(block,k):
   if not r and rr==2:continue
   new=(pp or p,int(bool(r or rr)))
   M[ST.index(new)][j]+=1
 return M

def mm(A,B):return [[sum(x*y for x,y in zip(row,col)) for col in zip(*B)] for row in A]
def mv(A,x):return [sum(a*b for a,b in zip(row,x)) for row in A]
v=[1,0,0,0,0,0]; ell=[1,1,1,1,0,0]
S=[matrix([0],b) for b in range(2)];T=[matrix([2,4,1],k) for k in range(4)]
def dvalue(bits,p):
 x=v[:];i=0
 while i<len(bits):
  if i==p:x=mv(T[2*bits[i]+bits[i+1]],x);i+=2
  else:x=mv(S[bits[i]],x);i+=1
 return sum(a*b for a,b in zip(ell,x))
if __name__=='__main__':
 import sympy as s
 for b in range(2): print('S',b,s.Matrix(S[b]))
 for k in (1,2):
  F=s.Matrix(T[k])-3*s.Matrix(S[k%2])*s.Matrix(S[k//2]); print('F',k,F,'rank',F.rank());print(F.rank_decomposition())
