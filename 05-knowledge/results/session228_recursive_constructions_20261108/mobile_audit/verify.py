#!/usr/bin/env python3
"""Independent literal DP and six-state cover audit; no producer imports."""
from itertools import product, permutations
from pathlib import Path
import hashlib,json

def ck(v,m):
 if not v: raise RuntimeError(m)
def graph(w,p):
 groups=[];outs=[];i=0
 while i<len(w):
  size=3 if i==p else 1
  groups.append(size);outs.append(2*w[i]+w[i+1] if size==3 else w[i]);i+=2 if size==3 else 1
 n=1+sum(groups);a=[[False]*n for _ in range(n)];offset=1
 for size,k in zip(groups,outs):
  for j in range(size):
   u=offset+j;a[0][u]=j<k;a[u][0]=j>=k
   for v in range(offset+size,n):a[u][v]=True
  if size==3:
   for j in range(3):a[offset+j][offset+(j+1)%3]=True
  offset+=size
 ck(all(a[i][j]!=a[j][i] for i in range(n) for j in range(i)), 'tournament')
 return a

def profile(a):
 n=len(a);full=(1<<n)-1;dp=[[0]*n for _ in range(full+1)];st=[[0]*n for _ in range(full+1)]
 for v in range(n):dp[1<<v][v]=1
 st[1][0]=1
 for mask in range(1,full+1):
  for v in range(n):
   if mask>>v&1:
    for u in range(n):
     if not(mask>>u&1) and a[v][u]:dp[mask|1<<u][u]+=dp[mask][v];st[mask|1<<u][u]+=st[mask][v]
 h=[sum(t) for t in dp];h[0]=1
 return h[full],dp[full][0],sum(sum(st[s])*h[full^s] for s in range(full+1)),h[full^1]

# State: last pre-root endpoint category empty/in/out; post-root path begun.
states=list(product(range(3),range(2)))
def matrix(n,k):
 a=[[i!=j and (n==1 or (i+1)%3==j) for j in range(n)] for i in range(n)]
 M=[[0]*6 for _ in range(6)]
 for assignment in product(range(2),repeat=n):
  for P in permutations([i for i in range(n) if assignment[i]==0]):
   for Q in permutations([i for i in range(n) if assignment[i]==1]):
    if any(not a[u][v] for path in (P,Q) for u,v in zip(path,path[1:])):continue
    for col,(last,begun) in enumerate(states):
     if Q and not begun and Q[0]>=k:continue
     newlast=(1 if P[-1]>=k else 2) if P else last
     row=states.index((newlast,int(begun or bool(Q))))
     M[row][col]+=1
 return M
def mm(A,B):return [[sum(a*b for a,b in zip(row,col)) for col in zip(*B)] for row in A]
def mv(A,v):return [sum(a*b for a,b in zip(row,v)) for row in A]
def sub(A,B,c=1):return [[a-c*b for a,b in zip(ar,br)] for ar,br in zip(A,B)]
I=[[int(i==j) for j in range(6)] for i in range(6)];S=[matrix(1,k) for k in range(2)];T=[matrix(3,k) for k in range(4)]
v=[1,0,0,0,0,0];ell=[1,1,1,1,0,0]
def count(w,p):
 x=v;i=0
 while i<len(w):
  x=mv(T[2*w[i]+w[i+1]] if i==p else S[w[i]],x);i+=2 if i==p else 1
 return sum(a*b for a,b in zip(ell,x))
def formula(r,p):return (4**r+6*r-1)//3+(0 if p==0 else (-2)**p-2)
def main():
 literal=trailing=0;boundary=[]
 for m in range(2,8):
  for w in product(range(2),repeat=m):
   for p in range(m-1):
    D,E,U,d=profile(graph(w,p));literal+=1
    ck(D==count(w,p),'literal D versus independent cover transfer')
    ck(U==3*(1+int(''.join(map(str,w)),2)) and d==3,'binary U and deletion')
    if w[-1] and p<m-2:ck(E==0,'trailing singleton E');trailing+=1
    if m==2:boundary.append([list(w),D,E,U,d])
 A=mm(S[1],S[0]);poly=mm(mm(mm(A,sub(A,I,4)),sub(A,I)),sub(A,I))
 ck(not any(map(any,poly)),'annihilator')
 # 4x4 grid in independent exponents: establishes recurrence proof.
 grids={}
 for parity in (0,1):
  grid=[]
  for j in range(4):
   row=[]
   for k in range(4):
    r=j+k+2;p=2*j+parity;D=count([0,1]*r,p)
    ck(D==formula(r,p),'initial grid');row.append(D)
   grid.append(row)
  grids[str(parity)]=grid
 large=0
 for r in range(2,45):
  ds=[count([0,1]*r,p) for p in range(2*r-2)]
  ck(ds==[formula(r,p) for p in range(2*r-2)],'large formulas')
  ck(len(set(ds))==2*r-2,'distinctness');large+=len(ds)
 print(json.dumps({'status':'PASS','literal_all_words_lengths':[2,7],'literal_profiles':literal,'trailing_singleton_checks':trailing,'terminal_triangle_boundary_profiles':boundary,'annihilator':'A(A-4I)(A-I)^2 = 0','initial_grids':grids,'large_transfer_checks':large,'max_r':44,'S':S,'T':T,'scope':'Physical fibers only; no interval coverage or recursive radix claim.'},sort_keys=True,indent=2))
if __name__=='__main__':main()
