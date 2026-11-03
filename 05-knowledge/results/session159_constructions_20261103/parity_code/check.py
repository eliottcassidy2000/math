from itertools import combinations, permutations
from math import comb

def check(A,B,C):
 n=A+B+C+1; a=0;c=A;b=A+B;d=A+B+C
 adj=[[i<j for j in range(n)] for i in range(n)]
 for i,j in [(a,b),(c,d)]:adj[i][j]=False;adj[j][i]=True
 h=[[0]*n for _ in range(1<<n)]
 for v in range(n):h[1<<v][v]=1
 for mask in range(1,1<<n):
  for v in range(n):
   if h[mask][v]:
    for u in range(n):
     if not(mask>>u&1) and adj[v][u]: h[mask|1<<u][u]+=h[mask][v]
 H=sum(h[-1]); cycles=[]
 for size in range(3,n+1,2):
  for S in combinations(range(n),size):
   v=S[0]
   for tail in permutations(S[1:]):
    seq=(v,)+tail
    if all(adj[seq[i]][seq[(i+1)%size]] for i in range(size)):
     cycles.append(sum(1<<x for x in S))
 N=len(cycles); P=sum((u&v)==0 for i,u in enumerate(cycles) for v in cycles[i+1:])
 formulaP=(2**(A+C-2)*3**(B-1)-(C==1)*2**(A-1)-(A==1)*2**(C-1)+(A==C==1)*(-1)**(B-1))//4
 formulaH=1+2**(A+B-1)+2**(B+C-1)-(C==1)*2**(A-1)-(A==1)*2**(C-1)+2**(A+C-1)*3**(B-1)
 if (P,H)!=(formulaP,formulaH): raise ValueError((A,B,C,P,formulaP,H,formulaH))
 if H!=1+2*N+4*P: raise ValueError((A,B,C,H,N,P))
 return N,P,H
count=0
for A in range(1,7):
 for B in range(1,7):
  for C in range(1,7):
   if A+B+C+1<=9:
    N,P,H=check(A,B,C);count+=1
    if (A,B,C) in [(1,2,1),(1,3,1),(2,2,2)]:print((A,B,C),'N,P,H=',N,P,H)
print('PASS triples',count,'max vertices 9')
