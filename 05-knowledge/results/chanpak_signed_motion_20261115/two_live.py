from check_motion import dp, Q, TT, tournament, isstrong, path_dp, lift, grades
from collections import defaultdict

def mask_shelf(A):
 n=len(A); f=(1<<n)-1; e=dp(A); s=dp(tuple(zip(*A))); h=[1]+list(map(sum,e[1:]))
 q=sum(h[x]*h[f^x] for x in range(f+1))
 rows=[]
 for m in range(f+1):
  sel=[v for v in range(n) if m>>v&1]
  mE=[sum(row[v] for v in sel) for row in e];mL=[sum(row[v] for v in sel) for row in s]
  rows.append((sum(mE[x]*h[f^x] for x in range(f+1)),sum(mL[x]*h[f^x] for x in range(f+1)),sum(mE[x]*mL[f^x] for x in range(f+1)),mE[f],mL[f]))
 return h[f],q,rows

def sweep(d,e,c):
 h1,q1,left=mask_shelf(d);h2,q2,right=mask_shelf(e);S=2**c
 shelves=defaultdict(dict)
 for i,(E,L,J,T,A) in enumerate(left):
  for j,(EE,LL,JJ,TT,AA) in enumerate(right):
   B=h2*E+S*q1*EE-S*q2*L-h1*LL
   C=h1*h2-h2*J-h1*JJ-S*EE*L-(T*AA if c==0 else 0)
   shelves[B][C]=(i,j)
 longest=0;best=None
 for B,row in shelves.items():
  for C in row:
   if C-2 in row:continue
   j=1
   while C+2*j in row:j+=1
   if j>longest:longest=j;best=(B,C,[row[C+2*i] for i in range(j)])
 return longest, best
if __name__=='__main__':
 for c in range(4):
  A=tournament(7,1464080);B=tournament(8,45442625)
  print(c,sweep(A,B,c))
