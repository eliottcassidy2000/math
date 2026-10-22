from itertools import combinations,product
from math import factorial
roots=[]
for i,j in combinations(range(8),2):
 for a,b in product([-2,2],repeat=2):
  r=[0]*8;r[i]=a;r[j]=b;roots.append(tuple(r))
roots += [r for r in product([-1,1],repeat=8) if sum(v<0 for v in r)%2==0]
def count(R):
 adj=[{j for j in range(i+1,len(R)) if sum(a*b for a,b in zip(R[i],R[j]))==4} for i in range(len(R))]
 def cliq(c,k):
  if k==1:return len(c)
  return sum(cliq(c&adj[v],k-1) for v in c)
 return cliq(set(range(len(R))),5)
for name,R in [('E6',[r for r in roots if r[6]+r[7]==0 and r[5]-r[6]==0]),('E7',[r for r in roots if r[6]+r[7]==0]),('E8',roots)]:
 c=count(R);print(name,'roots',len(R),'unordered star frames',c,'A5 subsystems',c//12,flush=True)
