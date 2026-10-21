import json, math, itertools
from collections import defaultdict
from verify_pathcovers import need,subset_c,brute_fixed

def circulate(n,S):return [[int((j-i)%n in S) for j in range(n)] for i in range(n)]
def walk_contract(A,c):
 q=len(A);m=len(c)-1;states={((0,)*q,-1):1};ans=0
 for length in range(q*m+1):
  nxt=defaultdict(int)
  for (counts,last),v in states.items():
   if last>=0:ans+=v*math.prod(c[r] for r in counts)
   for j in range(q):
    if counts[j]<m and (last<0 or A[last][j]):
     cc=list(counts);cc[j]+=1;nxt[tuple(cc),j]+=v
  states=nxt
 return ans

def fixed_contract(A,c,d):
 q=len(A);m=len(c)-1;p=(q-1)//2;pair=[0]+[min(i,q-i) for i in range(1,q)];caps=[(m-1)//2]+[m]*p
 states={((0,)*(p+1),-1):1};ans=0
 for length in range(sum(caps)+1):
  nxt=defaultdict(int)
  for (counts,last),v in states.items():
   if last<0 or A[last][0]:ans+=v*d[2*counts[0]+1]*math.prod(c[r] for r in counts[1:])
   for j in range(q):
    t=pair[j]
    if counts[t]<caps[t] and (last<0 or A[last][j]):
     cc=list(counts);cc[t]+=1;nxt[tuple(cc),j]+=v
  states=nxt
 return ans

def main():
 rows=[]
 for n,S,a0,e0 in [(3,{1},3,1),(5,{1,2},5,1),(7,{1,2,4},21,3)]:
  A=circulate(n,S);c=subset_c(A);d=brute_fixed(A,[(-i)%n for i in range(n)])
  aut=[p for p in itertools.permutations(range(n)) if all(A[p[i]][p[j]]==A[i][j] for i in range(n) for j in range(n))]
  central=[p for p in aut if all(p[(-i)%n]==(-p[i])%n for i in range(n))]
  need(len(aut)==a0 and len(central)==e0,'independent base groups')
  need(walk_contract(A,[0,1])==c[1],'base H');need(fixed_contract(A,[0,1],[0,1])==d[1],'base d')
  H=walk_contract(A,[0,3,6,6]);D=fixed_contract(A,[0,3,6,6],[0,1,0,2]);a=a0*3**n;e=e0*3**((n-1)//2)
  need(H%a==0 and D%e==0,'seed integrality');h=H//a;b=D//e;need((h-b)%2==0,'seed parity')
  rows.append(dict(base=n,base_c=c,base_d=d,n=3*n,H=H,d1=D,aut=a,centralizer=e,h=h,b=b,paired_black_degree=(h-b)//2))
 print(json.dumps(rows,indent=2))
if __name__=='__main__':main()
