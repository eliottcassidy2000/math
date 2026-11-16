from itertools import combinations
from oracle import path_dp,lift
def require(q):
 if not q: raise AssertionError
def energy(x):
 B=[i for i in range(x.bit_length()) if x>>i&1]
 return x*x+(x&1)+2*sum((1<<(i-j-1))+(i==j+1) for j,i in combinations(B,2))
def TT(n): return [[int(i<j) for j in range(n)] for i in range(n)]
def main():
 for n in range(1,15):
  E=[energy(x) for x in range(1<<n)]
  for x in range(len(E)-1):
   r=((x+1)&-(x+1)).bit_length()-1
   g=x>>(r+1)
   require(E[x+1]-E[x]==2*(1+3*g+(g&1)))
  require([x for x in range(len(E)-1) if E[x+1]-E[x]==2]==[(1<<r)-1 for r in range(n)])
 for n in range(1,5):
  for x in range(1<<n):
   V={v for v in range(n) if x>>(n-1-v)&1}
   W=lift([TT(1),TT(n),TT(n),TT(1)],[set(),V,{n-1-v for v in V},set()])
   require(path_dp(W)==2**(2*n+2)+1-energy(x))
 for n in range(1,12):
  def K(P):return (energy(P)-P*P-(P&1))//2
  def odd(P,t):return P*P+K(P)+t*(P+(P&1))
  for P in range((1<<n)-1):
   r=((P+1)&-(P+1)).bit_length()-1;g=P>>(r+1)
   require(odd(P+1,0)-odd(P,1)==1+3*g+(g&1))
  energies=sorted({odd(P,t) for P in range(1<<n) for t in (0,1)})
  edges=[i for i,(x,y) in enumerate(zip(energies,energies[1:])) if y-x==1]
  require(len(edges)==n and all(b>a+1 for a,b in zip(edges,edges[1:])))
 print('32,752 even carry checks, odd carry graph through 11 mirror pairs, and 30 literal tournament counts verified.')
if __name__=='__main__': main()
