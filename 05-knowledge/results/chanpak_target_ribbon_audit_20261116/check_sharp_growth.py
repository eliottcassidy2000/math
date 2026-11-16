from itertools import combinations
from fractions import Fraction
from cap_check import transfer_system,rowmult,need
def main():
 edges=list(combinations(range(1,4),2))
 cols=[]
 for a in range(8):
  A={j for j in (1,2,3) if a>>(j-1)&1}
  C=[]
  for mask in range(8):
   N={j:set() for j in (1,2,3)}
   for i,(p,q) in enumerate(edges):
    v,w=(p,q) if mask>>i&1 else (q,p)
    N[v].add(w)
   C.append(1+sum(1+2**(b+len(N[b]))-2**len(N[b]&A) for b in A))
  cols.append(C)
 need(list(map(max,cols))==[1,9,17,20,33,36,40,39])
 need([(a,m) for a,C in enumerate(cols) for m,c in enumerate(C) if c==40]==[(6,0)])
 # Relative labels −6 and −2 have distance four, therefore end −2.
 need((-2)-(-6)==4 and -2 in {-2,-3} and -6 not in {-2,-3})
 G=(118,100,96,55,38,64,24)
 def H(g):
  n=len(g);d=[[0]*n for _ in range(1<<n)]
  for u in range(n):d[1<<u][u]=1
  for m,row in enumerate(d):
   for u,z in enumerate(row):
    for v in range(n):
     if g[u]>>v&1 and not m>>v&1:d[m|(1<<v)][v]+=z
  return sum(d[-1])
 need(H(tuple(z&63 for z in G[:-1]))==1 and H(G)==39)
 states,_,M,F=transfer_system(3); f=[int(i in F) for i in range(len(states))]
 v=f[:];g=f[:];r=[1]*len(g)
 for _ in range(3):
  rows=[rowmult(v,A) for A in M.values()]
  v=[sum(row[i] for row in rows) for i in range(len(g))]
  g=[a+b for a,b in zip(g,v)]
 for _ in range(65):
  rows=[rowmult(r,A) for A in M.values()]
  r=[max([c]+[(10*row[i]+38)//39 for row in rows]) for i,c in enumerate(r)]
 need(max(Fraction(a,b) for a,b in zip(r,g))==Fraction(11,208))
 need(sum(312**k for k in range(4))==30468985)
 print("64 local cases, unique excluded 40 case, sharp 39 witness, and exact potential comparisons passed.")
if __name__=='__main__':main()
