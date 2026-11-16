from verify import count, cyc, flip, h_word
from itertools import combinations

A=(160,145,259,519,205,94,1295,876,1083,375,703)
B=(12,145,1282,102,173,327,1175,813,603,127,955)
def main():
 for t,h in ((A,30011),(B,30013)):
  assert count(t)==h
  assert all(not(v>>i&1) for i,v in enumerate(t))
  assert all((t[i]>>j&1)+(t[j]>>i&1)==1 for i,j in combinations(range(11),2))
  assert len(cyc(t))==38
 scores=tuple(v.bit_count() for v in A)
 assert scores==tuple(v.bit_count() for v in B)==(2,3,3,4,5,5,6,6,6,7,8)
 s=tuple(range(14))+tuple(k+14 for k in scores)
 assert len(s)==25 and sum(s)==300
 assert sum(k*k for k in s)-70**2==-76
 edges=[]
 for w in ((0,5,3),(6,3,9,8,5,4),(0,7,6,8,10,2)):
  for u,v in zip(w,w[1:]+w[:1]):
   assert A[u]>>v&1
   edges.append((u,v))
 assert len({frozenset(e) for e in edges})==15
 assert flip(A,*edges)==B
 for j in range(2,31):
  k=2*j-4
  x=(5**j-1)//2
  U=format(x,'b');V=U[::-1].translate(str.maketrans('01','10'))
  W='0'+'1'*k+'0'
  lo=h_word(U+'01'+W+'10'+V);hi=h_word(U+'10'+W+'01'+V)
  offset=2*h_word(U)-x-1-(k==0)-(1<<k)
  assert lo==10**(2*j)+offset and hi==lo+2
  if j==2:assert lo==10007
  if j==3:assert lo==999999
 print('Different Hamiltonian counts, identical 25-score lattice vector:',s+(70,))
if __name__=='__main__':main()
