from itertools import combinations
from math import isqrt
A=(578,1164,1449,673,1039,659,1086,1105,1787,1238,41)
B=(578,1164,1449,689,1063,643,1086,1105,1787,1238,41)
def count(R):
 n=len(R);N=1<<n
 dp=[[0]*n for _ in range(N)]
 for i in range(n):dp[1<<i][i]=1
 for m,v in enumerate(dp):
  for i,x in enumerate(v):
   if not x:continue
   e=R[i]&~m
   while e:
    b=e&-e;e^=b
    dp[m|b][b.bit_length()-1]+=x
 return sum(dp[-1])
def flip(R,*edges):
 r=list(R)
 for i,j in edges:r[i]^=1<<j; r[j]^=1<<i
 return tuple(r)
def cyc(r):
 n=len(r)
 return [((1<<i)|(1<<j)|(1<<k)) for i,j,k in combinations(range(n),3)
         if ((r[i]>>j&1)==(r[j]>>k&1)==(r[k]>>i&1))]
def contract(R,chain):
 u,v,w=chain
 assert R[u]>>v&1 and R[v]>>w&1
 keep=[j for j in range(len(R)) if j not in chain]+[u]
 return [sum((((R[w] if j==u else R[j])>>k)&1)<<q
             for q,k in enumerate(keep) if k!=j) for j in keep]

def contract(r, chain):
 vs=[v for v in range(len(r)) if v not in chain]
 assert all(r[a]>>b&1 for a,b in zip(chain,chain[1:]))
 out=[]
 for a in vs:
  out.append(sum(((r[a]>>b)&1)<<j for j,b in enumerate(vs)) |
             (((r[a]>>chain[0])&1)<<len(vs)))
 out.append(sum(((r[chain[-1]]>>b)&1)<<j for j,b in enumerate(vs)))
 return tuple(out)

def h_word(w):
 h,x,p=1,0,0
 for b in map(int,w):
  if b==0:h+=x+p
  p=b;x=2*x+b
 return h

def word_rows(w):
 n=1+len(w)
 d=[((1<<n)-1)^((1<<(j+1))-1) for j in range(n)]
 for j,b in enumerate(w,1):
  if b=='0':
   d[0]^=1<<j;d[j]^=1
 return tuple(d)

def join(x,y,cut=None):
 n,m=len(x),len(y)
 d=[r|(((1<<m)-1)<<n) for r in x]+[r<<n for r in y]
 if cut:
  u,v=cut
  d[u]^=1<<(n+v);d[n+v]^=1<<u
 return tuple(d)

def main():
 for R,t in ((A,30011),(B,30013)):
  assert count(R)==t
  assert all(t%d for d in range(2,isqrt(t)+1))
  assert all(not ((r>>i)&1) for i,r in enumerate(R))
  assert all((R[i]>>j&1)+(R[j]>>i&1)==1 for i,j in combinations(range(11),2))
  c=cyc(R)
  print(t,'minimum deletions to transitive', next(d for d in range(12) if any(all(s&v for v in c) for s in range(1<<11) if s.bit_count()==d)))
 e,f=(3,4),(4,5)
 assert flip(A,e,f)==B
 assert count(contract(A,(5,4,3)))==1007
 assert count(contract(B,(3,4,5)))==1485
 print('individual toggles and both:',*(count(flip(A,*es)) for es in ((e,),(f,),(e,f))))
 assert count(contract(A,(5,4,3)))==1007
 assert count(contract(B,(3,4,5)))==1485
 print('chain terms:',1007,1485,'mixed:',2492)
 missing={29989,30011,30013}
 for n in range(2,17):
  for v in range(1<<(n-2)):
   w='1'+format(v,'0'+str(n-2)+'b')+'0' if n>2 else '10'
   assert h_word(w) not in missing
 print('All orders: none of 29989, 30011, 30013 has a one-apex witness.')
 T=(2,4,1)
 words=[('1000011001100',0,30011),('10001000110110',1,30013)]
 for w,y,h in words:
  assert count(join(T,word_rows(w),(0,y)))==h
 print('Crossed-triangle certificates:',[h for _,_,h in words])
 U='1100';V=U[::-1].translate(str.maketrans('01','10'))
 assert [count(join(T,word_rows(U+W+V))) for W in ('010010','100001')]==[30021,30027]
 print('Ordinal-triangle scaling gives 30021/30027.')
 adj=[(count(flip(A,e)),e) for e in combinations(range(11),2)]
 assert all(h!=30013 for h,e in adj)
 print('No single reverse hits 30013 at this witness; closest', sorted(adj,key=lambda x:abs(x[0]-30013))[:5])
if __name__=='__main__':main()
