import math,itertools
def cycles(M):
 found=set()
 for a in range(M):
  trail=[];x=a
  while x not in trail:trail.append(x);x=2*x%M
  found.add(frozenset(trail[trail.index(x):]))
 return found
for M in range(2,201):
 r=M-1;k=2*r+1;owner={v:a for a in range(1,M) for v in (k-a,r-a)}
 for a in range(1,M):assert owner.get(k-2*a,0)==2*a%M
 cs=cycles(M); s=0;u=M
 while u%2==0:s+=1;u//=2
 assert sum(len(c) for c in cs)==u
 total=0
 for d in range(3,u+1,2):
  if u%d:continue
  phi=sum(math.gcd(v,d)==1 for v in range(d))
  period=next(t for t in range(1,d+1) if pow(2,t,d)==1)
  assert phi%period==0;total+=phi//period
 assert total==len(cs)-1
 for a in range(M):
  x=a;depth=0
  while x not in set().union(*cs):x=2*x%M;depth+=1
  val=s if a==0 else (a&-a).bit_length()-1
  assert depth==max(0,s-val)
 assert (len(cs)==1)==(u==1)
for M,C in [(21,5),(63,12),(341,37)]:assert len(cycles(M))-1==C
# Exhaustively inspect every order and every reversal assignment for small
# abstract incidence systems, using actual path-index triples. Labels are
# arbitrary fixed binary labels, as the ordering lemma allows.
for M in range(2,5):
 r=M-1;k=2*r+1; triples=[(k-a,r-a,k-2*a) for a in range(1,M)]
 C=len(cycles(M))-1;compatible=incompatible=0
 for order in itertools.permutations(range(k)):
  ranks={v:i for i,v in enumerate(order)}
  tri=[sorted(t,key=ranks.get) for t in triples]
  num=0
  for bits in range(1<<k):
   good=True
   for h,t in enumerate(tri):
    values=[((bits>>v)&1)^((v+h)%2) for v in t]
    if values[0]!=values[2] or values[0]==values[1]:good=False;break
   num+=good
  assert num in (0,2**(C+1))
  compatible+=num>0;incompatible+=num==0
 assert compatible>0 and (incompatible==0)==(M&(M-1)==0)
 print('all-order enumeration',M,'compatible',compatible,'incompatible',incompatible,'solutions',2**(C+1))
print('PASS: parent indices and CRT/depth through M200; examples M21,63,341; exhaustive order/reversal controls M2,3,4')
