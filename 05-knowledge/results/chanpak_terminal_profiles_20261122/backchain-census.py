from collections import defaultdict,Counter
from sys import argv
N=int(argv[1]) if len(argv)>1 else 100001
def odds(n):return (n&1)==1
def sieve(n):
 P=bytearray(b'\x01')*(n+1);P[0:2]=b'\x00\x00'
 for p in range(2,int(n**.5)+1):
  if P[p]:P[p*p::p]=b'\x00' * ((n-p*p)//p+1)
 return P
seen=set()
todo=[]
w=2
while w+1<=N:todo.append((w,0,1));w*=2
Fib=defaultdict(set);Full=defaultdict(int)
FiPrime=Counter();P=sieve(N)
while todo:
 d,u,v=todo.pop()
 s=(d,u,v)
 if s in seen:continue
 seen.add(s)
 h=d+u+v
 assert odds(h)
 Full[h]+=1
 if not u: Fib[h].add(v)
 z=d+2*u,v,d+u
 if sum(z)<=N:todo.append(z)
 a=h+u
 w=2
 while h+w*a<=N:
  todo.append((w*a,0,h));w*=2
Full[1]+=1
print('cutoff',N,'all states',len(seen)+1,'heights',len(Full),'positive states',sum(map(len,Fib.values())),'positive heights',len(Fib))
print('full fiber maxima',sorted(Full.items(),key=lambda x:x[1],reverse=True)[:8])
print('positive fiber maxima',[(h,len(gs),sorted(gs)[:20]) for h,gs in sorted(Fib.items(),key=lambda x:len(x[1]),reverse=True)[:8]])
print('prime heights',sum(P[h] for h in Full),'of','all primes',sum(P))
