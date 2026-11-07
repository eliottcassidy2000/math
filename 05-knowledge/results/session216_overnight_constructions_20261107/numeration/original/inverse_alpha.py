from pathlib import Path
import json,sys
C=[tuple(map(int,l.split())) for l in (Path(__file__).parent/'CORE_STATES.txt').read_text().splitlines()]
def f(n):
 z=0
 while n:
  if n&1:z+=(n>>1)+((n>>1)&1)
  n>>=1
 return z
def find(N):
 tested=0
 for h in (5,9):
  for alpha in range(max(1,(N-h+1)//2),N-h+1):
   seed=alpha;b=0;r=0;rev=0
   while seed>0:
    for c,(hh,q,A,G,T,S,E) in enumerate(C):
     if hh!=h or seed<A or (seed-A)%q:continue
     p=(seed-A)//q
     if not p and not A:continue
     tau=p*G+h*f(p)+T+(p&1)*S+seed*rev+h*f(r)+E*((r>>(b-1))&1 if b else 0)
     tested+=1
     if 2*alpha+h-tau==N:return dict(status='FOUND',N=N,core=c,P=p,R=r,a=p.bit_length(),b=b,alpha=alpha,tau=tau,tested=tested)
    eps=seed&1;r|=eps<<b;rev=(rev<<1)|eps;seed=(seed-h*eps)//2;b+=1
 return dict(status='ABSENT_FROM_GRAMMAR',N=N,tested=tested)
if __name__=='__main__':
 for N in map(int,sys.argv[1:] or ['23','47','63','65','265','287','1237']):print(json.dumps(find(N)))
