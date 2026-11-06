import itertools,json,sys,pathlib
sys.path.insert(0,str(pathlib.Path(__file__).resolve().parents[1]))
from counting import nine,J,gate,labelled
from compiler import move

def brute(a):
 n=len(a);N=1<<n;F=N-1;h=[0]*N;end=[0]*N;start=[0]*N;h[0]=1
 for k in range(1,n+1):
  for p in itertools.permutations(range(n),k):
   if not all(a[u]>>v&1 for u,v in zip(p,p[1:])):continue
   m=sum(1<<i for i in p);h[m]+=1
   if p[-1]==0:end[m]+=1
   if p[0]==0:start[m]+=1
 beta=[0]*N
 for m in range(N):
  x=m
  while True:
   beta[m]+=h[x]*h[m^x]
   if not x:break
   x=(x-1)&m
 L=sum(end[m]*h[F^m]for m in range(N));R=sum(start[m]*h[F^m]for m in range(N));T=sum(start[m]*beta[F^m]for m in range(N))
 return h[F],start[F],end[F],L,R,h[F^1],beta[F],beta[F^1],T

def child(a,op):
 n=len(a)
 if op=='O':return tuple([x|1<<n for x in a]+[0])
 if op=='Z':return tuple([a[0]]+[x|1<<n for x in a[1:]]+[1])
 if op=='J':return J(a)
 return gate(a,int(op[-1]))
ct=0;exceptions=[]
for n in range(1,5):
 for a in labelled(n):
  p=brute(a);assert p==nine(a);ct+=1
  for op in ('O','Z','J','K1','K2'):
   b=child(a,op);q=brute(b);assert q==nine(b);ct+=1
   expect=move(p,op)
   if n==1 and op in ('Z','K1'):
    exceptions.append((op,expect[1],q[1]));expect=expect[:1]+(q[1],)+expect[2:]
   assert expect==q,(a,op,p,q,expect)
   h,s,e,L,R,d,C,bb,T=p;gq=q[6]-3*q[3]-q[0];g=C-3*L-h;tau=2*T+2*L-C-2*R-bb+2*d
   predicted={'O':2*g+3*L+h,'J':tau+(R-h+e-d)+(bb-L-d),'Z':tau+(R-h+e-d)+(bb-L-d)+3*(C-2*L-bb),'K1':tau+2*(C-L)+2*T-R-h+e,'K2':2*tau+5*(C-2*L)+12*T-3*h-2*R-2*d+2*e}[op]
   assert gq==predicted
print(json.dumps(dict(status='PASS',method='Independent enumeration of all directed vertex permutations, then literal bin convolution',scope='All labeled rooted inputs through n4 and all five literal children, including singleton exceptions',profiles=ct,singleton_exceptions=exceptions),indent=2))
