#!/usr/bin/env python3
"""Literal rooted paths and ordered two/three-path covers; standalone."""
import itertools,json
COORDS=['h','s','e','L','R','d','C','b','T']
def need(ok,why):
 if not ok:raise RuntimeError(why)
def profile(a):
 n=len(a);N=1<<n;F=N-1;end=[[0]*n for _ in range(N)];start=[[0]*n for _ in range(N)];H=[0]*N;H[0]=1
 for i in range(n):end[1<<i][i]=start[1<<i][i]=1
 for S in range(1,N):
  for i in range(n):
   if not(S>>i&1):continue
   rem=S^(1<<i)
   if rem:
    for j in range(n):
     if rem>>j&1:
      if a[j]>>i&1:end[S][i]+=end[rem][j]
      else:start[S][i]+=start[rem][j]
  H[S]=sum(end[S])
 beta=[0]*N
 for S in range(N):
  A=S
  while True:
   beta[S]+=H[A]*H[S^A]
   if not A:break
   A=(A-1)&S
 L=sum(end[S][0]*H[F^S]for S in range(N))
 R=sum(start[S][0]*H[F^S]for S in range(N))
 T=sum(start[S][0]*beta[F^S]for S in range(N))
 return(H[F],start[F][0],end[F][0],L,R,H[F^1],beta[F],beta[F^1],T)
def child(a,op):
 n=len(a);b=[x|(1<<n)for x in a]+[0]
 if op in'ZJ':b[0]^=1<<n;b[n]=1
 if op=='J':
  order=[n]+list(range(n));b=[sum(1<<j for j,v in enumerate(order)if b[u]>>v&1)for u in order]
 return b
def move(p,op):
 h,s,e,L,R,d,C,b,T=p
 if op=='O':return(h,s,0,L,2*R,d,2*C,2*b,3*T)
 if op=='Z':return(h+R-e,s,d,L+b,2*R-d,d,2*C+2*T-2*L,2*b,3*T-b)
 if op=='J':return(h+R-e,s,h-e,C-L,h+R,h,2*C+2*T-2*L,C,T+C)
 raise ValueError(op)
def graphs(n):
 es=list(itertools.combinations(range(n),2))
 for code in range(1<<len(es)):
  a=[0]*n
  for k,(i,j)in enumerate(es):
   if code>>k&1:a[i]|=1<<j
   else:a[j]|=1<<i
  yield a
counts={};trials=0;singleton_exception=None
for n in range(1,6):
 ct=0
 for a in graphs(n):
  p=profile(a);h,s,e,L,R,d,C,b,T=p
  need(C>=2*L and R>=d>=e and h>=s and b>=L and T>=b,'physical positivity')
  if n>=2:need(T>=b+s,'disjoint singleton and full-path cover classes')
  for op in'OZJ':
   expected=move(p,op);got=profile(child(a,op));trials+=1
   if n==1 and op=='Z':
    need(got[1]==0 and expected[1]==1,'singleton Z control')
    expected=expected[:1]+(0,)+expected[2:];singleton_exception=dict(before=p,after=got)
   need(got==expected,('literal transfer',n,a,op,p,expected,got))
   cv=(L,C-2*L,T-b,b-L)
   l,A,W,q=cv
   cm={'O':(l,2*l+2*A,3*W+l+q,l+2*q),'Z':(2*l+q,2*A+2*W,3*W,q),'J':(l+A,2*W+2*l+2*q,W+l+q,l)}[op]
   hh,ss,ee,ll,rr,dd,CC,bb,TT=got
   need(cm==(ll,CC-2*ll,TT-bb,bb-ll),'positive cover block')
   F=2*l+A+W+q;FF=2*cm[0]+cm[1]+cm[2]+cm[3]
   need(FF>=2*F,'factor-two positive cover growth')
   if n>=2 and s>0:
    B=h+L-s;BB=got[0]+got[3]-got[1]
    need(BB==B if op=='O'else BB>=B+2,'positive-start B growth')
  ct+=1
 counts[n]=ct
# Independently count one three-bin coordinate by assigning vertices to bins and permutation.
a=[6,4,0];p=profile(a);T_direct=0
for bins in itertools.product(range(3),repeat=3):
 if bins[0]!=0:continue
 vs=[[v for v in range(3)if bins[v]==b]for b in range(3)]
 def paths(vs,rooted=False):
  return sum((not rooted or w and w[0]==0)and all(a[u]>>v&1 for u,v in zip(w,w[1:]))for w in itertools.permutations(vs))
 T_direct+=paths(vs[0],True)*paths(vs[1])*paths(vs[2])
need(T_direct==p[-1]==9,'literal labeled three-bin control')
print(json.dumps(dict(status='PASS',coordinate_order=COORDS,labeled_root_zero_graphs_by_order=counts,literal_child_checks=trials,three_bin_permutation_control=T_direct,singleton_Z_exception=singleton_exception,scope='Exact all-input transfer checks through order5; all-order proof separate.'),sort_keys=True,indent=2))
