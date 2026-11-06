#!/usr/bin/env python3
from itertools import product
import json

def need(test,why):
 if not test:raise RuntimeError(why)

def core(n,code):
 r=[0]*n;k=0
 for i in range(n):
  for j in range(i+1,n):
   a,b=(i,j) if code>>k&1 else (j,i);r[a]|=1<<b;k+=1
 return r

def induced(r,vs):
 return [sum(1<<j for j,v in enumerate(vs) if r[u]>>v&1) for u in vs]

def counts(r):
 n=len(r);N=1<<n;d=[[0]*n for _ in range(N)];h=[0]*N;h[0]=1
 for i in range(n):d[1<<i][i]=1
 for S in range(1,N):
  for v in range(n):
   if not d[S][v]:continue
   p=r[v]&~S
   while p:
    b=p&-p;p-=b;d[S|b][b.bit_length()-1]+=d[S][v]
  h[S]=sum(d[S])
 return h,d

def hp(r):return counts(r)[0][-1]

def start(r,q=0):
 n=len(r);d=[[0]*n for _ in range(1<<n)];d[1<<q][q]=1
 for S in range(1,1<<n):
  for v in range(n):
   p=r[v]&~S
   if not d[S][v]:continue
   while p:
    b=p&-p;p-=b;d[S|b][b.bit_length()-1]+=d[S][v]
 return sum(d[-1])

def source(bits):
 n=len(bits)+1;r=[0]*n
 for i,b in enumerate(bits,1):
  if b:r[0]|=1<<i
  else:r[i]|=1
  for j in range(1,i):r[j]|=1<<i
 return r

def context(K,M,U,flip=0):
 n=len(K);r=K[:]+[0]*(len(U)-1)
 for u in range(1,len(U)):
  x=n+u-1
  if ((U[0]>>u)&1)^flip:r[0]|=1<<x
  else:r[x]|=1
  for c in range(1,n):
   if M>>c&1:r[x]|=1<<c
   else:r[c]|=1<<x
  for v in range(1,u):
   y=n+v-1
   if U[v]>>u&1:r[y]|=1<<x
   else:r[x]|=1<<y
 return r

def small_positive(U):
 # Root0; edge context root->a and triangle a->root->b->a.
 return context([2,0],0,U),context([4,1,2],0,U)

def main():
 sources=[source([1,a,b,0]) for a,b in product(range(2),repeat=2)]
 hs=[hp(U) for U in sources]
 need(hs==[9,9,11,9],('sources',hs))
 for U in sources:
  need(start(U)==1,'unit start')
  need(all(U[a]>>b&1 for a,b in zip([0,1,2,3,4],[1,2,3,4,0])),'spanning cycle')
 contexts=0;checks=0;quotients=set()
 for n in range(1,5):
  for code in range(1<<(n*(n-1)//2)):
   K=core(n,code)
   for M0 in range(1<<(n-1)):
    M=M0<<1;contexts+=1
    for flip in range(2):
     Ts=[context(K,M,U,flip) for U in sources]
     h=[hp(T) for T in Ts]
     x,y=n+1,n+2
     T=Ts[0]
     for z in range(len(T)):
      if z in (0,x,y):continue
      need((T[x]>>z&1)==(T[y]>>z&1),'literal twins')
     Q=induced(T,[z for z in range(len(T)) if z not in (0,y)])
     qh=hp(Q);quotients.add(qh)
     need(qh%2==1,'quotient parity')
     delta=h[3]-h[2]-h[1]+h[0]
     need(delta==-2*qh,('twin curvature',n,code,M,flip,h,qh))
     checks+=1
 # Independent exact positive-constructor verification.
 positive=0;maximum_order=0
 for n in range(1,6):
  for code in range(1<<(n*(n-1)//2)):
   U=core(n,code);hh,end=counts(U);h=hh[-1];s=start(U);full=(1<<n)-1
   L=sum(end[S][0]*hh[full^S] for S in range(1,1<<n) if S&1)
   P,C=small_positive(U)
   need(hp(P)==h+L-s,('edge',n,code))
   need(hp(C)==h+2*L,('triangle',n,code))
   need(hp(C)-2*hp(P)==2*s-h,('correction',n,code))
   positive+=1;maximum_order=max(maximum_order,len(C))
 print(json.dumps({'status':'PASS','four_source_counts_00_01_10_11':hs,'all_four_sources_strong_and_root_start_one':True,'labelled_contexts_core_orders_1_through_4':contexts,'curvature_checks_including_both_polarities':checks,'distinct_contracted_path_counts':sorted(quotients),'positive_identity_rooted_inputs_orders_1_through_5':positive,'positive_identity_max_output_order':maximum_order,'scope':'fixed one-copied-root uniform module contexts, arbitrary core size in proof; no unrestricted spectrum exclusion'},indent=2,sort_keys=True))
if __name__=='__main__':main()
