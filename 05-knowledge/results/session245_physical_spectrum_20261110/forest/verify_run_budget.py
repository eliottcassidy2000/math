#!/usr/bin/env python3
"""Independent finite controls for the run-budget deduction. No density extrapolation."""
import itertools,json,math,sys

def need(ok,why):
 if not ok: raise RuntimeError(why)

def runs(w): return 0 if not w else 1+sum(x!=y for x,y in zip(w,w[1:]))

def hp(a):
 n=len(a);dp=[[0]*n for _ in range(1<<n)]
 for i in range(n):dp[1<<i][i]=1
 for mask in range(1,1<<n):
  for i in range(n):
   if mask>>i&1:
    p=mask^(1<<i)
    if p:dp[mask][i]=sum(dp[p][j] for j in range(n) if a[j]>>i&1)
 return sum(dp[-1]) if n else 1

def graph(edge,w):
 n=2+len(w);a=[0]*n;a[0 if edge else 1]|=1<<(1 if edge else 0)
 for t,s in enumerate(w,2):
  for r in range(2):
   if s>>r&1:a[r]|=1<<t
   else:a[t]|=1<<r
  for j in range(2,t):a[j]|=1<<t
 return a

def reduce_word(edge,w):
 a=graph(edge,w);reach=[x|(1<<i) for i,x in enumerate(a)];n=len(a)
 for v in range(n):
  for i in range(n):
   if reach[i]>>v&1:reach[i]|=reach[v]
 return tuple(s for t,s in enumerate(w,2) if any(i!=t and (reach[t]>>i&1) and (reach[i]>>t&1) for i in range(n)))

def L(h):
 result=0;a=0;v=1
 while v<=h:
  b=0;u=v
  while u<=h:result=max(result,2*a+3*b);b+=1;u*=5
  a+=1;v*=3
 return result

def main():
 if '--fail' in sys.argv:need(False,'intentional hostile sentinel')
 count_checks=0
 for q in range(2,5):
  for length in range(8):
   actual=[0]*(length+1)
   for n in range(length+1):
    for w in itertools.product(range(q),repeat=n):actual[runs(w)]+=1
   expected=[1]+[q*(q-1)**(j-1)*math.comb(length,j) for j in range(1,length+1)]
   need(actual==expected,('word count',q,length,actual,expected));count_checks+=1
 kernel_checks=0;deletions=0
 for edge in (False,True):
  for n in range(6):
   for w in itertools.product(range(4),repeat=n):
    reduced=reduce_word(edge,w);h=hp(graph(edge,w));hr=hp(graph(edge,reduced))
    need(h==hr,('SCC deletion',edge,w,reduced));need(len(reduced)<=L(h),('inherited kernel control',edge,w,h,reduced))
    need(runs(reduced)<=runs(w) and set(reduced)<=set(w),('run/alphabet preservation',w,reduced))
    kernel_checks+=1;deletions+=len(w)-len(reduced)
 print(json.dumps({'status':'PASS','exact_word_count_cases':count_checks,'literal_fixed_apex_kernel_cases':kernel_checks,'singleton_chain_vertices_deleted':deletions,'scope':'finite controls; universal count bound proved in RUN_BUDGET.md; inherited kernel credited'},indent=2,sort_keys=True))
if __name__=='__main__':main()
