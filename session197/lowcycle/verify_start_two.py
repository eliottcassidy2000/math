#!/usr/bin/env python3
"""Standalone canonical graph controls for START_TWO.md. Explicit -O-safe gates."""
import itertools as it,json,hashlib
from pathlib import Path
P=Path(__file__).resolve().parent

def need(ok,why):
 if not ok:raise RuntimeError(why)
def relabel(r,order):return tuple(sum(1<<j for j,v in enumerate(order) if r[u]>>v&1) for u in order)
def strong(r):
 n=len(r);reach=[x|1<<u for u,x in enumerate(r)]
 for v in range(n):
  for u in range(n):
   if reach[u]>>v&1:reach[u]|=reach[v]
 return all(x==(1<<n)-1 for x in reach)
def words(r,root,cap=3):
 n=len(r);full=(1<<n)-1;out=[]
 def visit(w,used):
  if used==full:out.append(tuple(w));return
  todo=r[w[-1]]&~used
  while todo and len(out)<cap:
   bit=todo&-todo;todo-=bit;visit(w+[bit.bit_length()-1],used|bit)
 visit([root],1<<root);return out

def construct(U,entry,w,bits):
 m=len(U);n=m+w+1;b=n-1;r=list(U)+[0]*(w+1)
 for u in range(m):
  r[u]|=((1<<w)-1)<<m
  if u in entry:r[b]|=1<<u
  else:r[u]|=1<<b
 for i in range(w):
  v=m+i;r[v]|=((1<<(w-i-1))-1)<<(v+1)
  if bits>>i&1:r[b]|=1<<v
  else:r[v]|=1<<b
 return tuple(r),b

def cycle_key(r):
 n=len(r)
 return min(relabel(r,tuple(range(i,n))+tuple(range(i))) for i in range(n))
def rooted_two_key(r,root,ws=None):
 if ws is None:ws=words(r,root)
 need(len(ws)==2,'not exactly two rooted Hamilton paths')
 return min(relabel(r,w) for w in ws)
def labelled(n):
 pairs=list(it.combinations(range(n),2))
 for bits in range(1<<len(pairs)):
  r=[0]*n
  for j,(u,v) in enumerate(pairs):
   if bits>>j&1:r[u]|=1<<v
   else:r[v]|=1<<u
  yield tuple(r)

def main():
 F=[0,1]
 for _ in range(205):F.append(F[-1]+F[-2])
 rootone={1:{(0,)},2:set()};unroot={};pairs={};one_stats=[];one_physical=0
 for n in range(3,10):
  bank=set()
  for m in range(1,n):
   w=n-m-1
   if m==1 and w==0:continue
   for U in rootone[m]:
    for bits in range(1<<max(0,w-1)):
     r,b=construct(U,(0,),w,bits);r=relabel(r,(b,)+tuple(range(b)))
     need(r not in bank,'root-one decomposition collision');bank.add(r)
     ws=words(r,0);need(len(ws)==1 and ws[0]==tuple(range(n)) and r[n-1]&1,'literal root-one cycle');one_physical+=1
  rootone[n]=bank;need(len(bank)==F[2*n-4],'root-one enumeration')
  unroot[n]={cycle_key(r) for r in bank};need(len(unroot[n])==(1 if n==3 else F[2*n-6]),'unrooted C1 count')
  hist={};pairbank=[]
  for U in sorted(unroot[n]):
   S={u for u in range(n) if len(words(U,u))==1};low={u for u in range(n) if U[u].bit_count()==1};t=len(S)
   if n==3:need(t==3 and len(low)==3,'triangle boundary')
   else:
    need(t in (2,3),'start-one block size')
    starts=[i for i in range(n) if S=={(i+j)%n for j in range(t)}]
    need(len(starts)==1,'nonconsecutive start-one set');i=starts[0]
    need(low=={(i+j)%n for j in range(t-1)},'block outdegree pattern')
    need(all(relabel(U,tuple(range(j,n))+tuple(range(j)))!=U for j in range(1,n)),'C1 rigidity')
   hist[t]=hist.get(t,0)+1
   marks=list(it.combinations(sorted(S),2))
   if n==3:marks=marks[:1]
   pairbank.extend((U,pair) for pair in marks)
  pairs[n]=pairbank;d=1 if n==3 else F[2*n-8]+3*F[2*n-7]
  need(len(pairbank)==d,'unordered-pair bank')
  if n>=4:need(hist.get(2,0)==F[2*n-8] and hist.get(3,0)==F[2*n-7],'refined Fibonacci counts')
  one_stats.append(dict(order=n,rooted_start_one=len(bank),unrooted_C1=len(unroot[n]),start_one_size_histogram=hist,unordered_start_one_pairs=d))
 # All-order arithmetic identities, independent recurrence versus denominator.
 d=[0]*101;d[3]=1
 for n in range(4,101):d[n]=F[2*n-8]+3*F[2*n-7]
 B=[0]*101
 for n in range(4,101):
  B[n]=sum((1 if n-m==1 else 2**(n-m-2))*(B[m]+d[m]) for m in range(3,n))
  if n>=8:need(B[n]==6*B[n-1]-11*B[n-2]+6*B[n-3]-B[n-4],'start-two rational recurrence')
 need(B[4:11]==[1,5,18,60,191,589,1775],'initial start-two coefficients')
 root_two={n:set() for n in range(1,11)};two_stats=[];two_physical=0
 for n in range(4,11):
  for m in range(3,n):
   w=n-m-1;inputs=[(U,(0,)) for U in root_two[m]]+pairs[m]
   for U,entry in inputs:
    for bits in range(1<<max(0,w-1)):
     r,b=construct(U,entry,w,bits);ws=words(r,b);key=rooted_two_key(r,b,ws)
     need(key not in root_two[n],'start-two canonical decomposition collision');root_two[n].add(key)
     need(strong(r),'start-two construction not strong');need(sum(bool(r[w[-1]]>>b&1) for w in ws) in (1,2),'cycle count outside1/2');two_physical+=1
  need(len(root_two[n])==B[n],'start-two enumeration')
  hist={1:0,2:0}
  for r in root_two[n]:hist[sum(bool(r[w[-1]]&1) for w in words(r,0))]+=1
  digest=hashlib.sha256(json.dumps(sorted(root_two[n]),separators=(',',':')).encode()).hexdigest()
  two_stats.append(dict(order=n,rooted_start_two=len(root_two[n]),cycle_count_histogram=hist,adjacency_bank_sha256=digest))
 complete=[];labelled_count=0
 for n in range(1,7):
  actual=set()
  for r in labelled(n):
   labelled_count+=1
   if not strong(r):continue
   for root in range(n):
    ws=words(r,root)
    if len(ws)==2:actual.add(rooted_two_key(r,root,ws))
  need(actual==root_two[n],'complete labelled start-two inventory')
  complete.append(dict(order=n,rooted_start_two=len(actual)))
 result=dict(status='PASS',root_one_physical_graphs=one_physical,root_one_inventory=one_stats,root_two_physical_graphs=two_physical,root_two_inventory=two_stats,complete_labelled_tournaments=labelled_count,complete_labelled_max_order=6,complete_root_two_inventory=complete,recurrence_checked_through_order=100)
 data=json.dumps(result,sort_keys=True,indent=2)+'\n';(P/'start_two_certificate.json').write_text(data);print(data,end='')
if __name__=='__main__':main()
