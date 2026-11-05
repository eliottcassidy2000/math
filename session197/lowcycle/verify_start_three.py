#!/usr/bin/env python3
"""Physical canonical start-three grammar and independent small-universe controls."""
from pathlib import Path
import hashlib,json,itertools as it
P=Path(__file__).resolve().parent
needhash='ca788fa6489cbca7bed8f9f9a84a553ff9f5de7502833544f9b51b9f059142ae'
if hashlib.sha256((P/'verify_start_two.py').read_bytes()).hexdigest()!=needhash:raise RuntimeError('start-two dependency changed')
from verify_start_two import need,relabel,strong,words,construct,cycle_key,labelled

def key3(r,root):
 ws=words(r,root,4);need(len(ws)==3,'not exactly three starts')
 return min(relabel(r,w)for w in ws)
def tail_triangle(U,i,j,smask,tmask):
 m=len(U);k=i+3+j;n=m+k+1;root=n-1;r=list(U)+[0]*(k+1)
 for u in range(m):
  r[u]|=((1<<k)-1)<<m
  if u==0:r[root]|=1<<u
  else:r[u]|=1<<root
 tri=m+i
 for a in range(m,m+k):
  for b in range(a+1,m+k):
   if tri<=a<b<tri+3:
    if (a-tri,b-tri)in ((0,1),(1,2)):r[a]|=1<<b
    else:r[b]|=1<<a
   else:r[a]|=1<<b
 singles=list(range(m,tri))+list(range(tri+3,m+k))
 for h,a in enumerate(singles):
  if smask>>h&1:r[root]|=1<<a
  else:r[a]|=1<<root
 for h in range(3):
  if tmask>>h&1:r[root]|=1<<(tri+h)
  else:r[tri+h]|=1<<root
 return tuple(r),root

def whole_key(r):
 n=len(r);cycles=[];full=(1<<n)-1
 def visit(w,mask):
  if mask==full:
   if r[w[-1]]&1:cycles.append(tuple(w))
   return
  todo=r[w[-1]]&~mask
  while todo:
   b=todo&-todo;todo-=b;visit(w+[b.bit_length()-1],mask|b)
 visit([0],1);need(len(cycles)==3,'whole graph not C3')
 return min(relabel(r,w[i:]+w[:i])for w in cycles for i in range(n))

def main():
 F=[0,1]
 for _ in range(210):F.append(F[-1]+F[-2])
 def fib(k):return 1 if k==-1 else F[k]
 rootone={1:{(0,)},2:set()};unroot={};mixed={};triples={};clone_stats=[]
 for n in range(3,10):
  bank=set()
  for m in range(1,n):
   w=n-m-1
   if m==1 and w==0:continue
   for U in rootone[m]:
    for bits in range(1<<max(0,w-1)):
     r,b=construct(U,(0,),w,bits);r=relabel(r,(b,)+tuple(range(b)));bank.add(r)
  rootone[n]=bank;unroot[n]={cycle_key(r)for r in bank}
  mixed[n]=[];triples[n]=[];root2_byS={2:0,3:0};clone_count=0
  for U in unroot[n]:
   ws=[words(U,u,4)for u in range(n)];S1=[u for u in range(n)if len(ws[u])==1];S2=[u for u in range(n)if len(ws[u])==2]
   if n>=4:root2_byS[len(S1)]+=len(S2)
   mixed[n].extend((U,(a,b))for a in S1 for b in S2)
   if len(S1)==3:triples[n].append((U,tuple(S1)))
   # Every cycle-consecutive unique-start pair is a physical transitive clone.
   seen=set()
   for a in S1:
    b=(a+1)%n
    if b not in S1:continue
    r,root=construct(U,(a,b),0,0);need(strong(r)and len(words(r,root,4))==2,'clone start-two')
    need(sum(bool(r[w[-1]]>>root&1)for w in words(r,root,4))==1,'clone C1')
    need(U[a].bit_count()==1,'pair first not degree1')
    for x in range(n):
     if x!=a:need(bool(r[root]>>x&1)==bool(U[a]>>x&1),'clone external incidence')
    k=min(relabel(r,w)for w in words(r,root,4));seen.add(k)
   clone_count+=len(seen)
  if n>=4:
   need(root2_byS=={2:F[2*n-8],3:fib(2*n-9)},'refined C1 start-two count')
   need(len(mixed[n])==2*F[2*n-8]+3*fib(2*n-9),'mixed pair inventory')
  clone_stats.append({'order':n,'mixed_pairs':len(mixed[n]),'start_two_by_start_one_size':root2_byS,'consecutive_pair_clone_outputs':clone_count})
 # Exact coefficient recurrence from the canonical grammar, checked independently against the claimed denominator.
 B1=[0]*101;B1[1]=1
 for n in range(3,101):B1[n]=F[2*n-4]
 M=[0]*101;T=[0]*101;T[3]=1
 for n in range(4,101):M[n]=2*F[2*n-8]+3*fib(2*n-9);T[n]=F[2*n-7]
 B3=[0]*101
 for n in range(4,101):
  B3[n]=sum((1 if n-m==1 else 2**(n-m-2))*(B3[m]+M[m]+T[m])for m in range(3,n))-(n==4)
  B3[n]+=sum(B1[m]*(2*(n-m-1)-3)*2**(n-m-4)for m in range(1,n-3))
 for n in range(101):
  lhs=B3[n]-sum(c*B3[n-j]for j,c in ((1,6),(2,-11),(3,6),(4,-1))if n>=j)
  need(lhs==({5:8,6:-17,7:8}.get(n,0)),'start-three rational GF')
 need(B3[5:9]==[8,31,106,343],'initial start-three coefficients')
 banks={n:set()for n in range(1,11)};inventory=[];branches=[0,0]
 for n in range(4,11):
  for m in range(3,n):
   w=n-m-1
   entries=[(U,(0,))for U in banks[m]]+mixed[m]+triples[m]
   for U,entry in entries:
    if m==3 and len(entry)==3 and w==0:continue
    for bits in range(1<<max(0,w-1)):
     r,b=construct(U,entry,w,bits);need(strong(r),'entry-three strongness');key=key3(r,b)
     need(key not in banks[n],'entry-three canonical collision');banks[n].add(key);branches[0]+=1
  for m in range(1,n-3):
   k=n-m-1
   for U in rootone[m]:
    for i in range(k-2):
     j=k-3-i
     for smask in range(1<<(i+j-(1 if j else 0))):
      for tmask in ((0,1,3)if j==0 else(0,1,3,7)):
       r,b=tail_triangle(U,i,j,smask,tmask);need(strong(r),'triangle-tail strongness');key=key3(r,b)
       need(key not in banks[n],'triangle-tail canonical collision');banks[n].add(key);branches[1]+=1
  need(len(banks[n])==B3[n],('start-three inventory',n,len(banks[n]),B3[n]))
  hist={1:0,2:0,3:0}
  for r in banks[n]:hist[sum(bool(r[w[-1]]&1)for w in words(r,0,4))]+=1
  inventory.append({'order':n,'rooted_start_three':len(banks[n]),'cycle_histogram':hist,'adjacency_sha256':hashlib.sha256(json.dumps(sorted(banks[n]),separators=(',',':')).encode()).hexdigest()})
 bad=[0]*101
 for n in range(5,101):bad[n]=F[2*n-7]
 for n in range(4,11):
  actual=sum(not any(r[u].bit_count()==n-2 for u in range(1,n))for r in banks[n])
  need(actual==bad[n],('exceptional rooted start3 count',n,actual,bad[n]))
 c3stats=[]
 for n in range(5,11):
  marked=set();out=set()
  for m in range(5,n):
   for U in banks[m]:
    r,b=construct(U,(0,),n-m-1,0);marked.add(key3(r,b));out.add(whole_key(r))
  for m in range(1,n-3):
   k=n-m-1
   for U in rootone[m]:
    for i in range(k-2):
     r,b=tail_triangle(U,i,k-3-i,0,0);marked.add(key3(r,b));out.add(whole_key(r))
  expected=B3[n-1]+sum(B1[1:n-3]);need(len(out)==expected,'unmarked C3 with outdegree1')
  only=[r for r in out if not any(x.bit_count()==n-2 for x in r)]
  need(len(only)==bad[n-1],'one-sided C3 class count')
  allc=set(out)
  for r in out:
   rev=tuple(sum(1<<v for v in range(n)if r[v]>>u&1)for u in range(n));allc.add(whole_key(rev))
  need(len(allc)==expected+bad[n-1],'total C3 union')
  if c3stats:need(len(marked)==len(out)+c3stats[-1]['marked_outdegree_one'],'two-degree-one count bijection')
  c3stats.append({'order':n,'marked_outdegree_one':len(marked),'has_outdegree_one':len(out),'outdegree_one_but_no_indegree_one':len(only),'total_C3':len(allc)})
 A=[0]*101;C=[0]*101
 for n in range(5,101):A[n]=B3[n-1]+sum(B1[1:n-3]);C[n]=A[n]+bad[n-1]
 for n in range(101):
  lhs=C[n]-sum(c*C[n-j]for j,c in ((1,6),(2,-11),(3,6),(4,-1))if n>=j)
  need(lhs==({5:1,6:5,7:-17,8:11,9:-1}.get(n,0)),'complete C3 generating function')
 complete=[];labelcount=0
 for n in range(1,7):
  actual=set()
  for r in labelled(n):
   labelcount+=1
   if not strong(r):continue
   for root in range(n):
    if len(words(r,root,4))==3:actual.add(key3(r,root))
  need(actual==banks[n],'complete labelled root3 comparison');complete.append({'order':n,'classes':len(actual)})
 result={'status':'PASS','C3_physical_class_inventory':c3stats,'exceptional_start_three_counts':bad[1:11],'C3_generating_function':'z^5(1+5z-17z^2+11z^3-z^4)/(1-3z+z^2)^2','refined_C1_controls':clone_stats,'rooted_start_three_inventory':inventory,'physical_branch_counts':branches,'complete_labelled_tournaments':labelcount,'complete_labelled_inventory':complete,'recurrence_checked_through_order':100,'generating_function':'z^5(8-17z+8z^2)/(1-3z+z^2)^2'}
 data=json.dumps(result,sort_keys=True,indent=2)+'\n';(P/'start_three_certificate.json').write_text(data);print(data,end='')
if __name__=='__main__':main()
