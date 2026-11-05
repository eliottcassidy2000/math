#!/usr/bin/env python3
"""Actual bridge-C2 class checks, with pinned start-two graph helpers."""
import json,hashlib,itertools as it
from pathlib import Path
from verify_start_two import construct,relabel,cycle_key,rooted_two_key,words,strong,labelled,need
P=Path(__file__).resolve().parent
need(hashlib.sha256((P/'verify_start_two.py').read_bytes()).hexdigest()=='ca788fa6489cbca7bed8f9f9a84a553ff9f5de7502833544f9b51b9f059142ae','dependency hash')

def indegree_one(r):return {v for v in range(len(r)) if sum(x>>v&1 for x in r)==1}
def outdegree_one(r):return {v for v,x in enumerate(r) if x.bit_count()==1}
def induced(r,vs):return relabel(r,vs)
def cycle_words(r,cap=3):
 n=len(r);full=(1<<n)-1;out=[]
 if n<3:return out
 def visit(w,mask):
  if mask==full:
   if r[w[-1]]&1:out.append(tuple(w))
   return
  todo=r[w[-1]]&~mask
  while todo and len(out)<cap:
   bit=todo&-todo;todo-=bit;visit(w+[bit.bit_length()-1],mask|bit)
 visit([0],1);return out

def unrooted_key(r,cs):
 return min(relabel(r,w[i:]+w[:i]) for w in cs for i in range(len(r)))

def main():
 one={1:{(0,)},2:set()};pairs={};pairchecks=0
 for n in range(3,10):
  one[n]=set()
  for m in range(1,n):
   w=n-m-1
   if m==1 and w==0:continue
   for U in one[m]:
    for bits in range(1<<max(0,w-1)):
     r,b=construct(U,(0,),w,bits);one[n].add(relabel(r,(b,)+tuple(range(b))))
  unroot={cycle_key(r) for r in one[n]};pairs[n]=[]
  for U in unroot:
   S=[v for v in range(n) if len(words(U,v))==1];ins=indegree_one(U)
   if n>=4:need(ins-set(S),'C1 lacks indegree-one vertex outside start-one set')
   choices=list(it.combinations(S,2))
   for pair in choices:need(ins-set(pair),'pair covers all indegree-one vertices');pairchecks+=1
   if n==3:choices=choices[:1]
   pairs[n].extend((U,pair) for pair in choices)
 B=[0]*101;D=[0]*101;F=[0,1]
 for _ in range(205):F.append(F[-1]+F[-2])
 D[3]=1
 for n in range(4,101):D[n]=F[2*n-8]+3*F[2*n-7]
 for n in range(4,101):B[n]=sum((1 if n-m==1 else 2**(n-m-2))*(B[m]+D[m]) for m in range(3,n))
 two={n:set() for n in range(1,11)};bridge={n:set() for n in range(1,11)};c=[0]*11;stats=[];indchecks=autchecks=inversechecks=0
 for n in range(4,11):
  for m in range(3,n):
   w=n-m-1
   for U,entry in [(r,(0,)) for r in two[m]]+pairs[m]:
    for bits in range(1<<max(0,w-1)):
     r,b=construct(U,entry,w,bits);key=rooted_two_key(r,b)
     need(key not in two[n],'rooted decomposition collision');two[n].add(key)
  need(len(two[n])==B[n],'rooted bank size')
  for r in two[n]:
   need(indegree_one(r)-{0},'start-two indegree lemma');indchecks+=1
   ws=words(r,0);cs=[w for w in ws if r[w[-1]]&1]
   if len(cs)==2:
    need(outdegree_one(r) and indegree_one(r),'degree-one equivalences')
    bridge[n].add(unrooted_key(r,cs))
   if r[0].bit_count()==1:c[n]+=1
  need(c[n]==sum(B[:n]),'marked outdegree-one count')
  need(len(bridge[n])==B[n-1],'unrooted bridge count')
  double=0;inverse=set()
  for r in bridge[n]:
   cs=cycle_words(r);need(len(cs)==2,'unrooted C2')
   automorphisms=sum(relabel(r,w[j:]+w[:j])==r for w in cs for j in range(n));need(automorphisms==1,'bridge C2 rigidity');autchecks+=1
   lows=outdegree_one(r);need(len(lows) in(1,2),'low degree alternatives')
   if len(lows)==2:
    double+=1;u,v=sorted(lows)
    if not r[u]>>v&1:u,v=v,u
    vs=[x for x in range(n) if x!=u];S=induced(r,vs);root=vs.index(v)
    need(strong(S),'inverse graph not strong');key=rooted_two_key(S,root)
    need(key in two[n-1] and key[0].bit_count()==1,'inverse marked graph');inverse.add(key);inversechecks+=1
  target={r for r in two[n-1] if r[0].bit_count()==1}
  need(inverse==target,'two-low-degree inverse bijection')
  need(double==c[n-1],'double low count')
  stats.append(dict(order=n,rooted_start_two=len(two[n]),marked_outdegree_one=c[n],unrooted_bridge_C2=len(bridge[n]),two_outdegree_one=double))
 complete=[];total=two_strong_roots=0
 for n in range(1,7):
  actual_bridge=set();atoms=set()
  for r in labelled(n):
   total+=1
   if n<3 or not strong(r):continue
   is2=all(strong(induced(r,[u for u in range(n) if u!=v])) for v in range(n))
   if is2:
    for v in range(n):need(len(words(r,v))>=3,'2-strong root lower bound');two_strong_roots+=1
   cs=cycle_words(r)
   if len(cs)!=2:continue
   need(bool(outdegree_one(r))==bool(indegree_one(r)),'complete degree equivalence')
   key=unrooted_key(r,cs)
   if outdegree_one(r):actual_bridge.add(key)
   else:
    atoms.add(key)
    for u in range(n):
     for v in range(n):
      if r[u]>>v&1:
       x=list(r);x[u]^=1<<v;need(strong(x),'degree-free C2 unexpectedly has bridge')
  need(actual_bridge==bridge[n],'complete labelled bridge C2 inventory')
  need(len(atoms)==(1 if n==5 else 0),'small residual atoms')
  complete.append(dict(order=n,bridge_classes=len(actual_bridge),bridgefree_classes=len(atoms)))
 need(5**41<=2**96<5**42,'residual order arithmetic')
 result=dict(status='PASS',C1_pair_indegree_controls=pairchecks,start_two_indegree_controls=indchecks,bridge_C2_rigidity_checks=autchecks,two_vertex_inverse_checks=inversechecks,physical_inventory=stats,complete_labelled_tournaments=total,complete_labelled_max_order=6,complete_unrooted_C2_inventory=complete,two_strong_root_lower_bound_checks=two_strong_roots,formula_checked_through_order=100,exact_total_C2_formula_from_order=74)
 data=json.dumps(result,sort_keys=True,indent=2)+'\n';(P/'bridge_two_certificate.json').write_text(data);print(data,end='')
if __name__=='__main__':main()
