#!/usr/bin/env python3
"""Exact Artin-D4 word controls and conditional low-degree actual passports.
The motivating geometric paths remain heuristic; no braid certificate here.
"""
from itertools import combinations,product
import hashlib,json
GATES=0

def check(value,label):
 global GATES
 GATES+=1
 if not value:raise RuntimeError(label)
def inv_word(w):return tuple(-i for i in w[::-1])
def word(*ws):
 out=[]
 for w in ws:
  for i in w:
   if out and out[-1]==-i:out.pop()
   else:out.append(i)
 return tuple(out)
def comm(a,b):return word(a,b,inv_word(a),inv_word(b))
def braid(a,b):return word(a,b,a,inv_word(b),inv_word(a),inv_word(b))
def conjugate(a,b):return word(a,b,inv_word(a))
def free_action(w,row=None):
 row=list(row or [(1,),(2,),(3,),(4,)])
 for j in w:
  k=abs(j)-1;a,b=row[k:k+2]
  row[k:k+2]=[conjugate(a,b),a] if j>0 else[b,conjugate(inv_word(b),a)]
 return row
WORDS={
 'cusp_plus':[-1,-2,1,1,1,2,1],
 'cusp_minus':[2,1,3,3,3,-1,-2],
 'cusp_two':[2,1,1,1,-2],
 'node0':[2,-1,2,2,1,-2],
 'node1':[2,1,2,3,3,-2,-1,-2],
 'node2':[2,3,2,2,-3,-2]}
SPLITS={
 'cusp_plus':([-1,-2],[1,1,1],0),
 'cusp_minus':([2,1],[3,3,3],2),
 'cusp_two':([2],[1,1,1],0),
 'node0':([2,-1],[2,2],1),
 'node1':([2,1,2],[3,3],2),
 'node2':([2,3],[2,2],1)}

def algebra_controls():
 a,b,c,d=[(i,) for i in range(1,5)];e=conjugate(b,c)
 expected=[(b,c),(b,d),(a,e),(conjugate(inv_word(e),a),b),(a,d),(e,conjugate(b,d))]
 for (name,(prefix,core,k)),pair in zip(SPLITS.items(),expected):
  check(WORDS[name]==prefix+core+list(inv_word(prefix)),'complete literal conjugated word')
  check(tuple(free_action(prefix)[k:k+2])==pair,'literal common-access core pair')
 R=braid(b,c)
 check(word(conjugate(e,b),inv_word(c))==R,'e b e inverse differs from c by exactly braid relator')
 check(word(conjugate(c,e),inv_word(b))==inv_word(R),'c e c inverse differs from b by inverse braid relator')
 check(conjugate(e,comm(conjugate(inv_word(e),a),b))==comm(a,conjugate(e,b)),
       'node0 conjugation is an exact free identity')
 check(conjugate(c,braid(a,e))==braid(conjugate(c,a),conjugate(c,e)),
       'third cusp conjugation is an exact free identity')

def mul(p,q):return tuple(p[q[i]] for i in range(len(p)))
def inverse(p):return tuple(p.index(i) for i in range(len(p)))
def conj(p,q):return mul(mul(p,q),inverse(p))
def commute(p,q):return mul(p,q)==mul(q,p)
def braided(p,q):return mul(mul(p,q),p)==mul(mul(q,p),q)
def support(p):return set(i for i,j in enumerate(p) if i!=j)
def fixed(p):return set(range(len(p)))-support(p)
def action(w,row):
 row=list(row)
 for j in w:
  k=abs(j)-1;a,b=row[k:k+2]
  row[k:k+2]=[conj(a,b),a] if j>0 else[b,conj(inverse(b),a)]
 return row

def permutation(n,*cycles):
 p=list(range(n))
 for cycle in cycles:
  for a,b in zip(cycle,cycle[1:]+cycle[:1]):p[a]=b
 return tuple(p)
def group(gens):
 n=len(gens[0]);identity=tuple(range(n));found={identity};todo=[identity]
 for p in todo:
  for g in gens:
   q=mul(p,g)
   if q not in found:found.add(q);todo.append(q)
 return found

def transposition_controls():
 counts=[]
 for n in range(2,8):
  ts=[permutation(n,list(c)) for c in combinations(range(n),2)]
  b=ts[0];leaves=[t for t in ts if braided(b,t)];valid=0;max_support=0
  for a,c,d in product(leaves,repeat=3):
   if not(commute(a,c) and commute(a,d) and commute(c,d)):continue
   valid+=1;size=len(support(a)|support(b)|support(c)|support(d))
   check(size<=4,'every transposition D4 image moves at most four labels')
   if a==b or c==b or d==b:check(a==b==c==d,'central leaf equality forces all leaves equal')
   else:check(len({a,c,d})<=2,'at most two distinct noncentral commuting leaves')
   max_support=max(max_support,size)
  counts.append([n,valid,max_support])
 return counts

def scalar_controls():
 rows=[]
 for d in [5,6,7]:
  for N in [3,4,5]:
   aa=set();nr=0
   for a in range(1,d):
    delta=d-a;q=d-2*a
    for ns in product(range(a+1),repeat=3):
     if any(2*n<3*a-d for n in ns):continue
     W=2*a+1-sum(ns)
     parity=any(2*n==3*a-d for n in ns)
     allowed=[k for k in range(max(0,q),delta+1) if not parity or k%2==0]
     sums={0}
     for _ in range(N):sums={s+k for s in sums for k in allowed}
     if W not in sums:continue
     nr+=1;aa.add(a)
     check(a in({3} if d==5 else {3,4} if d==6 else {4}),'conditional scalar and saturated-parity a bounds')
   check(aa==({3} if d==5 else {3,4} if d==6 else {4}),'complete named scalar universe')
   rows.append([d,N,sorted(aa),nr])
 return rows

def three_cycle_controls():
 rows=[]
 for d,a in [(6,3),(7,4)]:
  ps=[permutation(d,list(c)) for c in combinations(range(d),3)]
  ps+= [inverse(p) for p in ps]
  b=ps[0]
  for c in ps:
   if commute(b,c):check(len(support(b)&support(c)) in [0,3],'commuting three-cycle intersection')
   n=len(fixed(b)&fixed(c))
   check(n==d-6+len(support(b)&support(c)),'actual full-fixed cusp intersection')
   if braided(b,c):check(2*n>=3*a-d,'local injection on full-fixed actual sets')
  for c,dd in product(ps,repeat=2):
   if support(c)&support(dd):continue
   lower=(3*a-d+1)//2
   if len(fixed(b)&fixed(c))<lower or len(fixed(b)&fixed(dd))<lower:continue
   check(False,'one three-cycle cannot meet two disjoint three-sets twice')
  rows.append([d,a,len(ps)])
 return rows

def named_control(n,gens,expected_order,expected_nc,expected_omega,expected_chi):
 a,b,c,d=gens;e=conj(b,c);f=conj(inverse(e),a)
 for name,w in WORDS.items():check(action(w,gens)==list(gens),'full literal word on named control '+name)
 GG=group(gens);check(len(GG)==expected_order,'named image order')
 check({g[0] for g in GG}==set(range(n)),'named action transitive')
 kept=len(fixed(b));check(all(len(fixed(g))==kept for g in gens),'same fixed-letter number')
 cusp_pairs=[(b,c),(b,d),(a,e)];node_pairs=[(f,b),(a,d),(c,d)]
 nc=[];ww=[]
 for sigma,tau in cusp_pairs:
  A=fixed(sigma);B=fixed(tau)
  check({mul(sigma,tau)[i] for i in A}==B,'actual full-fixed re-access on named cusp')
  nc.append(len(A&B))
 for sigma,tau in node_pairs:
  check(commute(sigma,tau),'named node commute')
  ww.append(len(support(sigma)&support(tau)))
  check(len(fixed(sigma)&fixed(tau))==2*kept-n+ww[-1],'actual node count')
 chi=-2*kept+sum(nc)+sum(ww)
 check(nc==expected_nc and ww==expected_omega and chi==expected_chi,'named exact Euler ledger')
 return dict(degree=n,kept=kept,group_order=len(GG),cusps=nc,nodes=ww,chi=chi)

def main():
 algebra_controls();trans=transposition_controls();scalars=scalar_controls();three=three_cycle_controls()
 p=lambda *cycles:permutation(4,*cycles)
 s4=named_control(4,[p([2,3]),p([1,2]),p([0,1]),p([2,3])],24,[1,1,1],[0,2,0],1)
 # Labels 0..3 are +1..+4 and 4..7 are -1..-4.
 p=lambda *cycles:permutation(8,*cycles)
 wd4=named_control(8,[p([0,1],[4,5]),p([1,2],[5,6]),p([2,3],[6,7]),p([2,7],[6,3])],192,[2,2,2],[0,0,4],2)
 data=dict(transposition_universe=trans,scalar_universe=scalars,three_cycle_controls=three,S4_actual_ledger_control=s4,W_D4_stopping_control=wd4)
 print('Exact abstract Artin-D4 and conditional mapping-degree5/6/7 passport gates PASS')
 print('GEOMETRIC WORDS REMAIN HEURISTIC: no actual three-cusp complement theorem asserted')
 print(json.dumps(data,sort_keys=True,indent=2))
 print('Always-active gates:',GATES)
 print('Semantic SHA256:',hashlib.sha256(json.dumps(data,sort_keys=True,separators=(',',':')).encode()).hexdigest())
if __name__=='__main__':main()
