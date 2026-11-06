#!/usr/bin/env python3
"""Exact -O-safe controls; universal justification is the one-port basis lemma."""
import hashlib,json,random
from fractions import Fraction
from pathlib import Path
import importlib.util
from counting import profiles,nine,labelled,gate,literal,J,relabel,check_tournament
from compiler import move,base,inventory,compile_target
P=Path(__file__).resolve().parent
MATRICES={1:{1:[[1,-1,3,-2],[0,0,1,-1],[2,-1,4,-1],[1,-1,1,0]],2:[[4,-6,10,-4],[1,-2,2,-1],[3,-2,6,-1],[2,-2,2,0]]},2:{1:[[3,-2,6,-2],[1,0,1,-1],[3,-1,8,-1],[1,0,2,0]],2:[[12,-12,28,-8],[3,-2,4,-2],[7,-4,18,-2],[4,-2,6,0]]}}
def need(x,msg):
 if not x:raise RuntimeError(msg)
def solve(X,Y):
 a=[[Fraction(x)for x in u+v]for u,v in zip(X,Y)]
 for j in range(4):
  pivot=next(i for i in range(j,4)if a[i][j]);a[j],a[pivot]=a[pivot],a[j]
  c=a[j][j];a[j]=[x/c for x in a[j]]
  for i in range(4):
   if i==j:continue
   c=a[i][j];a[i]=[x-c*y for x,y in zip(a[i],a[j])]
 return [[a[j][i+4]for j in range(4)]for i in range(4)]
def op_graph(a,op):
 if op=='K1':return gate(a,1)
 if op=='K2':return gate(a,2)
 n=len(a)
 if op=='O':return tuple([x|(1<<n)for x in a]+[0])
 if op=='Z':return tuple([a[0]]+[x|(1<<n)for x in a[1:]]+[1])
 return J(a)
def check_one(a):
 check_tournament(a);p=nine(a);n=len(a)
 for op in('O','Z','J','K1','K2'):
  out=op_graph(a,op);check_tournament(out);actual=nine(out);pred=move(p,op)
  if n==1 and op=='Z':pred=pred[:1]+(0,)+pred[2:]
  if n==1 and op=='K1':pred=pred[:1]+(1,)+pred[2:]
  need(actual==pred,('transfer',a,op,p,actual,pred))
  if n>=2 and p[1]>0:
   h,s,e,L,R,d,C,b,T=p
   need(R>=d+s and d>=e and h>=s and C>=2*L and b>=L and T>=b,'cone')
   if op=='O':need(base(actual)==base(p),'O base')
   else:need(base(actual)>=base(p)+2,'strict growth')
   if op=='K1':need(base(actual)>=base(p)+2*s and actual[3]>=2**(n-1),'K1 cap')
   if op=='K2':need(actual[3]>=2**(n+1),'K2 cap')
   if op in('Z','J'):need(base(actual)>=2**(n-1),'Z/J cap')
 return 5

def canonical_inventory(N=8):
 spec=importlib.util.spec_from_file_location('inherited_start_two',P/'dependencies/verify_start_two.py');old=importlib.util.module_from_spec(spec);spec.loader.exec_module(old)
 rootone={1:{(0,)},2:set()};pairs={};two={n:set()for n in range(1,N+1)}
 for n in range(3,N+1):
  bank=set()
  for m in range(1,n):
   w=n-m-1
   for U in rootone[m]:
    for bits in range(1<<max(0,w-1)):
     a,r=old.construct(U,(0,),w,bits);bank.add(old.relabel(a,(r,)+tuple(range(r))))
  rootone[n]=bank;unroot={old.cycle_key(a)for a in bank};marked=[]
  for a in unroot:
   S=[i for i in range(n)if len(old.words(a,i,3))==1]
   import itertools
   marks=list(itertools.combinations(S,2))
   if n==3:marks=marks[:1]
   marked.extend((a,p)for p in marks)
  pairs[n]=marked
  for m in range(3,n):
   w=n-m-1
   for U,entry in [(U,(0,))for U in two[m]]+pairs[m]:
    for bits in range(1<<max(0,w-1)):
     a,r=old.construct(U,entry,w,bits);two[n].add(old.rooted_two_key(a,r))
 # Independently generate the new finite word grammar on literal graphs.
 banks={n:{1:set(),2:set()}for n in range(2,N+1)};banks[2][1].add((2,0))
 for n in range(2,N+1):
  for s in(1,2):
   for a in banks[n][s]:
    for op in('O','Z','J','K1','K2')if s==1 else('O','Z','J'):
     out=op_graph(a,op);nn=len(out)
     if nn>N:continue
     ws=old.words(out,0,3);ss=s*(2 if op.startswith('K')else 1)
     need(len(ws)==ss,'word start')
     key=min(old.relabel(out,w)for w in ws);banks[nn][ss].add(key)
 stats=[]
 for n in range(2,N+1):
  strongs={a for a in banks[n][2]if old.strong(a)}
  need(strongs==two[n],('canonical source discrepancy',n,len(strongs),len(two[n])))
  stats.append(dict(order=n,canonical_strong_start_two=len(two[n]),gate_strong_start_two=len(strongs),gate_weak_or_strong_start_two=len(banks[n][2])))
 # All B<=31 have a representative through order 7, by the proved cap.
 actual={base(nine(a))for n in range(2,8)for a in two[n]if base(nine(a))<=31}
 info,w=inventory(31)
 need(info['complete']and set(w)==actual,'independent low-target support')
 return stats,sorted(actual)

def main():
 basis=[(0,),(2,0),(4,1,2),(6,0,2)];reconstructed={}
 for k in(1,2):
  X=[list(profiles(a,2)[0][k-1])for a in basis]
  for g in(1,2):
   Y=[list(profiles(gate(a,g),2)[0][k-1])for a in basis]
   M=solve(X,Y);need(M==MATRICES[g][k],'basis matrix mismatch')
   reconstructed[f'K{g}_level{k}']=[[int(x)for x in row]for row in M]
 literal_checks=0;inputs=0
 for n in range(1,6):
  for a in labelled(n):inputs+=1;literal_checks+=check_one(a)
 random.seed(203)
 for n in(6,7):
  for _ in range(10):
   a=[0]*n
   for i in range(n):
    for j in range(i+1,n):
     if random.getrandbits(1):a[i]|=1<<j
     else:a[j]|=1<<i
   inputs+=1;literal_checks+=check_one(tuple(a))
 stats,small=canonical_inventory()
 good=compile_target(137);need(good['status']=='PHYSICAL_START_TWO_WITNESS','137 positive')
 need(nine(tuple(good['rows']))==tuple(good['profile'])and good['profile'][1]==2 and base(good['profile'])==137,'137 literal')
 bad=compile_target(139);need(bad['status']=='ABSENT_START_TWO_PROFILE'and bad['complete'],'139 scope')
 incomplete=compile_target(139,1);need(incomplete['status']=='INCOMPLETE'and not incomplete['complete'],'resource guard')
 for target in(0,2,-1):
  try:compile_target(target)
  except ValueError:pass
  else:raise RuntimeError('accepted invalid target')
 # Minimum target and identity/singleton boundaries.
 need(compile_target(1)['status']=='ABSENT_START_TWO_PROFILE','B=1 control')
 need(nine(gate((0,),1))[1]==1 and nine(gate((0,),2))[1]==2,'singleton starts')
 need(nine(gate((0,),2))==nine(gate((2,0),1)),'singleton K2/edge K1 boundary')
 result=dict(status='PASS',rooted_inputs=inputs,literal_children_counted=literal_checks,complete_labelled_max_order=5,random_extra_orders=[6,7],basis_matrices=reconstructed,canonical_inventory=stats,complete_start_two_bases_through31=small,positive_B137_order=good['order'],negative_B139_states=bad['states'],resource_limit_status=incomplete['status'],dependencies={p.name:hashlib.sha256(p.read_bytes()).hexdigest()for p in sorted((P/'dependencies').glob('*'))if p.is_file()})
 print(json.dumps(result,sort_keys=True,indent=2))
if __name__=='__main__':main()
