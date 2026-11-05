from pathlib import Path
from itertools import combinations
import json
from core import require,inventories,cycles
from pair_families import construct as family
from structure import *
P=Path(__file__).parent

def relabel(T,order):return tuple(induced(T,list(order)))
def words(T,r,cap=3):
 full=(1<<len(T))-1;found=[]
 def visit(p,mask):
  if mask==full:found.append(tuple(p));return
  bits=T[p[-1]]&~mask
  while bits and len(found)<cap:
   bit=bits&-bits;bits-=bit;visit(p+[bit.bit_length()-1],mask|bit)
 visit([r],1<<r);return found

def canonical(T,cs=None):
 if cs is None:cs=[p for p in words(T,0,10**12)if T[p[-1]]&1]
 return min(relabel(T,p[k:]+p[:k])for p in cs for k in range(len(T)))
def rootkey(T,r):
 ws=words(T,r);require(len(ws)==2,'two actual paths')
 return min(relabel(T,p)for p in ws)
def one_key(T):
 n=len(T);return min(relabel(T,tuple(range(k,n))+tuple(range(k)))for k in range(n))
def attach(U,entries,w,bits):
 m=len(U);n=m+w+1;r=list(U)+[0]*(w+1);b=n-1
 for x in range(m):
  for y in range(m,b):r[x]|=1<<y
  r[b if x in entries else x]|=1<<(x if x in entries else b)
 for i in range(w):
  x=m+i
  for y in range(x+1,b):r[x]|=1<<y
  r[b if bits>>i&1 else x]|=1<<(x if bits>>i&1 else b)
 return r,b

N=10;one={1:{(0,)},2:set()};pairs={};two={n:set()for n in range(1,N+1)};bank={n:set()for n in range(1,N+1)}
for n in range(3,N+1):
 one[n]=set()
 for m in range(1,n):
  w=n-m-1
  if m==1 and w==0:continue
  for U in one[m]:
   for bits in range(1<<max(0,w-1)):
    T,r=attach(U,(0,),w,bits);one[n].add(relabel(T,(r,)+tuple(range(r))))
 pairs[n]=[]
 for U in {one_key(T)for T in one[n]}:
  starts=[r for r in range(n)if len(words(U,r))==1];choices=list(combinations(starts,2))
  if n==3:choices=choices[:1]
  pairs[n].extend((U,pair)for pair in choices)
 for m in range(3,n):
  w=n-m-1
  for U,entry in [(T,(0,))for T in two[m]]+pairs[m]:
   for bits in range(1<<max(0,w-1)):
    T,r=attach(U,entry,w,bits);two[n].add(rootkey(T,r))
 for T in two[n]:
  cs=[p for p in words(T,0)if T[p[-1]]&1]
  if len(cs)==2:bank[n].add(canonical(T,cs))
R5=(6,12,24,17,3);bank[5].add(canonical(R5))
image={n:set()for n in range(1,N+1)};present={n:0 for n in range(1,N+1)};bykind={n:{k:0 for k in 'NMPQ'}for n in range(1,N+1)};redundancy=0
for m in range(1,N+1):
 for c in range(1,m+1):
  for U in one.get(c,()):
   V=list(U)
   for _ in range(m-c):
    n=len(V);V=[x|(1<<n)for x in V]+[0]
   if m<=5:
    for v in range(1,17):redundant_P(V,0,v);redundancy+=1
   for kind,gate in [('N',2),('M',3),('P',4),('Q',3)]:
    if m==1 and kind!='P':continue
    for size in range(N-m-gate+1):
     if kind=='N'and size==0:continue
     choices=[1]if size==0 else range(2**(size-1)+1,2**size+1)
     for v in choices:
      T,r=family(kind,V,0,v);n=len(T)
      require(n==m+gate+size,'presentation size')
      cs=[p for p in words(T,0,10**12)if T[p[-1]]&1];require(len(cs)==2,'Gamma C2')
      image[n].add(canonical(T,cs));present[n]+=1;bykind[n][kind]+=1
      rec=compile_root(T,r);require(rec is not None,'designated root compiler')
      if kind!='P':
       sourcekey=relabel(rec['source'],words(rec['source'],rec['source_root'],2)[0])
       require(rec['family']==kind and rec['parameter']==v and sourcekey==tuple(V),'canonical presentation inverse')
checks=roots=0;stats=[]
for n in range(5,N+1):
 detected=set();outsideboth=qualifying=0
 for T in bank[n]:
  result=compile_graph(T);checks+=1
  if result is not None:detected.add(T)
  elif compile_graph(converse(T))is None:outsideboth+=1
  for rec in root_records(T):
   result=compile_root(T,rec['root'])
   require((result is not None)==(rec['initial_cycles']==1),'every qualifying root');roots+=1;qualifying+=result is not None
 require(detected==image[n],'complete graph criterion vs complete physical image')
 require(qualifying==sum(bykind[n][k]for k in 'NMQ'),'exact rooted presentation count')
 stats.append(dict(order=n,all_C2=len(bank[n]),presentations=present[n],Gamma_graphs=len(image[n]),outside_both_orientations=outsideboth,canonical_rooted_presentations=qualifying))
require([r['all_C2']for r in stats]==[2,5,18,60,191,589],'complete C2 counts')
require([r['Gamma_graphs']for r in stats]==[1,5,18,57,169,482],'image counts')
hostiles=[]
for h,T,param in[(77,[2,4,9,19,103,79,143,63],8),(289,[2,4,9,83,39,79,151,319,639,255],26)]:
 I=inventories(T);require(I['H']==h and cycles(T)==2 and compile_graph(T)is None,'explicit hostile')
 fullpaths=[p for r in range(len(T))for p in words(T,r,10**12)]
 require(len(fullpaths)==h and sum(T[p[-1]]&1 for p in fullpaths if p[0]==0)==2,'literal H/C hostile')
 op=converse(T)
 if h==289:require(compile_graph(op)is None,'two-sided hostile')
 A,r=family('P',[0],0,param)
 require(len(A)==len(T)and inventories(A)['H']==h and cycles(A)==2,'same-order scalar rebase')
 hostiles.append(dict(H=h,rows=T,root_tests=root_records(T),converse_root_tests=root_records(op),converse_compilation=compile_graph(op),repair=dict(family='P',source=[0],source_root=0,parameter=param,rows=A,root=r)))
receipt=dict(status='PASS',complete_C2_graphs=checks,start_two_root_checks=roots,P_redundancy_checks=redundancy,inventory=stats,hostiles=hostiles)
text=json.dumps(receipt,sort_keys=True,indent=2)+'\n';(P/'VALIDATION.json').write_text(text);print(text,end='')
