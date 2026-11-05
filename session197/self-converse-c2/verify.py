import hashlib,json
from pathlib import Path
from converse import *
from core import tournament
from bilinear import source_bank
P=Path(__file__).parent
EXPECTED={'BILINEAR_C2.md':'5eab9deb2871819d8a1b5684c4987f2f4da807adc55063a30675d531c5ad3813','C2_COMPLETE.md':'4cc0bf11efc435d7a2746daf6d8f1a58c49c94e57e608911b24464facdd1057d','FIBONACCI_CLASS.md':'d95a133e360b96589009ebd8a1c8449525c4a2847a88b33a4b3f2b839739e63e'}
for f,h in EXPECTED.items():require(hashlib.sha256((P/'dependencies'/f).read_bytes()).hexdigest()==h,'dependency '+f)

def hamilton(A):
 n=len(A);N=1<<n;D=[[0]*n for _ in range(N)]
 for v in range(n):D[1<<v][v]=1
 for M in range(1,N):
  for v in range(n):
   if not D[M][v]:continue
   rem=A[v]&((N-1)^M)
   while rem:
    z=rem&-rem;rem-=z;D[M|z][z.bit_length()-1]+=D[M][v]
 return sum(D[-1])

def rooted_key(A,c):
 q=inventories(A);require(q['s'][c]==1,'unique rooted path');M=(1<<len(A))-1;word=[c]
 while len(word)<len(A):
  M^=1<<word[-1]
  nxt=[x for x in range(len(A))if M>>x&1 and A[word[-1]]>>x&1 and q['start'][M][x]]
  require(len(nxt)==1,'unique next vertex');word.append(nxt[0])
 return relabel(A,word)

src={};source_counts={};seenprofiles={};hostile=None
for n in range(2,8):
 keys=set();source_counts[n]=len(sources(n))
 for A,c in sources(n):
  key=rooted_key(A,c);require(key not in keys,'duplicate canonical source');keys.add(key)
  v=profile(A,c);require(v[0]%2==1 and(v[2]-v[1])%2==0,'source parity');src[A,c]=v
  if v in seenprofiles and rooted_key(*seenprofiles[v])!=key and hostile is None:hostile=(seenprofiles[v],(A,c),v)
  seenprofiles[v]=(A,c)
require(list(source_counts.values())==[1,2,5,13,34,89],'source GF coefficients')
classes={n:set()for n in range(5,11)};self_counts={n:0 for n in classes};swaps=0;direct=0;diagonal=0
for n in range(2,8):
 for m in range(2,8):
  for kind in('N','M'):
   order=n+m+(1 if kind=='N'else 2)
   if order>10:continue
   for V,c in sources(n):
    for W,a in sources(m):
     T,_=construct(kind,V,c,W,a);U,_=construct(kind,W,a,V,c)
     p=swap_map(kind,n,m);require(adjacency_map(converse(T),U,p),'literal converse swap');swaps+=1
     C=cycle_words(T);require(len(C)==2,'C2 physical graph')
     key=canonical(T,C);require(key not in classes[order],'duplicate ordered presentation');classes[order].add(key)
     v=src[V,c];w=src[W,a];h=hamilton(T);require(h==value(kind,v,w),'literal scalar count');direct+=1
     maps=anti_maps(T,C);same=(n==m and rooted_key(V,c)==rooted_key(W,a))
     require(bool(maps)==same,'self-converse iff source isomorphism')
     if maps:
      require(len(maps)==1 and all(maps[0][maps[0][x]]==x for x in range(order)),'unique involution')
      require(sum(maps[0][x]==x for x in range(order))==(kind=='N'),'fixed central vertices')
      require(h%4==1,'diagonal mod4');self_counts[order]+=1;diagonal+=1
R5=[6,12,24,17,3];require(hamilton(R5)==15 and len(cycle_words(R5))==2 and len(anti_maps(R5))==5,'regular exception');self_counts[5]+=1
require([len(classes[n])for n in classes]==[1,5,18,60,191,589],'full graph class coefficients')
require(list(self_counts.values())==[2,1,2,2,5,5],'self-converse GF coefficients')
require(hostile is not None,'quadruple collision')
(V,c),(W,a),v=hostile;hostile_profile=v
T,_=construct('N',V,c,W,a);U,_=construct('N',V,c,V,c)
h=hamilton(T);require(h==hamilton(U)==value('N',v,v),'same-profile graph counts')
require(not anti_maps(T)and len(anti_maps(U))==1,'equal response is not source isomorphism')
# An exact bounded scalar inverse, using only actual source states and the proved diagonal cutoff.
bank=source_bank(100);image={15}
for v in bank:
 if v==(1,1,1,1):continue
 for kind in('N','M'):
  h=value(kind,v,v)
  if h<=10000:image.add(h)
require(all(h==15 or h%4==1 for h in image),'bounded diagonal residue')
for target in range(1,501,2):
 result=diagonal_inverse(target);require(bool(result)==(target in image),'complete inverse control')
 if result:require(hamilton(result['rows'])==target and anti_maps(result['rows']),'physical inverse witness')
receipt=dict(status='PASS',dependencies=EXPECTED,source_classes=source_counts,source_parity_controls=len(src),literal_swaps=swaps,full_C2_classes={n:len(classes[n])for n in classes},literal_graph_counts=direct,diagonal_graphs=diagonal,self_converse_classes=self_counts,regular_anti_isomorphisms=5,hostile=dict(source_left=list(V),root_left=c,source_right=list(W),root_right=a,quadruple=hostile_profile,target_rows=T,diagonal_rows=U,H=value('N',hostile_profile,hostile_profile)),bounded_self_converse_support_10000=sorted(image),exact_inverse_targets=250)
(P/'VALIDATION.json').write_text(json.dumps(receipt,indent=2,sort_keys=True)+'\n');print(json.dumps({k:v for k,v in receipt.items()if k not in('hostile','bounded_self_converse_support_10000')},sort_keys=True))
