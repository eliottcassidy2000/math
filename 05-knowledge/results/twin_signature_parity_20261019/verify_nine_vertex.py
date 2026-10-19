#!/usr/bin/env python3
"""Exact finite controls for the hand-proved n8 padding/n9 criterion."""
from verify_twin_signature import *
from math import comb
import pathlib,json

def split(a,u,v):
 parts=[]
 for start,banned in ((u,v),(v,u)):
  vv=[start];seen={start,banned}
  for x in vv:
   for y in a[x]:
    if y not in seen:seen.add(y);vv.append(y)
  ix={x:i for i,x in enumerate(vv)}
  parts.append([[ix[y] for y in a[x] if y in ix] for x in vv])
 return sorted(parts,key=len)

def criterion(a):
 if aut(a)==1:return False
 for u in range(9):
  for v in a[u]:
   if u<v:
    parts=split(a,u,v)
    if len(parts[0]) in (1,2) and aut(parts[1])==1:return False
 return True

def extensions(n,out):
 pred=[sum(1<<u for u in range(n) if out[u]>>v&1) for v in range(n)]
 dp=[0]*(1<<n);dp[0]=1
 for mask in range(1<<n):
  if not dp[mask]:continue
  for v in range(n):
   if not(mask>>v&1) and not(pred[v]&~mask):dp[mask|(1<<v)]+=dp[mask]
 return dp[-1]

def main():
 ts=trees(9);good=[];bad=[]
 for e in ts[9]:
  (good if criterion(adj_of(e,9)) else bad).append(e)
 require(len(good)==36 and len(bad)==11,'criterion split')
 require(sum(not hereditary(adj_of(e,9)) for e in good)==12,'nonhereditary count')
 # Classify the twelve new good shapes using their now-proved rigidity.
 new_good_summary={'automorphism_forced_even':0,'additional_two_cherry':0,'additional_all_even':0,'has_odd_orientations':0}
 for e in good:
  a=adj_of(e,9)
  if hereditary(a):continue
  parities=set()
  for mask in range(256):
   _,out=orient_data(e,9,mask);parities.add(extensions(9,out)%2)
  if 1 in parities:new_good_summary['has_odd_orientations']+=1
  elif afree(a)==0:new_good_summary['automorphism_forced_even']+=1
  elif sum(sum(len(a[w])==1 for w in a[v])>=2 for v in range(9))>=2:new_good_summary['additional_two_cherry']+=1
  else:new_good_summary['additional_all_even']+=1
 require(new_good_summary=={'automorphism_forced_even':3,'additional_two_cherry':3,'additional_all_even':4,'has_odd_orientations':2},'new good decomposition')
 # Every two-edge deletion has zero parity by a structural union certificate.
 from itertools import combinations
 from math import factorial
 second_cuts=0
 for e in ts[9]:
  for i,j in combinations(range(8),2):
   kept=[x for k,x in enumerate(e) if k not in(i,j)];aa=adj_of(kept,9);seen=set();parts=[]
   for v in range(9):
    if v in seen:continue
    vv=[v];seen.add(v)
    for x in vv:
     for y in aa[x]:
      if y not in seen:seen.add(y);vv.append(y)
    ix={x:k for k,x in enumerate(vv)};parts.append([[ix[y] for y in aa[x]] for x in vv])
   sizes=sorted(map(len,parts));require(len(parts)==3,'second cut components')
   if sizes[-1]>=7:require(sizes==[1,1,7],'large second-cut branch')
   else:
    require(all(hereditary(x) for x in parts),'small second-cut component')
    coefficient=factorial(9)
    for size in sizes:coefficient//=factorial(size)
    require(coefficient%2==0,'odd three-component shuffle')
   second_cuts+=1
 require(second_cuts==1316,'second-cut count')
 # The complete parity derivative is supported exactly on critical leaf edges.
 D=((0,1),(0,2),(0,3),(1,4),(1,5),(2,6),(4,7))
 P8=tuple((i,i+1) for i in range(7))
 odd_padded_shapes={canon(adj_of(D,8)),canon(adj_of(P8,8))}
 selector_masks=0
 for e in good:
  a=adj_of(e,9);critical=0
  for i,(u,v) in enumerate(e):
   pp=split(a,u,v)
   if len(pp[0])==1 and canon(pp[1]) in odd_padded_shapes:critical|=1<<i
  _,out=orient_data(e,9,0);base=extensions(9,out)%2
  for mask in range(256):
   _,out=orient_data(e,9,mask)
   require(extensions(9,out)%2==(base^((mask&critical).bit_count()%2)),'critical-leaf selector')
   selector_masks+=1
 # Eight-tree cut universe: all 23 shapes, seven edges each.
 cuts8=0
 for e in ts[8]:
  a=adj_of(e,8)
  for u,v in e:
   p=split(a,u,v);s=len(p[0]);cuts8+=1
   require(1<=s<=4,'cut size')
   if s>1:
    require(all(hereditary(x) for x in p),'small cut not good')
    require((9*comb(8,s))%2==0,'padding shuffle not even')
 require(cuts8==161,'eight-tree cut count')
 # All eleven rejected classes have explicit finite nonconstancy witnesses.
 wdata=json.loads(pathlib.Path(__file__).with_name('nine_vertex_bad_witnesses.json').read_text())
 expected={canon(adj_of(e,9)) for e in bad};observed=set()
 for w in wdata['witnesses']:
  e=tuple(map(tuple,w['edges']));m=w['orientation_mask'];h=w['host_out_bits']
  require(all(not(h[i]>>i&1) for i in range(9)),'host loop')
  require(all(((h[i]>>j&1)+(h[j]>>i&1))==1 for i in range(9) for j in range(i+1,9)),'host not tournament')
  b=embedding_count(e,9,m,transitive(9));c=embedding_count(e,9,m,h)
  require(b==w['transitive_count'] and c==w['host_count'] and (b-c)%2,'bad-class count')
  _,out=orient_data(e,9,m)
  require(b==extensions(9,out),'linear-extension independence control')
  observed.add(canon(adj_of(e,9)))
 require(expected==observed,'eleven-class witness coverage')
 # D8, its isolated padding, and both concrete leaf extensions.
 D=((0,1),(0,2),(0,3),(1,4),(1,5),(2,6),(4,7))
 require(aut(adj_of(D,8))==2 and not hereditary(adj_of(D,8)),'D8 structure')
 padding_counts=[]
 for mask in range(128):
  _,out=orient_data(D,9,mask)
  c=extensions(9,out);padding_counts.append(c)
  require(c%2==1,'D8 padding not odd')
 result=[]
 for new_edge,pair_edge in [((0,8),2),((2,8),5)]:
  e=D+(new_edge,);a=adj_of(e,9)
  require(aut(a)==2 and afree(a)==64 and not hereditary(a) and criterion(a),'extension structure')
  iso=set();odd_iso=set();odd=zero=0
  for mask in range(256):
   aa,out=orient_data(e,9,mask);code,g=ocanon_aut(aa,out);iso.add(code)
   count=extensions(9,out)
   expected_parity=((mask>>pair_edge)^(mask>>7))&1
   require(count%2==expected_parity,'exact cherry selector')
   require(g==(1 if expected_parity else 2),'oriented group selector')
   if expected_parity:odd+=1;odd_iso.add(code)
   else:zero+=1
  require((odd,zero,len(odd_iso),len(iso))==(128,128,64,192),'orientation class totals')
  result.append({'new_edge':new_edge,'odd_labelled_masks':odd,'zero_labelled_masks':zero,'odd_isomorphism_classes':len(odd_iso),'all_isomorphism_classes':len(iso)})
 print(json.dumps({'status':'PASS','unlabelled_nine_trees':47,'good_nine_trees':36,'nonhereditary_good_nine_trees':12,'new_good_decomposition':new_good_summary,'bad_class_witnesses':len(observed),'eight_tree_edge_cuts':cuts8,'nine_tree_two_edge_cuts':second_cuts,'critical_leaf_selector_masks':selector_masks,'D8_padding_orientation_masks':len(padding_counts),'D8_padding_transitive_count_range':[min(padding_counts),max(padding_counts)],'extensions':result},indent=2))
if __name__=='__main__':main()
