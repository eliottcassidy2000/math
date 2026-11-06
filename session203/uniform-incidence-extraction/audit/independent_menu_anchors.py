#!/usr/bin/env python3
from itertools import product
from math import comb
from random import Random
import json

def need(c,m):
 if not c:raise RuntimeError(m)
graphs=paths=0
for t in range(5):
 for mid in product(range(2),repeat=t):
  bits=[0,0,1,1]+list(mid)+[0,0,1,1];m=len(bits)
  for pa,pc,ac in product(range(2),repeat=3):
   P=range(4);S=range(m-4,m)
   pM,pR=next((i,j) for i in P for j in P if i!=j and bits[i]^pa and bits[j]^pc)
   sL,sM=next((i,j) for i in S for j in S if i!=j and not(bits[i]^pa) and not(bits[j]^pc))
   a=[0]*(m+2)
   for i in range(m):
    for j in range(i+1,m):a[i]|=1<<j
   a[m if ac else m+1]|=1<<(m+1 if ac else m)
   for v,p in ((m,pa),(m+1,pc)):
    for i,b in enumerate(bits):
     if b^p:a[v]|=1<<i
     else:a[i]|=1<<v
   words=set()
   for ass in product(range(3),repeat=t):
    L=[i for i in P if i not in (pM,pR)]+[4+i for i,b in enumerate(ass) if b==0]+[sL]
    M=[pM]+[4+i for i,b in enumerate(ass) if b==1]+[sM]
    R=[pR]+[4+i for i,b in enumerate(ass) if b==2]+[i for i in S if i not in (sL,sM)]
    word=tuple(sorted(L)+[m]+sorted(M)+[m+1]+sorted(R))
    need(sorted(word)==list(range(m+2)),'vertex partition')
    need(all(a[word[i]]>>word[i+1]&1 for i in range(len(word)-1)),'literal injected path')
    words.add(word);paths+=1
   need(len(words)==3**t,'injectivity');graphs+=1
rows=[];histograms=0
for q in range(1,9):
 counts=[]
 def rec(left,k,acc):
  global histograms
  if k==1:
   v=acc+[left];histograms+=1;counts.append(sum(n*(n-1)//2 for n in v));return
  for n in range(left+1):rec(left-n,k-1,acc+[n])
 rec(2*q,q,[])
 need(min(counts)==q,'sharp column pair floor')
 rows.append({'q':q,'histograms':len(counts),'minimum_equal_color_pairs':min(counts),'grid_rows':2*q,'grid_columns':q*(2*q-1)+1,'free_coordinates':(q+1)*(2*q-1)})
need(3**20==3486784401 and 6*2**29+1==3221225473,'m28 constants')
for m in range(28,101):need(3**(m-8)>6*2**(m+1)+1,'all tested growth thresholds')
rng=Random(203);grids=0
for q in range(1,11):
 for trial in range(100):
  labels={};found=False
  for col in range(q*(2*q-1)+1):
   color=[rng.randrange(q) for _ in range(2*q)]
   for i in range(2*q):
    for j in range(i+1,2*q):
     if color[i]==color[j]:
      key=(i,j,color[i])
      if key in labels and labels[key]!=col:found=True
      labels[key]=col
  need(found,'monochromatic rectangle');grids+=1
print(json.dumps({'status':'PASS','literal_graphs':graphs,'explicit_distinct_path_injections':paths,'free_middle_bits_max':4,'color_histograms':histograms,'column_controls':rows,'random_colored_grids':grids,'m28_path_floor':3**20,'m28_plus_target_upper_bound':6*2**29+1,'scope':'independent finite controls; all-order claims rest on the explicit injection and pigeonhole argument'},indent=2,sort_keys=True))
