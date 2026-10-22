#!/usr/bin/env python3
import itertools as it,json,math
from pathlib import Path
from check_model import data,tournament,hp,factor

def centralizer(m,q,lab,odd=False,mate_flips=()):
 T,s=tournament(m,q,lab,odd)
 for v in mate_flips:
  a,b=2*v,2*v+1;T[a]^=1<<b;T[b]^=1<<a
 autos=[]
 # Fold-compatible automorphisms commuting with sigma cannot flip pair signs.
 for p in it.permutations(range(m)):
  lift=tuple([2*p[v//2]+v%2 for v in range(2*m)]+([2*m] if odd else []))
  if hp.auto(T,lift):autos.append(p)
 f=hp.dp(T,s)['fixed_hp'];assert f%len(autos)==0
 return {'m':m,'q':q,'lab':lab,'odd':odd,'mate_flips':list(mate_flips),'fixed_hp':f,'centralizer':len(autos),'blue':f//len(autos),'pair_automorphisms':autos,'outmasks':T,'sigma':s}

def binary(h):
 assert h>0 and h%2
 m=(h+1).bit_length()-1
 if 2**m<h+1:m+=1
 es=list(it.combinations(range(m),2));q=(1<<len(es))-1;digits=(h-1)//2
 lab=sum((1<<k) for k,(u,v) in enumerate(es) if digits>>(v-1)&1)
 return m,q,lab

def main():
 records=[]
 for m in range(1,6):
  for h in range(1,2**m,2):
   # use requested m even if the minimal binary construction has fewer pairs
   q=(1<<(m*(m-1)//2))-1;digits=(h-1)//2
   lab=sum(1<<k for k,(u,v) in enumerate(it.combinations(range(m),2)) if digits>>(v-1)&1)
   for odd in (False,True):
    r=centralizer(m,q,lab,odd);assert r['fixed_hp']==h and r['centralizer']==1
 for h in (7,21,31,63,127):
  m,q,lab=binary(h)
  for odd in (False,True):
   r=centralizer(m,q,lab,odd);assert r['fixed_hp']==h and r['centralizer']==1;records.append(r)
 records.extend([centralizer(5,0,150),centralizer(5,2,330),centralizer(5,2,330,True),centralizer(5,2,330,True,(0,))])
 Path(__file__).with_name('blue_checks.json').write_text(json.dumps(records,indent=2)+'\n')
 print([(r['m'],r['fixed_hp'],r['centralizer'],r['blue'],r['odd']) for r in records])
if __name__=='__main__':main()
