#!/usr/bin/env python3
import itertools,math,json
from pathlib import Path

def necklace(d):return min(d[k:]+d[:k] for k in range(len(d)))
def displacement(q):return tuple((q[(i+1)%len(q)]-q[i])%len(q) for i in range(len(q)))
def converse(q):n=len(q);return tuple((-q.index((-i)%n))%n for i in range(n))
def gauge(q):
 n=len(q);return min(tuple((q[(i+a)%n]+b)%n for i in range(n)) for a in range(n) for b in range(n))
rows=[]
for p in(3,5,7):
 classes={};valid=0;wordclasses=set();wordchecks=0
 for q in itertools.permutations(range(p)):
  key=necklace(displacement(q));g=gauge(q)
  if key in classes and classes[key]!=g:raise RuntimeError('necklace class collision')
  classes[key]=g
 for d in itertools.product(range(1,p),repeat=p):
  wordchecks+=1;partials=[0]
  for x in d:partials.append((partials[-1]+x)%p)
  if partials[-1]!=0 or len(set(partials[:-1]))!=p:continue
  q=tuple(partials[:-1]);valid+=1;wordclasses.add(necklace(d))
  if displacement(q)!=d:raise RuntimeError('reconstruction')
 if valid!=math.factorial(p-1) or wordclasses!=set(classes):raise RuntimeError('realizable necklace census')
 fixed=sum(len(set(d))==1 for d in wordclasses)
 if fixed!=p-1:raise RuntimeError('fixed words')
 rows.append({'p':p,'all_nonzero_words_tested':wordchecks,'realizable_words_q0_zero':valid,'necklace_classes':len(classes),'constant_necklaces':fixed})
q=(0,2,4,1,3)
if necklace(displacement(q))==necklace(displacement(converse(q))):raise RuntimeError('converse pair merged ordinarily')
if necklace(displacement(q)[::-1])!=necklace(displacement(q)):raise RuntimeError('constant word reversal')
out={'status':'PASS','rows':rows,'converse_not_word_reversal_witness':{'q':q,'word':displacement(q),'converse_word':displacement(converse(q))}}
Path(__file__).with_name('displacement-necklaces-audit.json').write_text(json.dumps(out,indent=2)+'\n');print(out)
