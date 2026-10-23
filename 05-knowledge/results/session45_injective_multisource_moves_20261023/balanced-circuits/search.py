from itertools import combinations
from random import Random
from functools import reduce
from operator import xor
import json
from pathlib import Path
HERE=Path(__file__).resolve().parent
D=json.load((HERE/'pool.json').open())
labels={int(s):tuple(v) for s,v in D['labels'].items()};ids=sorted(labels)
def total(seq):return reduce(xor,(labels[s][0] for s in seq),0),sum(labels[s][1] for s in seq)%128
def primitive(A,B):
 for j in range(1,len(A)//2+1):
  if {total(c) for c in combinations(A,j)} & {total(c) for c in combinations(B,j)}:return False
 return True
rng=Random(4501);found={};attempts={}
for k in range(3,9):
 count=0
 for trial in range(300):
  z=ids[:];rng.shuffle(z);left=z[:10];right=z[10:];bank={}
  for a in combinations(left,k):bank.setdefault(total(a),[]).append(a)
  for b in combinations(right,k):
   for a in bank.get(total(b),[]):
    count+=1
    if primitive(a,b):found[k]={'plus':sorted(a),'minus':sorted(b),'sum':total(a),'slots':2*k,'owner_clock_cost':4*k,'partition_trial':trial,'matching_candidates_tested':count};break
   if k in found:break
  if k in found:break
 attempts[k]=count
print(json.dumps({'found':found,'matching_candidates_tested':attempts},indent=2))
json.dump({'found':found,'matching_candidates_tested':attempts},(HERE/'candidates.json').open('w'),indent=2)
