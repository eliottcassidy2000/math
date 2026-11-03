from itertools import product
from fractions import Fraction as F
import json
from test import cert, require

def cycle(B,d):
 e=len(B)
 total=sum(map(len,B)); overlaps=[len(B[i]&B[(i+1)%e]) for i in range(e)]
 full=len(set.intersection(*B)); union=len(set.union(*B))
 bound=total-sum(overlaps)+full
 defect=0
 for q in range(d):
  s=[q in b for b in B]
  if any(s) and not all(s):
   runs=sum(s[i] and not s[(i-1)%e] for i in range(e))
   defect+=runs-1
 require(bound-union==defect>=0, "cycle defect")
 return dict(sum_sizes=total,overlaps=overlaps,full_intersection_size=full,cycle_bound=bound,union_size=union,defect=defect)
sets=[{q for q in range(4) if mask>>q&1} for mask in range(16)]
checks=0
for bs in product(sets,repeat=4):cycle(bs,4);checks+=1
examples=[]
for d in [8,9,44,527]:
 e=[4,2,5,3,7,13];c=cert(d,[d*i for i in range(1,9)],e,F(1,15),F(1,d*9))
 c['cycle']=cycle(list(map(set,c['bad_sets'])),d)
 examples.append(c)
c=cert(9,[9*i for i in range(1,7)],[1,3,5,7,13,44,527,8],F(1,15),F(1,63))
# Exact conditional translation 5+8=13 on the eight-sheet clock.
c8=cert(8,[8*i for i in range(1,8)],[1,3,5,7,13,44,527],F(1,15),F(1,64))
b5=set(c8['bad_sets'][2]);b13=set(c8['bad_sets'][4]);require(b13=={(q-5)%8 for q in b5}, "phase translation")
print(json.dumps(dict(exhaustive_four_mask_families=checks,cycle_examples=examples,nine_sheet_eight_exception_example=c,translation_check=True),indent=2))
