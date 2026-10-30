from collections import Counter
from itertools import product
import json
from verify_norms import table,rectangle,line_equation,require

def main():
 ls=[(1,2,3),(1,4,5),(1,6,7),(2,4,6),(2,5,7),(3,4,7),(3,5,6)]
 rows=[((2,4,6),(2,5,7),(3,4,7),(3,5,6)),((1,4,5),(1,6,7),(2,4,6),(2,5,7)),((1,2,3),(1,6,7),(2,4,6),(3,4,7))]
 expected=Counter({((),0):70})
 for bits in product(range(2),repeat=3):
  if not any(bits):continue
  supp=set();rhs=0
  for bit,row in zip(bits,rows):
   if bit:supp.symmetric_difference_update(row);rhs^=1
  expected[(tuple(sorted(supp)),rhs)]=6
 actual=Counter();rects=[]
 for t in range(1,8):
  for a in range(8):
   if a>(a^t):continue
   for b in range(8):
    if b>(b^t):continue
    rects.append((t,a,b));s,r=line_equation(8,t,a,b);actual[(tuple(sorted(s)),r)]+=1
 require(actual==expected,'symbolic rectangle multiplicity')
 hist=Counter()
 for z in product(range(2),repeat=7):
  f=table(8,z);bad=sum(not rectangle(f,*r) for r in rects);hist[bad]+=1
  syndrome=tuple(sum(z[ls.index(l)] for l in row)%2^1 for row in rows)
  require(bad==(24 if any(syndrome) else 0),'defect by syndrome')
  # Independently expand the full polynomial and sum squared coefficients.
  coeff=Counter()
  for a,b,c,d in product(range(8),repeat=4):
   if a^b==c^d:coeff[tuple(sorted((a,c)))+tuple(sorted((b,d)))]+=(-1)**(f[a][b]+f[c][d])
  for a,b in product(range(8),repeat=2):coeff[(a,a,b,b)]-=1
  require(sum(v*v for v in coeff.values())==16*bad,'mean-square coefficient sum')
 require(hist=={0:16,24:112},'unexpected histogram')
 print(json.dumps({'status':'PASS','symbolic_rectangles':112,'tautologies':70,'nonzero_syndrome_functionals':7,'multiplicity_each':6,'tables':128,'defect_histogram':dict(hist),'invalid_mean_square':384,'normalized_rademacher_mean_square':'3/32'},sort_keys=True,indent=2))
if __name__=='__main__':main()
