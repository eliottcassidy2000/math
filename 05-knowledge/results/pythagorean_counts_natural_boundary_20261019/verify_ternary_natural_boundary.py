#!/usr/bin/env python3
"""Exact controls for the two-wall ternary count and bounded defect theorem.

No floating arithmetic. The natural-boundary assertion is proved in the
companion note using the cited Pólya--Carlson theorem, not inferred by tests.
"""
from fractions import Fraction as Q
from collections import Counter
import argparse, hashlib, json

CHECKS=0
def check(ok,msg):
 global CHECKS
 CHECKS+=1
 if not ok: raise RuntimeError(msg)
MAPS={'L':(2,-1,1,0),'M':(2,1,1,0),'R':(1,2,0,1)}
IDENTITY=(1,0,0,1)
def mul(x,y):
 a,b,c,d=x;e,f,g,h=y
 return (a*e+b*g,a*f+b*h,c*e+d*g,c*f+d*h)
def val(mat,t):
 a,b,c,d=mat
 return (a*t+b)/(c*t+d)
def det(mat):
 a,b,c,d=mat
 return a*d-b*c
def quartic(t):
 n,d=t.numerator,t.denominator
 return n**4-8*n**3*d+2*n*n*d*d+8*n*d**3+d**4
def cmp_mu(t):
 # -1 for t<mu, +1 for t>mu. mu lies in (1,2), and is
 # the unique quartic root there; rational equality is impossible.
 if t<=1:return -1
 if t>=2:return 1
 q=quartic(t)
 check(q!=0,'a rational sample cannot equal the quartic wall')
 return -1 if q>0 else 1
def below_mu(t):return cmp_mu(t)<0

def compile_counts(depth):
 mat=IDENTITY;orientation=1;index=0;counts=[0];digits=[];corrections=[];letters=[]
 for n in range(1,depth+1):
  old=mat
  if orientation*cmp_mu(val(mat,Q(2)))>0:letter='L'
  elif orientation*cmp_mu(val(mat,Q(3)))>0:letter='M'
  else:letter='R'
  digit='LMR'.index(letter)
  if orientation<0:digit=2-digit
  index=3*index+digit
  mat=mul(mat,MAPS[letter]);orientation*=det(MAPS[letter])
  check(det(mat)==orientation,'unimodular branch orientation')
  e=int(below_mu(val(mat,Q(2))))+int(below_mu(val(mat,Q(3))))
  count=2*index+e
  # The only parent cylinder that changes the two-root lower-cut count.
  local_delta=sum(below_mu(val(mul(old,s),t)) for s in MAPS.values() for t in (Q(2),Q(3)))-3*sum(below_mu(val(old,t)) for t in (Q(2),Q(3)))
  check(count-3*counts[-1]==local_delta,'exact one-cylinder defect recurrence')
  check(abs(local_delta)<=2,'sharp signed defect bound')
  check(0<=count<=3**n,'count within level size')
  counts.append(count);digits.append(digit);corrections.append(e);letters.append(letter)
 return counts,digits,corrections,letters,index

def direct_tree(depth):
 # Independent chronological Euclid-pair recursion, not cylinder indices.
 nodes=[(2,1)];counts=[0];level_sizes=[1]
 for n in range(1,depth+1):
  nodes=[child for a,b in nodes for child in ((2*a-b,a),(2*a+b,a),(a+2*b,b))]
  c=0
  for a,b in nodes:
   # Direct Pythagorean a_leg=a^2-b^2, b_leg=2ab, hyp=a^2+b^2.
   x,y,z=a*a-b*b,2*a*b,a*a+b*b
   check(x*x+y*y==z*z,'Pythagorean triple')
   check(z*z!=4*x*y,'strict positive-label boundary')
   c+=z*z>4*x*y
  counts.append(c);level_sizes.append(len(nodes))
 return counts,level_sizes

def run(depth=2048,tree_depth=10):
 allc,alld,alle,allletters,_=compile_counts(depth+1)
 c,d,e,letters=allc[:-1],alld[:-1],alle[:-1],allletters[:-1]
 index=0
 for x in d:index=3*index+x
 check(all(e[n-1]==alld[n] for n in range(1,depth+1)),'boundary correction equals next spatial ternary digit')
 check(all(c[n]-3*c[n-1]==alld[n]-alld[n-1] for n in range(1,depth+1)),'exact digit-difference defect formula')
 prefix=0
 for n in range(depth+1):
  check(c[n]-2*prefix==alld[n],'count tower reconstructs next digit')
  prefix+=c[n]
 direct,sizes=direct_tree(tree_depth)
 check(c[:tree_depth+1]==direct,'boundary compiler equals independent complete tree')
 check(c[:11]==[0,0,0,2,5,15,44,132,396,1188,3566],'inherited first ten counts')
 defects=[c[0]]+[c[n]-3*c[n-1] for n in range(1,len(c))]
 check(set(defects)=={-2,-1,0,1,2},'all five allowed defect values occur in declared prefix')
 # Exact telescoping: sum_{j<=n} defect_j/3^j = c_n/3^n.
 accum=Q(0)
 for n in range(min(depth,128)+1):
  accum+=Q(defects[n],3**n)
  check(accum==Q(c[n],3**n),'regularized series telescoping')
 # Every nonzero arithmetic progression has bounded recurrence defects.
 subsequences=[]
 for stride in range(1,9):
  for offset in range(stride):
   vals=c[offset::stride];base=3**stride
   dd=[vals[i]-base*vals[i-1] for i in range(1,len(vals))]
   check(all(abs(z)<=base-1 for z in dd),'decimated defect bound from geometric sum')
   subsequences.append({'stride':stride,'offset':offset,'checked_terms':len(vals),'bound':base-1})
 # Hostile: rational density 1/2 has bounded integer defects but rational GF.
 control=[(3**n-1)//2 for n in range(40)]
 check(all(control[n]-3*control[n-1]==1 for n in range(1,40)),'rational-density hostile')
 # A long zero-defect stretch does not imply eventual zero or rationality.
 check(all(defects[n]==0 for n in range(11,57)) and defects[57]!=0,'finite plateau hostile')
 raw=json.dumps(c,separators=(',',':')).encode()
 receipt={'status':'PASS','checks':CHECKS,'depth':depth,'complete_tree_depth':tree_depth,
  'complete_tree_total_nodes':sum(sizes),'first_counts':c[:21],
  'first_100_defects':defects[:100],
  'defect_histogram':{str(k):v for k,v in sorted(Counter(defects).items())},
  'counts_json_sha256':hashlib.sha256(raw).hexdigest(),
  'boundary_itinerary':''.join(letters),'ternary_digits':''.join(map(str,d)),
  'last_count':str(c[-1]),
  'limit_interval':{'lower_numerator':str(2*index),'upper_numerator':str(2*index+2),'denominator':str(3**depth)},
  'decimations':subsequences,
  'scope':'Finite exact controls only. Natural boundary, non-D-finiteness, and irrational density are proof statements in the companion note.'}
 return receipt

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--depth',type=int,default=2048);p.add_argument('--tree-depth',type=int,default=10)
 args=p.parse_args();print(json.dumps(run(args.depth,args.tree_depth),sort_keys=True,indent=2))
