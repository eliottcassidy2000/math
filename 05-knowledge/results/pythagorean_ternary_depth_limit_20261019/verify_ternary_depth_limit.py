#!/usr/bin/env python3
"""Exact cylinder, suffix-classifier, and ternary-limit controls."""
from fractions import Fraction as Q
import json
CHECKS=0
def check(condition,message):
 global CHECKS
 CHECKS+=1
 if not condition:raise RuntimeError(message)
MAPS={'L':(2,-1,1,0),'M':(2,1,1,0),'R':(1,2,0,1)}
def mul(x,y):
 a,b,c,d=x;e,f,g,h=y
 return a*e+b*g,a*f+b*h,c*e+d*g,c*f+d*h
def det(x):a,b,c,d=x;return a*d-b*c
def val(x,q):a,b,c,d=x;return (a*q+b)/(c*q+d)
def polynomial(q):
 n,m=q.numerator,q.denominator
 return n**4-8*n**3*m+2*n*n*m*m+8*n*m**3+m**4
def compare_mu(q):
 if q is None:return 1
 if q<=1:return -1
 if q>=2:return 1
 z=polynomial(q);check(z!=0,'mu is not a rational endpoint')
 return -1 if z>0 else 1
def compare_lambda(q):
 if q is None:return 1
 if q<=3:return -1
 if q>=8:return 1
 z=polynomial(q);check(z!=0,'lambda is not a rational endpoint')
 return 1 if z>0 else -1

def boundary_prefixes(compare,depth):
 matrix=(1,0,0,1);orientation=1;index=0;word='';digits='';data=[]
 for d in range(1,depth+1):
  if orientation*compare(val(matrix,Q(2)))>0:s='L'
  elif orientation*compare(val(matrix,Q(3)))>0:s='M'
  else:s='R'
  a='LMR'.index(s)
  if orientation<0:a=2-a
  index=3*index+a;word+=s;digits+=str(a)
  matrix=mul(matrix,MAPS[s]);orientation*=det(MAPS[s])
  check(det(matrix)==orientation,'orientation from matrix')
  data.append({'depth':d,'index':index,'word':word,'digits':digits,'matrix':matrix})
 return data

def cylinder_type(matrix):
 a,b,c,d=matrix
 x=Q(a+b,c+d);y=Q(a,c) if c else None
 if y is None:lo,hi=x,None
 else:lo,hi=min(x,y),max(x,y)
 if compare_mu(hi)<0 or compare_lambda(lo)>0:return 'positive'
 if compare_mu(lo)>0 and compare_lambda(hi)<0:return 'negative'
 return 'unknown'

def main():
 mu=boundary_prefixes(compare_mu,128);lam=boundary_prefixes(compare_lambda,128)
 for x,y in zip(mu,lam):
  d=x['depth']
  check(y['index']==3**d-1-x['index'],'reversed cylinder index symmetry')
  check(y['word']==x['word'].translate(str.maketrans('LMR','RML')),'L-R conjugacy word')
  check(all(int(a)+int(b)==2 for a,b in zip(x['digits'],y['digits'])),'complementary ternary digits')
 # Check the conjugacy itself as rational maps at many test inputs.
 S=(1,1,1,-1)
 for n in range(2,101):
  t=Q(n,1)
  for a,b in [('L','R'),('R','L'),('M','M')]:
   check(val(S,val(MAPS[a],t))==val(MAPS[b],val(S,t)),'exact involution conjugacy')
 levels=[];nodes=[(1,0,0,1)]
 expected={3:2,4:5,5:15,6:44,7:132,8:396,9:1188,10:3566}
 for depth in range(1,11):
  nodes=[mul(mat,MAPS[s]) for mat in nodes for s in 'LMR']
  counts={'positive':0,'negative':0,'unknown':0};actual=0
  for matrix in nodes:
   kind=cylinder_type(matrix);counts[kind]+=1
   node_positive=polynomial(val(matrix,Q(2)))>0
   if kind!='unknown':check(node_positive==(kind=='positive'),'sound cylinder classification')
   actual+=node_positive
  check(counts['unknown']==2,'exactly two ambiguous suffix cylinders')
  check(counts['positive']==2*mu[depth-1]['index'],'exact wholly-positive cylinder count')
  check(counts['positive']<=actual<=counts['positive']+2,'depth-density bracket')
  if depth in expected:check(actual==expected[depth],'earlier exact tree counts')
  levels.append({'depth':depth,'total':3**depth,'actual_positive_nodes':actual,
                 'wholly_positive_cylinders':counts['positive'],'negative_cylinders':counts['negative'],
                 'unknown_cylinders':counts['unknown']})
 D=128;index=mu[D-1]['index'];denominator=3**D;lower=2*index;upper=2*index+2
 scale=10**50;lo=lower*scale//denominator;hi=upper*scale//denominator
 check(lo==hi,'certified 50-digit decimal truncation')
 decimal=f'{lo//scale}.{lo%scale:050d}'
 receipt={'status':'PASS','checks':CHECKS,'levels':levels,
  'limit_certificate':{'depth':D,'lower_numerator':lower,'upper_numerator':upper,
       'common_denominator':denominator,'interval_width':f'2/{denominator}',
       'decimal_truncation_50_places':decimal,'mu_boundary_word':mu[-1]['word'],
       'mu_orientation_ternary_digits':mu[-1]['digits']},
  'scope':'Exact interval and finite controls; the proof of irrationality is in the text.'}
 print(json.dumps(receipt,indent=2,sort_keys=True))
if __name__=='__main__':main()
