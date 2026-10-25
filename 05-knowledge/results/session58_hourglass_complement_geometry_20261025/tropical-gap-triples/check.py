#!/usr/bin/env python3
"""Independent exact PG enumeration and matrix-word checks; no assert guards."""
from itertools import product, combinations
from fractions import Fraction
from math import comb, prod
from pathlib import Path
import json

def require(ok,msg):
 if not ok:raise RuntimeError(msg)
def formula(k):
 n=21**k-6*7**k-3*5**k+14*3**k-6
 require(n%48==0,'integrality');return n//48

def primitive(av):
 if not av:return 0
 k=len(av);num=prod(24*a-4 for a in av)-3*4**k-6*prod(4*a+2 for a in av)+14*2**k
 require(num%48==0,'primitive integrality');return num//48

def general(av):
 num=prod(12*a*a+8*a+1 for a in av)-3*prod(4*a+1 for a in av)-6*prod(2*a*a+4*a+1 for a in av)+14*prod(2*a+1 for a in av)-6
 require(num%48==0,'general integrality');return num//48

def points(k):
 return [v for v in product((0,1,2),repeat=k) if any(v) and next(x for x in v if x)==1]
def support(v):return sum((1<<i) for i,x in enumerate(v) if x)
def allowed_masks(a,b,c):return ((a^b^c)&~(a&b&c))==0

def main():
 direct=[]
 for k in range(6):
  ps=points(k);masks=list(map(support,ps));n=0;pn=0
  for a,b,c in combinations(masks,3):
   good=allowed_masks(a,b,c);n+=good;pn+=good and ((a|b|c)==(1<<k)-1)
  require(len(ps)==(3**k-1)//2,'PG point count')
  require(n==formula(k),'direct PG count k='+str(k));require(pn==primitive((1,)*k),'direct primitive PG count')
  direct.append({'k':k,'points':len(ps),'all_distinct_triples':comb(len(ps),3),'tropical_triples':n,'primitive_tropical_triples':pn})
 alphabet=[r for r in product(range(3),repeat=3) if sum(x!=0 for x in r)!=1]
 require(len(alphabet)==21,'row alphabet')
 matrices=[]
 for k in range(4):
  nz=pair01=rankone=distinct=0
  for rows in product(alphabet,repeat=k):
   cols=[tuple(r[j] for r in rows) for j in range(3)]
   if not all(any(v) for v in cols):continue
   nz+=1
   def eq(u,v):return u==v or u==tuple((-x)%3 for x in v)
   ab,ac,bc=eq(cols[0],cols[1]),eq(cols[0],cols[2]),eq(cols[1],cols[2])
   pair01+=ab;rankone+=ab and ac;distinct+=not(ab or ac or bc)
  require(nz==21**k-3*5**k+2,'nonzero column IE')
  require(pair01==2*(7**k-3**k),'pair event')
  require(rankone==4*(3**k-1),'rank-one intersection')
  require(distinct==48*formula(k),'free quotient')
  matrices.append({'k':k,'nonzero_column_matrices':nz,'first_pair_equal':pair01,'all_pairs_equal':rankone,'distinct_projective_column_matrices':distinct})
 bases=[21,7,5,3,1];coeff=[1,-6,-3,14,-6]
 def pmul(a,b):
  out=[0]*(len(a)+len(b)-1)
  for i,x in enumerate(a):
   for j,y in enumerate(b):out[i+j]+=x*y
  return out
 denominator=[1]
 for a in bases:denominator=pmul(denominator,[1,-a])
 numerator=[Fraction(0)]*5
 for a,c in zip(bases,coeff):
  term=[1]
  for b in bases:
   if a!=b:term=pmul(term,[1,-b])
  numerator=[x+Fraction(c*y,48) for x,y in zip(numerator,term)]
 require(denominator==[1,-37,422,-1982,3801,-2205],'denominator')
 require(numerator==[0,0,4,2,-126],'numerator')
 values=[formula(k) for k in range(41)]
 for k in range(5,41):require(sum(denominator[j]*values[k-j] for j in range(6))==0,'recurrence')
 require(formula(2)==4 and formula(3)==150,'known controls')
 # A forbidden singleton coordinate; a passing triple is not a numerical Schur certificate.
 require(not allowed_masks(1,2,2),'singleton hostile')
 require(allowed_masks(1,2,3),'two-prime positive support control')
 gaps=(4056,3000,3696)
 require(not any(gaps[i]==sum(gaps[j] for j in range(3) if j!=i) for i in range(3)),'projective/additive hostile')
 exponent_checks=[]
 vectors=[()]+[(a,) for a in range(1,8)]+[(1,2),(2,2),(3,3),(1,1,2),(1,2,3),(2,2,2)]
 for av in vectors:
  ls=[12*a*a+8*a+1 for a in av];us=[4*a+1 for a in av];vs=[2*a*a+4*a+1 for a in av];ws=[2*a+1 for a in av]
  predicted=(prod(ls)-3*prod(us)-6*prod(vs)+14*prod(ws)-6)//48
  box=[v for v in product(*(range(-a,a+1) for a in av)) if any(v) and next(x for x in v if x)>0]
  actual=0;primactual=0
  for t in combinations(box,3):
   good=all(sum(abs(v[i])==max(abs(w[i]) for w in t) for v in t)>=2 for i in range(len(av)))
   actual+=good;primactual+=good and all(max(abs(w[i]) for w in t)==av[i] for i in range(len(av)))
  require(actual==predicted,'general exponent box '+str(av));require(primactual==primitive(av),'primitive exponent box')
  require(sum(primitive(tuple(b for b in bv if b)) for bv in product(*(range(a+1) for a in av)))==actual,'hereditary maxima partition')
  require((prod(ls)-3*prod(us)-6*prod(vs)+14*prod(ws)-6)%48==0,'general integrality')
  exponent_checks.append({'exponents':av,'points':len(box),'candidate_triples':actual,'primitive_candidate_triples':primactual})
 for a in range(1,8):
  rows=list(product(range(-a,a+1),repeat=3))
  good=[r for r in rows if sum(abs(x)==max(map(abs,r)) for x in r)>=2]
  require(len(good)==12*a*a+8*a+1,'general local L')
  require(sum(r[0]==0 for r in good)==4*a+1,'general local U')
  require(sum(r[0]==r[1] for r in good)==2*a*a+4*a+1,'general local V')
  require(sum(r[0]==r[1]==r[2] for r in good)==2*a+1,'general local W')
 pbases=[20,6,4,2];pcoeff=[1,-6,-3,14];pden=[1]
 for a in pbases:pden=pmul(pden,[1,-a])
 pnum=[Fraction(-6*x,48) for x in pden]
 for a,c in zip(pbases,pcoeff):
  term=[1]
  for b in pbases:
   if a!=b:term=pmul(term,[1,-b])
  for i,x in enumerate(term):pnum[i]+=Fraction(c*x,48)
 require(pnum==[0,0,4,10,-120] and pden==[1,-32,284,-928,960],'primitive OGF')
 pv=[primitive((1,)*k) for k in range(41)]
 for k in range(5,41):require(sum(pden[j]*pv[k-j] for j in range(5))==0,'primitive recurrence')
 require(sum(pden[j]*pv[4-j] for j in range(5))==-120,'primitive recurrence starting-index hostile')
 for k in range(41):require(sum(comb(k,j)*pv[j] for j in range(k+1))==formula(k),'binomial transform')
 out={'status':'PASS','direct_projective_checks':direct,'independent_matrix_word_checks':matrices,'ogf_numerator':[int(x) for x in numerator],'ogf_denominator':denominator,'initial_values':values[:11],'recurrence_checked_through':40,'density_asymptotic':'(7/9)^k','general_exponent_checks':exponent_checks,'primitive_ogf_numerator':[int(x) for x in pnum],'primitive_ogf_denominator':pden,'primitive_initial_values':pv[:11],'primitive_recurrence_starts_at':5,'primitive_k4_residual':-120}
 Path(__file__).with_name('certificate.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
if __name__=='__main__':main()
