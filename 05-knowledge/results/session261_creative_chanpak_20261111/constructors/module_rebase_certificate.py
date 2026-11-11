#!/usr/bin/env python3
"""Exact explicit profile trade and finite all-length rooted certificates."""
from module_tools import *
from itertools import permutations
from math import factorial,gcd
import json
HERE=Path(__file__).resolve().parent

def literal_covers(T):
 n=len(T);c=[0]*n
 for p in permutations(range(n)):
  bad=sum(1<<i for i in range(n-1) if not T[p[i]]>>p[i+1]&1)
  for cuts in range(1<<(n-1)):
   if cuts&bad==bad:c[cuts.bit_count()]+=1
 return tuple(c)

def root_certificate(T,v,c1,c2):
 # n-1 exterior vertices => at most n module runs in a Hamilton path.
 cs=[]
 for r in range(len(T)+1):
  Z,_=m.inflate(T,0,v,r);rows=list(Z);start=len(T)-1;mask=((1<<r)-1)<<start
  for a in range(start,len(rows)):rows[a]&=~mask
  h=sum(m.endpoints(tuple(rows))[-1])
  if h%factorial(r):raise RuntimeError('unlabelled run divisibility')
  cs.append(h//factorial(r))
 expected=[cs[0],c1,c2]+[0]*(len(T)-2)
 if cs!=expected:raise RuntimeError(('root coefficient mismatch',cs,expected))
 return {'rows':T,'port':v,'run_coefficients':cs}

def main():
 A=(2,4,9,19,39,15);B=(2,4,57,3,11,27);expected=(17,190,828,1752,1800,720)
 for T in (A,B):
  if covers(T)!=expected or literal_covers(T)!=expected:raise RuntimeError('profile mismatch')
 F=(2,4,25,3,11);G=(2,4,9,19,7)
 if covers(F)!=covers(G) or gcd(*covers(F))!=3:raise RuntimeError('hostile mismatch')
 certs=[root_certificate(A,1,17,9),root_certificate(B,0,17,8)]
 # Concrete distinction when marks lie inside rather than outside the module.
 if m.data(B,6)!=(17,16,48):raise RuntimeError('selected interior control')
 if any(m.data(A,O)[0]-m.data(A,O)[1]==1 for O in range(64)):raise RuntimeError('interior hostile')
 out={'status':'PASS','six_profile':expected,'six_content':gcd(*expected),'root_certificates':certs,'five_hostile_profile':covers(F),'five_hostile_content':3,'scope':'All-host theorem is proved in MODULE_REBASE.md; finite skeleton coefficients certify the two all-length clone identities.'}
 HERE.joinpath('module_rebase_receipt.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
if __name__=='__main__':main()
