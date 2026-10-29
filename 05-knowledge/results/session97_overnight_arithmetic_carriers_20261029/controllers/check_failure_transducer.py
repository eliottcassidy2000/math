#!/usr/bin/env python3
from fractions import Fraction as F
from check_integer_termination import step,require
from check_exceptional_cylinders import WORDS,chart,fast_run
import json

def main():
 rows=[]
 for i,u in enumerate(WORDS):
  a=chart(u)
  for j,v in enumerate(WORDS):
   b=chart(v);w=u+v;c=chart(w)
   require(c['C']==3**b['J']*a['C']+(1<<a['D'])*b['C'],'composition carry')
   for k in [0,1,2,7,123]:
    n=c['residue']+k*c['modulus'];x,uw,_=fast_run(n);y,vw,_=fast_run(x)
    require(uw==u and vw==v and n<x<y,'guarded port transition')
    require(n%a['modulus']==a['residue'] and x%b['modulus']==b['residue'],'pulled guards')
   rows.append({'from':i,'to':j,**c})
 weights=[F(1,1<<chart(w)['D']) for w in WORDS]
 rho=sum(weights);qj=sum(p*len(w) for p,w in zip(weights,WORDS));qd=sum(p*sum(w) for p,w in zip(weights,WORDS))
 require(rho==F(1059,2048),'rho')
 require(sum(F(1,1<<r['D']) for r in rows)==rho*rho,'joint source mass')
 print(json.dumps({'rows':rows,'failure_mass':str(rho),
 'expected_failed_J_before_exit':str(qj/(1-rho)),
 'expected_failed_D_before_exit':str(qd/(1-rho)),
 'expected_attempts_to_syntactic_exit':str(1/(1-rho))},sort_keys=True))
if __name__=='__main__':main()
