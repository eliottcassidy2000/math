#!/usr/bin/env python3
"""Independent exhaustive symbolic and projective certificates for all 12 AP balances."""
import itertools,json
from pathlib import Path
import sympy as s
u,v=s.symbols('u v');O=Path(__file__).resolve().parent
A=[9*(18*u*u+50*u*v+37*v*v)**2,36*(u*u+13*u*v+16*v*v)**2,36*(9*u*u+53*u*v+64*v*v)**2,4*(3*u*u-25*u*v-48*v*v)**2]
def need(c,msg):
 if not c:raise ValueError(msg)
primes={(0,1,2):[13,13,13],(0,1,3):[7,7,7],(0,2,3):[11,11,11],(1,2,3):[13,17,13]}
cases=[]
for triple in itertools.combinations(range(4),3):
 for mid,p in zip(triple,primes[triple]):
  ends=[i for i in triple if i!=mid];bal=s.expand(2*A[mid]-sum(A[i] for i in ends));constant,factors=s.factor_list(bal)
  need(s.expand(constant*s.prod(f**e for f,e in factors)-bal)==0,'factorization failed')
  lines=[];residual=[]
  for f,e in factors:
   d=s.Poly(f,u,v).total_degree()
   if d==1:
    need(s.expand(f-(2*u+3*v))==0,'unexpected linear factor');lines.append(str(f))
   else:
    rows=[int(f.subs({u:t,v:1})%p) for t in range(p)];inf=int(f.subs({u:1,v:0})%p)
    need(all(rows) and inf!=0,'projective root found')
    need(d in [3,4],'unexpected residual degree')
    residual.append({'form':str(f),'degree':d,'exponent':e,'prime':p,'affine_residues':rows,'infinity_residue':inf})
  need(len(lines)==(1 if triple==(0,1,3) else 0),'unexpected number of linear factors')
  need(len(residual)==1,'unexpected residual count')
  cases.append({'triple_zero_based':triple,'middle_zero_based':mid,'balance':str(bal),'factorization':str(s.factor(bal)),'constant':int(constant),'linear_factors':lines,'residual_certificates':residual})
expected=[900,900,12996,900]
actual=[int(a.subs({u:-3,v:2})) for a in A];need(actual==expected,'boundary mismatch')
# Vanishing on the full rational line follows by homogeneity; check symbolically too.
t=s.symbols('t');need(all(s.expand(a.subs({u:-3*t,v:2*t})-n*t**4)==0 for a,n in zip(A,expected)),'full exceptional line mismatch')
cert={'source':'https://arxiv.org/html/2508.07806v1#S1.E1','indexing':'zero-based','cases':cases,'exceptional_line':'2u+3v=0','A_at_minus3_2':actual,'checks_passed':True}
(O/'ap-certificate.json').write_text(json.dumps(cert,indent=2)+'\n')
print(json.dumps({'checks_passed':True,'cases':len(cases),'A_at_minus3_2':actual}))
for c in cases:print(c['triple_zero_based'],c['middle_zero_based'],c['factorization'])
