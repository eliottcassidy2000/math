"""Exact symbolic trivalent recurrence audit, conditional on deg R6<=47."""
import sympy as s
import json
from pathlib import Path
z=s.symbols('z');p=1-z-z**3;q=1-z**2-z**6;D=p**2*q**2
# A generic center polynomial retains all coefficients through degree 47.
c=s.symbols('c0:48');R=sum(c[i]*z**i for i in range(48))
N=s.expand(z**5*q**2/2-z**5*p**2*q/2+R*p**2*q+z**10*p**2*(1+z+z**3))
checks=[]
def need(v,label):
 if not v:raise RuntimeError(label)
 checks.append(label)
need(s.gcd(p,s.diff(p,z))==1,'p squarefree')
need(s.gcd(q,s.diff(q,z))==1,'p2 squarefree')
need(s.gcd(p,q)==1,'coprime factors')
need(s.degree(D,z)==18,'denominator degree 18')
need(s.rem(q-2*z**4,p,z)==0,'p2=2z4 modulo p')
need(s.rem(p*(1+z+z**3)+2*z**4,q,z)==0,'p(1+z+z3)=-2z4 modulo p2')
need(s.rem(N-z**5*q**2/2,p,z)==0,'generic numerator residue modulo p')
need(s.rem(N-z**10*p**2*(1+z+z**3),q,z)==0,'generic numerator residue modulo p2')
need(s.gcd(z**5*q**2,p)==1,'numerator residue coprime to p')
need(s.gcd(z**10*p**2*(1+z+z**3),q)==1,'numerator residue coprime to p2')
need(s.degree(N,z)==59,'generic numerator degree 59')
need(s.expand(N).coeff(z,59)==-c[47],'top numerator coefficient -c47')
need(s.expand(D).coeff(z,0)==1,'normalized recurrence')
receipt={'status':'PASS','checks':len(checks),'checked':checks,'assumption':'R6 is a polynomial of degree at most 47, supplied by separately reviewed center bound','minimal_eventual_recurrence_order':18,'safe_recurrence_start_n':60,'denominator_coefficients_low_to_high':[int(s.expand(D).coeff(z,i)) for i in range(19)],'no_earliest_onset_claim':True}
Path('session22-fixed-degree-review/trivalent_recurrence_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps(receipt,indent=2))
