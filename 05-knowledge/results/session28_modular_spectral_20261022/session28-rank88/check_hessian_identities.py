"""Exact normalizations only; not a computation of the rank-88 angular spectrum."""
import sympy as s,json
from pathlib import Path
def check(x,label):
 if not x: raise AssertionError(label)
n,r2,T,S2=s.symbols('n r2 T S2',positive=True)
Q2=S2-T*r2/n
laplace=s.factor(8*S2-4/(n+4)*2*(n+4)*Q2-2*T/(n*(n+2))*4*(n+2)*r2)
check(laplace==0,'quartic harmonic decomposition Laplacian')
n=s.Integer(88);N=s.Integer(168498000);minimum=s.Integer(8);dimension=n*(n+1)/2-1
constant=2*minimum**2*N/(n*(n+2))
check(constant==2723200,'angular scalar constant')
check(N*minimum**2*(1-1/n)-dimension*constant==0,'angular operator has zero trace')
check(minimum*N/n==15318000,'inverse-temperature leading correction')
# Valence calculation for every residue class of even unimodular ranks.
rows=[]
for residue in [0,8,16]:
 rows.append({'rank_mod24':residue,'degree2_quotient_weight':residue//2+2-12,'degree4_quotient_weight':residue//2+4-12})
check([r['degree4_quotient_weight'] for r in rows]==[-8,-4,0],'quartic obstruction classification')
out={'status':'EXACT normalizations and identities; no angular eigenvalue computation','harmonic_laplacian':str(laplace),'V_dimension':int(dimension),'angular_constant':int(constant),'asymptotic_inverse_pi_t_coefficient':15318000,'residue_classes':rows}
print(json.dumps(out,indent=2));Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2)+'\n')
