import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from pathlib import Path
from fractions import Fraction
import json,math
from core import states,compose
root=Path(__file__).parent;one=[(1,2,b,b,0,b,b) for b in (0,1)]
def need(x,m):
 if not x:raise RuntimeError(m)
def fold(ps):
 s=ps[0]
 for z in ps[1:]:s=compose(s,z)
 return s
out={}
w=json.loads((root/'LOOP_WITNESS.json').read_text());core=states(w['adj'])[w['mask']];Y=fold([one[int(b)] for b in w['Y']]);s=fold([one[1],one[0],core]);start=s[:5];res=[]
for k in range(4):
 res.append((4*s[2]+s[0]-s[4])%190)
 if k<3:s=tuple(z%190 for z in compose(s,Y))
need(tuple(z%190 for z in start)==s[:5],'not period3');need(res==[149,169,119,149],'residues');need(all(math.gcd(r,190)==1 for r in res),'gcd')
out['H95_period']={'modulus':190,'period':3,'residues':res[:3],'gcd_with_step':1,'coordinates_checked':['h','q','A','G','T'],'Y_start':Y[5]}
w=json.loads((root/'UNIT_LOOP_WITNESS.json').read_text());core=states(w['adj'])[w['mask']];Y=fold([one[int(b)] for b in w['Y']]);B=Y[1];R=Y[2];V=Y[3];TY=Y[4];need(Y[5]==0,'Y starts nonzero');beta=Fraction(37*R*V,B-1)-37*TY
alpha=6996+Fraction(148*R-1749*V,B-1)-Fraction(37*R*V,(B-1)**2);gamma=Fraction(-939)-Fraction(148*R-1749*V,B-1)+Fraction(37*R*V,(B-1)**2)
need(alpha>0,'alpha');checks=[];s=fold([one[1],one[0],core])
for k in range(9):
 D=4*s[2]+37-s[4];need(alpha*B**k+beta*k+gamma==D,'baseline');checks.append(k);s=compose(s,Y)
out['unit_baseline']={'B_exponent':len(w['Y']),'alpha_positive':True,'alpha_floor':alpha.numerator//alpha.denominator,'beta_floor':beta.numerator//beta.denominator,'gamma_floor':gamma.numerator//gamma.denominator,'checked_k':checks,'formula':'alpha*B^k+beta*k+gamma','R':R,'V':V,'T_Y':TY}
print(json.dumps(out,indent=2,sort_keys=True))
