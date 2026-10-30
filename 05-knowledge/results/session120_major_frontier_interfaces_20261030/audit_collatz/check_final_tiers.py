"""Final independent exact tier checks; supersedes weaker count bounds."""
from fractions import Fraction as F
import json

def need(p,msg):
    if not p: raise RuntimeError(msg)
gamma=F(5878847,19531250)
kappa=F(17587861,48828125)
D=8*5**10
need(F(10,8)<1+gamma<F(11,8),'first possible core eighth')
need((4*F(11,8)+1+gamma)/5==1+kappa,'quantized core transfer')
T=2*5**9
S=sum(T//(2*5**(e-1)) for e in range(1,11))
need((F(5,8)-F(1,5*D))*T-S==F(6,25),'allones tier increment')
caps=[(j+8-11)//2 for j in range(3,10)]
need(caps==[0,0,1,1,2,2,3],'sharpened per-j caps')
families=[(0,1,2),(1,1,1),(2,1,0),(3,1,0),(0,2,1),(1,2,0)]
rows=[]
for a,b,K in families:
    base=F(a,5)+b*kappa
    need(base+F(6*K,25)<1,'permitted last integer tier')
    need(base+F(6*(K+1),25)>=1,'next tier paid')
    rows.append({'a12':a,'b11':b,'max_internal_K_for_unpaid_candidate':K,
                 'next_tier_excess':str(base+F(6*(K+1),25)-1)})
print(json.dumps({'status':'PASS','superseding_per_j_caps':caps,'core_price_bound':str(1+kappa),
 'tier_increment':'6/25','family_tiers':rows,'additional_families_closed':0},indent=2))
