"""Independent exact constants for n>=1500 eventual sharpness.
This checker certifies rational endpoints and monotonic derivative bounds;
it does not enumerate any weighted matrix or infer an all-order fact from samples.
"""
from fractions import Fraction as F
from math import factorial
from hashlib import sha256
import json
N=1500;rows=[]
def need(x,msg):
 if not x:raise RuntimeError(msg)
def positive(name,value,extension):
 need(value>0,name);rows.append({'name':name,'positive_margin':str(value),'uniform_extension':extension})
partial=sum((F(4)**j/F(factorial(j)) for j in range(8)),F(0))
need(partial==F(16319,315),'exp4 partial exact value')
positive('c0>49/200',partial-50,'Positive Taylor remainder gives exp4>partial>50; hence c0>(1-1/50)/4.')
positive('kappa<1/187',F(1,187)-F(8,N-1),'8/(n-1) decreases on n>=1500.')
positive('kappa<c0',F(49,200)-F(1,187),'Use the preceding strict upper/lower bounds.')
positive('heavy_selected_weights_exceed_kappa',F(1,44)-F(1,187),'eta=1/44; independent of n once kappa<1/187.')
positive('V_upper_n/7',F(1,7)-F(1,8)-F(4,N-1),'After dividing by n, 1/8+4/(n-1) decreases.')
# Lower mass after deleting all light selected edges; preserve exact n-1 factors.
slope=F(49,1600)-F(1,88)-F(1,120)
need(slope==F(577,52800),'heavy-mass slope exact value')
endpoint=F(49*(N-1),1600)-F(16*N,N-1)-F(N,88)-F(N,120)
need(endpoint==F(9253539,26382400),'heavy-mass endpoint exact value')
positive('retained_heavy_mass>n/120',endpoint,'f(n)=slope*n-49/1600-16n/(n-1); derivative=slope+16/(n-1)^2>0.')
positive('heavy_mass_monotonic_slope',slope,'Additional derivative term16/(n-1)^2 is nonnegative.')
# This endpoint is integer; extension uses the logarithmic derivative bound below.
positive('epsilon_endpoint',F(4**50-1_500_000*3**50),'n=1500; equivalent to n*(3/4)^(n/30)<1/1000.')
positive('epsilon_decreases_after_multiplying_n',F(1,120)-F(1,N),'log(4/3)>1/4 by integral; derivative log[n*(3/4)^(n/30)]<1/n-1/120<0.')
need(F(4,1000*7)==F(1,1750),'epsilon error bound')
eta=F(1,44);absorption=F(49,200)-4*eta-F(1,1750)
positive('light_mass_absorption>3/20',absorption-F(3,20),'Since L>=1, the fixed error<=1/1750 is absorbed into L/1750.')
need(absorption-F(3,20)==F(271,77000),'light absorption exact excess')
need(4/F(3,20)==F(80,3),'L<80D/3 coefficient')
linear=44*F(1,187);radicand=4*F(80,3)*F(1,187)
need(linear==F(4,17) and radicand==F(320,561),'contraction exact terms')
positive('contraction_remaining_positive',1-linear,'Squaring the remaining positive term is legitimate.')
positive('strict_contraction_square_gap',(1-linear)**2-radicand,'Thus4/17+sqrt(320/561)<1 without floating point.')
need((1-linear)**2-radicand==F(137,9537),'contraction exact squared gap')
# All applied matching coefficients have r_h>=4M_h>n/30>=50, and L>=1 forces
# r_h<=floor(n/2)-1, hence no empty/perfect-matching exceptional coefficient.
need(N>=4 and F(N,30)>=50,'nonempty matching and global theorem domain')
result={'status':'PASS','cutoff':N,'rational_certificates':len(rows),'certificates':rows,'scope':'Exact endpoints plus stated analytic monotonic extensions; no tournament/graph census; active checks under Python -O'}
result['semantic_sha256']=sha256(json.dumps(result,sort_keys=True).encode()).hexdigest()
print(json.dumps(result,indent=2,sort_keys=True))
