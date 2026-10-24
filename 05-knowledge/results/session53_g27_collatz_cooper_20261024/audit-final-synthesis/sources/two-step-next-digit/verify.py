import json,math
from pathlib import Path
D=json.loads((Path(__file__).resolve().parent/'design.json').read_text());N=D['coefficient_bound']
def require(ok,msg):
    if not ok:raise RuntimeError(msg)
c=[1,2]
for n in range(1,N):
    q,r=divmod((26*n*n+13*n+2)*c[n]+3*(3*n-1)*(3*n-2)*c[n-1],(n+1)**2)
    require(r==0,'recurrence integrality');c.append(q)
rows=[];total=0
for p in range(3,D['odd_prime_bound']+1,2):
    if any(p%d==0 for d in range(2,math.isqrt(p)+1)):continue
    q,r=divmod(c[p*p]-2,p*p)
    require(r==0,'calibration divisibility')
    a=q*pow(2,-1,p)%p
    counts={'all':0,'n_divisible_by_p':0,'c_n_divisible_by_p':0,'both_divisible_by_p':0}
    first_linear=None;witnesses={}
    for n in range(2,N//(p*p)+1):
        residual=(c[p*p*n]-(1+p*p*a*n*n)*c[n])%(p**3)
        require(residual==0,f'fixed quadratic prediction failed p{p} n{n}')
        total+=1;counts['all']+=1
        for key,condition in [('n_divisible_by_p',n%p==0),('c_n_divisible_by_p',c[n]%p==0),('both_divisible_by_p',n%p==0 and c[n]%p==0)]:
            if condition:
                counts[key]+=1
                if key not in witnesses:witnesses[key]=n
        linear=(c[p*p*n]-(1+p*p*a*n)*c[n])%(p**3)
        if linear and first_linear is None:first_linear={'n':n,'residual_mod_p3':linear}
    require(counts['all']>=1,'no holdout')
    rows.append({'p':p,'a_p_calibrated_only_at_n1':a,'holdout_counts':counts,'first_subgroup_indices':witnesses,'first_linear_model_failure':first_linear})
require(rows[0]['first_linear_model_failure']=={'n':2,'residual_mod_p3':18},'p3 hostile control')
print(json.dumps({'status':'FINITE-EXACT PASS; general formula remains OPEN','coefficient_bound':N,'total_holdouts':total,'rows':rows},indent=2))
