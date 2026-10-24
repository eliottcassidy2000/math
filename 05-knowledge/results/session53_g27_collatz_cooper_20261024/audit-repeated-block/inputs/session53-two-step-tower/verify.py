import json,math
N=20000

def require(ok,msg):
    if not ok: raise RuntimeError(msg)
def vp(d,p):
    require(d!=0,'unexpected zero difference')
    v=0
    while d%p==0: v+=1;d//=p
    return v,d%p
c=[1,2]
for n in range(1,N):
    q,r=divmod((26*n*n+13*n+2)*c[n]+3*(3*n-1)*(3*n-2)*c[n-1],(n+1)**2)
    require(r==0,'nonintegral recurrence')
    c.append(q)
primes=[p for p in range(2,math.isqrt(N)+1) if all(p%d for d in range(2,math.isqrt(p)+1))]
rows=[];total=0
for p in primes:
    r=2
    while p**r<=N:
        values=[];divisible=0
        for n in range(1,N//(p**r)+1):
            val,unit=vp(c[p**r*n]-c[p**(r-2)*n],p)
            a=0;m=n
            while m%p==0:a+=1;m//=p
            require(val>=2*r-2,f'failure p{p},r{r},n{n}')
            require(val>=2*(r+a)-2,f'shifted failure p{p},r{r},n{n}')
            divisible+=bool(a);values.append((val,n,unit));total+=1
        v,n,u=min(values)
        v1,n1,u1=values[0]
        rows.append({'p':p,'depth':r,'n_count':len(values),'p_divisible_n_count':divisible,'required_valuation':2*r-2,'minimum_valuation':v,'minimum_witness_n':n,'minimum_witness_unit_mod_p':u,'n1_valuation':v1,'n1_unit_mod_p':u1})
        r+=1
require(all(x['n1_valuation']==2*x['depth']-2 for x in rows if x['p']==3),'p3 sharpness')
require(all(x['minimum_valuation']==2*x['depth']-1 for x in rows if x['p']==2 and x['depth']>=3),'p2 finite extra digit')
print(json.dumps({'status':'FINITE-EXACT PASS; infinite-depth hypothesis remains OPEN','coefficient_bound':N,'all_primes_up_to':math.isqrt(N),'comparisons':total,'rows':rows},indent=2))
