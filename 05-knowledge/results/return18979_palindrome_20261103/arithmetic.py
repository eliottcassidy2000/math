"""Exact finite controls for the all-modulus palindrome clock and norm charts."""
from pathlib import Path
from math import gcd,isqrt
import json
from sympy import Matrix,ZZ
from sympy.matrices.normalforms import smith_normal_form
HERE=Path(__file__).resolve().parent
def need(b,s):
    if not b:raise ValueError(s)
def factors(n):
    out={};p=2
    while p*p<=n:
        while n%p==0:out[p]=out.get(p,0)+1;n//=p
        p+=1
    if n>1:out[n]=1
    return out
def order(b,m):
    if m==1:return 1
    x=b%m;t=1
    while x!=1:x=x*b%m;t+=1
    return t
def wing(b,x,m):return (b*x*x-x-1)%m
def snf(a):
    D=smith_normal_form(Matrix(a),domain=ZZ);return [abs(int(D[i,i])) for i in range(3)]
def main():
    clocks=0
    for b in range(2,13):
        need(snf([[-1,-b**j,b**(2*j+1)] for j in range(3)])==[1,b-1,b*b*(b-1)**2*(b+1)],'wing Smith formula')
        for M in range(2,301):
            h=0;unit=M
            for p,e in factors(M).items():
                if b%p==0:
                    f=0;q=b
                    while q%p==0:q//=p;f+=1
                    h=max(h,(e+f-1)//f);unit//=p**e
            period=order(b,unit);xs=[];seen={};x=1%M
            while x not in seen:seen[x]=len(xs);xs.append(x);x=x*b%M
            need((seen[x],len(xs)-seen[x])==(h,period),'direct eventual state orbit')
            pairs=[(wing(b,x,M),wing(b,b*x%M,M)) for x in xs]
            need(len(set(pairs))==h+period,'two reads separate every actual clock state');clocks+=1
    need(snf([[1,666*10**j,100*100**j] for j in range(3)])==[1,18,296703000],'Belphegor Smith')
    p=7417;need(all(p%d for d in range(2,isqrt(p)+1)),'trial-division primality7417')
    T=order(10,p);need(T==2472 and all(pow(10,T//q,p)!=1 for q in factors(T)),'full order certificate7417')
    hits=[];bhits=[];x=1
    for m in range(T):
        if wing(10,x,p)==0:hits.append(m)
        y=10*x%p
        if (100*y*y+666*y+1)%p==0:bhits.append(m)
        x=x*10%p
    need(hits==[1312,2395] and not bhits,'complete7417 phase chart')
    nonres=[]
    for q in (7,11,13,101):
        need(41%q not in {x*x%q for x in range(q)},'wing discriminant nonresidue');nonres.append(q)
    valuation_counts={0:0,2:0,3:0};x=10;M=7**4;P=order(10,M)
    for m in range(P):
        z=(100*x*x+666*x+1)%M;v=0
        while z and z%7==0:v+=1;z//=7
        need(z!=0,'no fourth7-adic layer');expected=0 if m%6!=4 else 3 if m%42==40 else 2
        need(v==expected,'exact Belphegor ramified layers');valuation_counts[v]+=1;x=x*10%M
    out=dict(status='FINITE_EXACT_CONTROLS_FOR_PROVED_CLOCK_AND_NORM_LAWS',clock_cases=clocks,bases=[2,12],moduli=[2,300],
        free_weight_smith_decimal=[1,9,89100],belphegor_smith=[1,18,296703000],
        prime7417=dict(order=T,wing_zero_phases=hits,belphegor_zero_phases=bhits,discriminant_legendre=pow(443156%p,(p-1)//2,p),supplied_record_phase=944264%T),
        wing_nonresidue_primes=nonres,belphegor_mod2401=dict(period=P,valuation_counts=valuation_counts),
        residue_ring_clocks={str(m):order(10,m) for m in (49,343,1001,7417)},
        scope='No primality claim for a whole family; no re-proving the published 1888529-digit certificate; residue-ring clocks are not finite-field Frobenius')
    (HERE/'arithmetic.json').write_text(json.dumps(out,indent=2)+'\n');print('PASS all-modulus clock controls, Smith forms and ramified norm charts',flush=True)
if __name__=='__main__':main()
