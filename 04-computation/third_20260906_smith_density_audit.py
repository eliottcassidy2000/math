#!/usr/bin/env python3
"""Independent exact unit-density and full-kernel audit of Smith Case III."""
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from math import gcd
import json

GATES=0
TRACE=sha256()

def gate(ok,label,data=None):
    global GATES
    if not ok:raise RuntimeError(f'{label}: {data!r}')
    GATES+=1
    TRACE.update((label+':'+repr(data)+'\n').encode())

def valuation(n,p):
    if n==0:return None
    n=abs(n);v=0
    while n%p==0:n//=p;v+=1
    return v

def capped(n,p,K):
    v=valuation(n,p)
    return K if v is None else min(K,v)

def expected(K,p,j):
    if K==0:return F(j==0)
    if j==0:return F(p-2,p-1)
    if j==K:return F(1,(p-1)*p**(K-1))
    return F(1,p**j)

def reference_checks(m,p,d,e,K,units):
    modulus=p**e;c=m//p**d
    samples=sorted({units[0],units[-1],units[len(units)//2],c*pow(m+1,-1,p**K)%(p**K)})
    for u in samples:
        for v in (1,m+1):
            a=p**(e+d)*u;b=p**e*v
            gate(valuation(a,p)==e+d and valuation(b,p)==valuation(b-a,p)==e,'fixed-weighted-metric',(m,p,e,u,v))
            tau=u*pow(v,-1,modulus)%modulus
            k=capped(c*v-(m+1)*u,p,K)
            for s in (-1,0,1,p-1,p+1):
                # w=(1,s): tau_w=(u/v)*(1-bs)/(1-as).
                gate(gcd(1-a*s,p)==gcd(1-b*s,p)==1,'admissible-unit-reference')
                tw=tau*(1-b*s)*pow(1-a*s,-1,modulus)%modulus
                gate(tw==tau and capped(c-(m+1)*tw,p,K)==k,'reference-change-preserves-capped-coordinate',(m,p,e,u,v,s))
            for lift in (-2,-1,0,1,3):
                gate(capped(c*v-(m+1)*(u+lift*p**K),p,K)==k,'unit-residue-lift-invariance')
            # Independently recover the largest Smith exponent from three
            # actual determinantal valuations, retaining both competing caps.
            L=m*b-(m+1)*a
            rho=min(m*valuation(a,p),valuation(m,p)+(m-1)*valuation(a,p),m*valuation(b,p))
            minor_vals=[2*m*valuation(a,p),m*valuation(a,p)+m*valuation(b,p)+valuation(b-a,p)]
            if L:minor_vals.append((m-1)*valuation(a,p)+m*valuation(b,p)+valuation(L,p))
            sigma=min(minor_vals)
            total=2*m*valuation(a,p)+m*valuation(b,p)+2*valuation(b-a,p)
            listed=((m-1)*e+K,(m+1)*e+m*d-K+k,(m+2)*e+m*d-k)
            gate((rho,sigma-rho,total-sigma)==listed,'literal-three-determinantal-valuations',(m,p,e,u,v))

def kernel_checks(m,p,d,e,K,hist,total):
    P=(m+1)*e+m*d-K;Q=(m+2)*e+m*d;fixed=(m-1)*e+K
    gate(Q-P==e+K>=2*K,'ordered-factor-room',(m,p,e,K))
    for N in range(1,Q+K+3):
        R=max(0,min(K,N-P,Q-N))
        base=sum(min(N,s) for s in (fixed,P,Q))
        exponent_differences=[]
        for j in range(K+1):
            spectrum=(fixed,P+j,Q-j)
            gate(spectrum==tuple(sorted(spectrum)),'all-factor-states-ordered',(m,p,e,j))
            full=sum(min(N,s) for s in spectrum)
            change=full-base
            gate(change==min(j,R),'full-kernel-exponent-difference',(m,p,e,N,j))
            exponent_differences.append(change)
        logmean=sum((F(hist[j],total)*exponent_differences[j] for j in range(K+1)),F(0))
        cardinalmean=sum((F(hist[j],total)*p**exponent_differences[j] for j in range(K+1)),F(0))
        gate(logmean==F(p,(p-1)**2)*(1-F(1,p**R)),'expected-log-kernel-ratio',(m,p,e,N))
        gate(cardinalmean==R+1,'expected-full-kernel-cardinality-ratio',(m,p,e,N))
        if p==2 and K:
            actualbase=exponent_differences[1]
            actualmean=sum((F(hist[j],total)*F(p)**(exponent_differences[j]-actualbase) for j in range(1,K+1)),F(0))
            gate(actualmean==(F(R+1,2) if R else F(1)),'attained-dyadic-baseline-kernel-ratio',(m,p,e,N))
    central=P+K
    gate(central<=Q-K and max(0,min(K,central-P,Q-central))==K,'nonempty-central-precision-band',(m,p,e,K))

def main():
    records=[];unit_cases=pair_cases=models=zero_models=0
    for p in (2,3,5,7):
        for m in range(2,21):
            d=valuation(m,p)
            if not d:continue
            # Every feasible cap<=5 is present. Above e=md the law repeats;
            # two extra e values challenge the same cap with different metrics.
            es=set(range(0,min(5,m*d)+1))
            if m*d<=5:es|={m*d+1,m*d+7}
            for e in sorted(es):
                K=min(e,m*d);models+=1
                if K==0:
                    hist=Counter({0:1});total=1
                    gate(expected(0,p,0)==1,'zero-cap-deterministic',(m,p))
                    kernel_checks(m,p,d,e,K,hist,total);zero_models+=1
                    continue
                modulus=p**K;units=[u for u in range(modulus) if u%p]
                c=m//p**d
                hist=Counter(capped(c-(m+1)*tau,p,K) for tau in units)
                total=len(units);unit_cases+=total
                gate(total==(p-1)*p**(K-1),'full-unit-universe-size',(m,p,e,K))
                for j in range(K+1):gate(F(hist[j],total)==expected(K,p,j),'exact-density',(m,p,e,K,j))
                for r in range(1,K+1):
                    tail=sum(n for j,n in hist.items() if j>=r)
                    gate(F(tail,total)==F(1,(p-1)*p**(r-1)),'exact-large-gain-tail',(m,p,e,r))
                mean=sum((F(n*j,total) for j,n in hist.items()),F(0))
                gate(mean==F(p,(p-1)**2)*(1-F(1,p**K))<2,'bounded-mean-precision-improvement',(m,p,e,K))
                if p==2:
                    gate(hist[0]==0 and mean-1==1-F(1,2**(K-1))<1,'dyadic-attained-baseline',(m,p,e,K))
                if modulus<=49:
                    pairs=Counter();ratios=Counter()
                    for v in units:
                        inv=pow(v,-1,modulus)
                        for u in units:
                            tau=u*inv%modulus
                            ratios[tau]+=1;pairs[capped(c*v-(m+1)*u,p,K)]+=1
                    gate(set(ratios)==set(units) and set(ratios.values())=={total},'independent-unit-pair-ratio-fibers',(m,p,e,K))
                    gate(all(pairs[j]==hist[j]*total for j in range(K+1)),'full-small-unit-pair-law',(m,p,e,K))
                    pair_cases+=total**2
                reference_checks(m,p,d,e,K,units)
                kernel_checks(m,p,d,e,K,hist,total)
                records.append({'p':p,'m':m,'d':d,'e':e,'K':K,'unit_count':total,'counts':[hist[j] for j in range(K+1)],'mean_kappa':str(mean)})
    print('FINITE_MODELS '+json.dumps(records,sort_keys=True,separators=(',',':')))
    print('UNIVERSE primes=2,3,5,7 m=2..20 divisible_by_p all_feasible_K=1..5 e=K_plus_saturated_extra_metrics; models='+str(models)+' zero_cap_models='+str(zero_models)+' full_unit_ratios='+str(unit_cases)+' small_full_unit_pairs='+str(pair_cases))
    print('CLAIM exact_uniform_unit_density_and_tail; mean_kappa<2; dyadic_optional_mean<1; expected_kernel_cardinality_ratio=R+1; attained_dyadic_ratio=(R+1)/2_when_Rpositive_else1')
    print('SCOPE fixed_CaseIII_weighted_metric_and_declared_uniform_unit_ratio; no_distribution_on_arbitrary_projective_observers')
    print('PASS explicit_gates='+str(GATES)+' semantic_sha256='+TRACE.hexdigest())

if __name__=='__main__':main()
