"""Integer partial-permutation weights for the inherited projective Burnside DP."""
from math import factorial,gcd
from collections import defaultdict
from fractions import Fraction
import importlib.util,sys,time,json
from support import HERE,ROOT,need,save
sys.set_int_max_str_digits(0)

def census(N,integral=True):
    lengths=list(range(N if N%2 else N-1,0,-2));states={(N,(0,)*len(lengths)):1 if integral else Fraction(1)}
    facts=[factorial(i) for i in range(N+1)];counts=[];negative=0;transitions=0
    for stage,p in enumerate(lengths):
        if p==1:
            total=sum(w*(1<<(r*(r-1)//2))/(1 if integral else facts[r]) for (r,s),w in states.items()) if not integral else sum(w<<(r*(r-1)//2) for (r,s),w in states.items())
            if integral:
                result,rem=divmod(total,facts[N]);need(rem==0,'final Burnside divisibility')
            else:need(total.denominator==1,'rational integer output');result=total.numerator
            return result,counts,negative,transitions
        remaining=lengths[stage+1:];nxt=defaultdict(int if integral else Fraction)
        for (r,s),weight in states.items():
            active=[j for j in lengths[stage:] if j<=r];S=dict(zip(active,s))
            for k in range(r//p+1):
                rr=r-k*p;exponent=k*(p-1)//2+p*k*(k-1)//2+k*S.get(p,0)+rr*k
                ss=tuple(S[j]+k*(gcd(p,j)-j) for j in remaining if j<=rr)
                if integral:
                    multiplicity,rem=divmod(facts[r],facts[rr]*p**k*facts[k]);need(rem==0,'partial cycle placement integer')
                    z=weight*multiplicity
                    if exponent>=0:z<<=exponent
                    else:
                        z,rem=divmod(z,1<<(-exponent));need(rem==0,'exact negative gauge shift');negative+=1
                else:
                    z=weight*Fraction(1<<exponent,1) if exponent>=0 else weight/Fraction(1<<(-exponent),1)
                    z/=p**k*facts[k]
                nxt[(rr,ss)]+=z;transitions+=1
        states=nxt;counts.append(len(states))
    return 1,counts,negative,transitions

def reference():
    path=ROOT/'05-knowledge/results/session98_tournament_growth_20261029/recurrence/verify.py'
    spec=importlib.util.spec_from_file_location('prior_census',path);mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod);return mod

def main():
    old=reference();small=[]
    for n in range(31):
        v,st,neg,tr=census(n);need(v==old.direct(n),'direct odd-partition control');small.append(v)
    need(small[14]==28401423719122304 and small[15]==31021002160355166848,'float-rounding hostiles')
    big=[];timings=[]
    targets=[60,100,150,200] if '--quick' not in sys.argv else [60,100]
    for n in targets:
        start=time.perf_counter();a,counts,neg,tr=census(n);ti=time.perf_counter()-start
        start=time.perf_counter();b,cc,_,_=census(n,False);tf=time.perf_counter()-start
        need(a==b and counts==cc,'same states independent rational scaling')
        if n<=100:
            c,oc=old.fixed_degree(n);need(a==c and counts==oc[:-1],'unchanged inherited projective carrier')
        big.append({'n':n,'value':str(a),'peak_states':max(counts),'stage_states':counts,'transitions':tr,'negative_shift_checks':neg})
        timings.append({'n':n,'integer_seconds':ti,'fraction_seconds':tf,'ratio_fraction_over_integer':tf/ti})
        print('census',n,'integer',round(ti,3),'fraction',round(tf,3),'negative shifts',neg,flush=True)
    save('integer_census.json',{'status':'EXACT_INTEGER_RECURRENCE_AND_FINITE_CONTROLS','small_values':small,'large_values':big,'scope':'same inherited fixed-degree projective response; new integer weights and direct unit completion; no new asymptotic complexity claim'})
    save('census_timings.json',{'status':'OBSERVED_WALL_CLOCK_NOT_A_THEOREM','runs':timings})
    print('PASS integer census',flush=True)
if __name__=='__main__':main()
