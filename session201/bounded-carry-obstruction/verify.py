"""Exact controls for bounded source-error sparsity; no target coverage scan."""
from fractions import Fraction as F
from math import comb
import json

def need(c,m):
    if not c:raise RuntimeError(m)

def error_runs(k,A):
    t=F(1,2);e=F(0);i=0;runs=[]
    while i<k:
        o=0
        while i<k and ((A>>i)&1)==0:t/=2;i+=1;o+=1
        if i==k:break
        start=t;r=0
        need(F(1,1<<(o+1))<=start<=F(1,2),'run-start lower and upper bounds')
        while i<k and ((A>>i)&1)==1:t=(t+1)/2;i+=1;r+=1
        eta=start+(1-start)/(1<<r);e+=eta
        lower=F(1,1<<(o+1))+F(1,1<<(r+1))
        need(eta>=lower,'run additive lower bound')
        need(eta*eta*(1<<(o+r))>=1,'run geometric-mean lower bound')
        runs.append((o,r,eta))
    R=len(runs)
    if R:
        # Exact rational form of e >= R*2^(-k/(2R)); no floating powers.
        need(e.numerator**(2*R)*(1<<k)>=R**(2*R)*e.denominator**(2*R),'Jensen lower bound')
    else:need(e==0,'zero runs exact zero error')
    return e,runs

def allowed(E,k,R):
    if R==0:return True
    return E.numerator**(2*R)*(1<<k)>=R**(2*R)*E.denominator**(2*R)

def main():
    thresholds=(F(0),F(1,2),F(1),F(2));details=[];words=0
    for k in range(15):
        run_counts=[0]*((k+1)//2+1);counts={E:0 for E in thresholds}
        for A in range(1<<k):
            e,runs=error_runs(k,A);R=len(runs);run_counts[R]+=1
            need(sum(o+r for o,r,_ in runs)<=k,'run lengths do not exceed word length')
            for E in thresholds:
                if e<=E:counts[E]+=1;need(allowed(E,k,R),'necessary run count condition')
            words+=1
        need(run_counts==[comb(k+1,2*R) for R in range(len(run_counts))],'exact word enumeration by Z runs')
        for E in thresholds:
            bound=sum(comb(k+1,2*R) for R in range(len(run_counts)) if allowed(E,k,R))
            need(counts[E]<=bound,'finite combinatorial upper bound')
            details.append({'k':k,'E':str(E),'actual_words':counts[E],'necessary_run_bound':bound})
    for m in range(1,81):
        k=2*m;A=((1<<(2*m))-1)//3;e,runs=error_runs(k,A)
        need(len(runs)==m and all(o==(0 if i==0 else 1) and r==1 for i,(o,r,_) in enumerate(runs)),'alternating-word runs')
        need(e==F(2*m,3)+F(1,9)*(1-F(1,1<<(2*m))),'explicit linear-growth hostile')
    return {'status':'PASS','source_words_checked':words,'thresholds':[str(E) for E in thresholds],
            'exact_bounded_error_counts':details,'alternating_word_controls_through_m':80,
            'all_order_statement':'For fixed E, bounded-error source indices through N are N^{o(1)}',
            'scope':'Both source errors bounded; no claim that either error must grow for every target',
            'universal_occupancy':'OPEN'}

if __name__=='__main__':print(json.dumps(main(),sort_keys=True,indent=2))
