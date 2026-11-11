"""Common Gaussian mode stripped from inherited odd-factorial recurrence."""
from collections import defaultdict
from fractions import Fraction
from math import factorial,gcd

def common_dp(bound):
    width=max(1,bound.bit_length()); digit=(1<<width)-1
    masks=[(1<<(width*c))-1 for c in range(bound//2+2)]
    facts=[factorial(i) for i in range(bound+1)]
    twos=[(f&-f).bit_length()-1 for f in facts]
    odds=[f>>v for f,v in zip(facts,twos)]
    states={(0,0):1}; rows=[]
    for ell in range(bound if bound%2 else bound-1,0,-2):
        idx=(ell-1)//2
        row=sum((1 if j==1 else gcd(ell,j)-1)<<(width*i) for i,j in enumerate(range(1,ell,2)))
        tables={}
        for degree,_ in states:
            if degree in tables:continue
            entries=[]
            for k in range((bound-degree)//ell+1):
                nd=degree+ell*k
                factor,rem=divmod(odds[nd],odds[degree]*ell**k*odds[k]);assert rem==0
                count=max(0,(min(ell-2,bound-nd)-1)//2+1)
                entries.append((k,nd,(ell-1)*k*k//2-twos[k],factor,masks[count],k*row))
            tables[degree]=entries
        new=defaultdict(int);trans=0
        for (degree,response),weight in states.items():
            K=response&digit
            charge=0 if ell==1 else (response>>(width*idx))&digit
            for k,nd,base,factor,mask,offset in tables[degree]:
                e=base+k*charge
                if not mask:e+=(K+k)*(K+k-1)//2
                assert e>=0
                new[nd,(response+offset)&mask]+=(weight*factor)<<e;trans+=1
        states=dict(new);rows.append(dict(cycle_length=ell,states=len(states),transitions=trans))
    ans=[Fraction(0)]*(bound+1)
    for (d,r),w in states.items():assert not r;ans[d]+=Fraction(w,odds[d])
    return ans,dict(rows=rows,transitions=sum(r['transitions'] for r in rows))

if __name__=='__main__':
    import importlib.util,time,json,pathlib
    root=pathlib.Path(__file__).resolve().parents[4]
    p=root/'05-knowledge/results/degree_local_tournament_enumeration_20261021/session22-enumeration/local_scale.py'
    spec=importlib.util.spec_from_file_location('baseline',p);b=importlib.util.module_from_spec(spec);spec.loader.exec_module(b)
    for n in (0,1,2,3,4,5,10,20,30,50,100,150,200):
        t=time.perf_counter();a,info=common_dp(n);gt=time.perf_counter()-t
        t=time.perf_counter();bb,bi=b.local_dp(n);bt=time.perf_counter()-t
        assert a==bb
        print(json.dumps(dict(N=n,gauge_seconds=gt,baseline_seconds=bt,speedup=bt/gt,exact_equal=True,transitions=info['transitions'])),flush=True)
